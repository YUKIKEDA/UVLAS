#include "case_dictionary.hpp"

#include "validate_case.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <string_view>
#include <utility>

namespace uvlas::python {
namespace {

namespace py = pybind11;

struct EncodingFailure {
    core::FailureCode code;
};

[[noreturn]] void invalid() {
    throw EncodingFailure{core::FailureCode::InvalidInput};
}

[[noreturn]] void exhausted() {
    throw EncodingFailure{core::FailureCode::ResourceExceeded};
}

class ManifestEncoder {
public:
    std::string encode(PyObject* dictionary) {
        if (!PyDict_CheckExact(dictionary)) {
            invalid();
        }
        value(dictionary, 0);
        return std::move(bytes_);
    }

private:
    void append(std::string_view text) {
        if (text.size() > app::case_manifest_max_bytes - bytes_.size()) {
            exhausted();
        }
        bytes_.append(text);
    }

    void string(PyObject* object) {
        if (!PyUnicode_CheckExact(object)) {
            invalid();
        }
        if (static_cast<std::size_t>(PyUnicode_GET_LENGTH(object)) >
            app::case_manifest_max_string_bytes) {
            exhausted();
        }
        Py_ssize_t size = 0;
        const char* text = PyUnicode_AsUTF8AndSize(object, &size);
        if (text == nullptr) {
            if (!PyErr_ExceptionMatches(PyExc_UnicodeError)) {
                throw py::error_already_set{};
            }
            PyErr_Clear();
            invalid();
        }
        if (static_cast<std::size_t>(size) > app::case_manifest_max_string_bytes) {
            exhausted();
        }
        append("\"");
        constexpr std::string_view hex = "0123456789abcdef";
        const std::string_view utf8{text, static_cast<std::size_t>(size)};
        std::size_t start = 0;
        for (std::size_t i = 0; i < utf8.size(); ++i) {
            const auto ch = static_cast<unsigned char>(utf8[i]);
            if (ch >= 0x20 && ch != '"' && ch != '\\') {
                continue;
            }
            append(utf8.substr(start, i - start));
            if (ch == '"' || ch == '\\') {
                const std::array escaped{'\\', static_cast<char>(ch)};
                append({escaped.data(), escaped.size()});
            } else {
                const std::array escaped{'\\', 'u', '0', '0', hex[ch >> 4], hex[ch & 15U]};
                append({escaped.data(), escaped.size()});
            }
            start = i + 1;
        }
        append(utf8.substr(start));
        append("\"");
    }

    void integer(PyObject* object) {
        // Use the pinned CPython 3.12 API to avoid calling Python while
        // borrowing dictionary entries. A decimal digit needs fewer than four bits.
        const auto bits = _PyLong_NumBits(object);
        if (bits > 4 * app::case_manifest_max_number_bytes) {
            if (PyErr_Occurred() != nullptr) {
                PyErr_Clear();
            }
            exhausted();
        }
        const auto text = py::reinterpret_steal<py::str>(PyObject_Str(object));
        if (!text) {
            if (!PyErr_ExceptionMatches(PyExc_ValueError)) {
                throw py::error_already_set{};
            }
            PyErr_Clear();
            exhausted();
        }
        Py_ssize_t size = 0;
        const char* data = PyUnicode_AsUTF8AndSize(text.ptr(), &size);
        if (data == nullptr) {
            throw py::error_already_set{};
        }
        const std::string_view digits{data, static_cast<std::size_t>(size)};
        if (digits.size() > app::case_manifest_max_number_bytes) {
            exhausted();
        }
        append(digits);
    }

    void decimal(double number) {
        if (!std::isfinite(number)) {
            invalid();
        }
        std::array<char, 32> buffer{};
        const auto converted =
            std::to_chars(buffer.data(), buffer.data() + buffer.size(), number,
                          std::chars_format::general);
        if (converted.ec != std::errc{}) {
            invalid();
        }
        const std::string_view digits{buffer.data(),
                                      static_cast<std::size_t>(converted.ptr - buffer.data())};
        append(digits);
        if (digits.find_first_of(".eE") == std::string_view::npos) {
            append(".0");
        }
    }

    void container(PyObject* object, std::size_t depth, bool dictionary) {
        if (depth == ancestors_.size()) {
            exhausted();
        }
        const auto parent = ancestors_.begin() + static_cast<std::ptrdiff_t>(depth);
        if (std::find(ancestors_.begin(), parent, object) != parent) {
            invalid();
        }
        *parent = object;
        if (dictionary) {
            append("{");
            Py_ssize_t position = 0;
            PyObject* key = nullptr;
            PyObject* item = nullptr;
            bool first = true;
            while (PyDict_Next(object, &position, &key, &item) != 0) {
                if (!std::exchange(first, false)) {
                    append(",");
                }
                string(key);
                append(":");
                value(item, depth + 1);
            }
            append("}");
        } else {
            append("[");
            const auto size = PyList_GET_SIZE(object);
            for (Py_ssize_t i = 0; i < size; ++i) {
                if (i != 0) {
                    append(",");
                }
                value(PyList_GET_ITEM(object, i), depth + 1);
            }
            append("]");
        }
    }

    void value(PyObject* object, std::size_t depth) {
        if (++nodes_ > app::case_manifest_max_nodes) {
            exhausted();
        }
        if (PyDict_CheckExact(object)) {
            container(object, depth, true);
        } else if (PyList_CheckExact(object)) {
            container(object, depth, false);
        } else if (PyUnicode_CheckExact(object)) {
            string(object);
        } else if (object == Py_None) {
            append("null");
        } else if (PyBool_Check(object)) {
            append(object == Py_True ? "true" : "false");
        } else if (PyLong_CheckExact(object)) {
            integer(object);
        } else if (PyFloat_CheckExact(object)) {
            decimal(PyFloat_AS_DOUBLE(object));
        } else {
            invalid();
        }
    }

    std::array<PyObject*, app::case_manifest_max_depth> ancestors_{};
    std::string bytes_;
    std::size_t nodes_ = 0;
};

} // namespace

core::Result<std::string> copy_case_dictionary(py::handle dictionary) {
    try {
        return ManifestEncoder{}.encode(dictionary.ptr());
    } catch (const EncodingFailure& failure) {
        return std::unexpected(core::Failure{failure.code, core::Operation::PrepareModel});
    }
}

} // namespace uvlas::python
