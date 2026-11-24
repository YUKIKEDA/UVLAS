#include "case_dictionary.hpp"

#include "validate_case.hpp"

#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string_view>
#include <utility>

namespace uvlas::python {
namespace {

namespace py = pybind11;
using boost::multiprecision::cpp_int;

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
        int overflow = 0;
        const auto small = PyLong_AsLongLongAndOverflow(object, &overflow);
        if (overflow == 0) {
            if (PyErr_Occurred() != nullptr) {
                throw py::error_already_set{};
            }
            append(std::to_string(small));
            return;
        }
        // Use the pinned CPython 3.12 API to copy the integer without calling
        // Python or relying on the interpreter's decimal conversion setting.
        constexpr auto bit_limit = 4 * app::case_manifest_max_number_bytes;
        const auto bits = _PyLong_NumBits(object);
        if (bits > bit_limit) {
            if (PyErr_Occurred() != nullptr) {
                PyErr_Clear();
            }
            exhausted();
        }
        std::array<unsigned char, bit_limit / 8 + 1> bytes{};
        const auto size = bits / 8 + 1;
        if (_PyLong_AsByteArray(reinterpret_cast<PyLongObject*>(object), bytes.data(), size, 1,
                                1) != 0) {
            throw py::error_already_set{};
        }
        cpp_int integer;
        import_bits(integer, bytes.data(), bytes.data() + size, 8, false);
        if (overflow < 0) {
            integer -= cpp_int{1} << (8 * size);
        }
        const auto digits = integer.str();
        if (digits.size() > app::case_manifest_max_number_bytes) {
            exhausted();
        }
        append(digits);
    }

    void decimal(double number) {
        if (!std::isfinite(number)) {
            invalid();
        }
        static_assert(std::numeric_limits<double>::is_iec559 &&
                      std::numeric_limits<double>::digits == 53);
        const auto bits = std::bit_cast<std::uint64_t>(number);
        if (std::signbit(number)) {
            append("-");
        }
        constexpr int fraction_bits = 52;
        constexpr int exponent_bias = 1023;
        constexpr std::uint64_t hidden = std::uint64_t{1} << fraction_bits;
        const auto biased = static_cast<int>((bits >> fraction_bits) & 0x7ff);
        auto significand = (bits & (hidden - 1)) + (biased == 0 ? 0 : hidden);
        if (significand == 0) {
            append("0.0");
            return;
        }
        const auto trailing = std::countr_zero(significand);
        significand >>= trailing;
        const auto exponent = std::max(1, biased) - exponent_bias - fraction_bits + trailing;
        cpp_int coefficient{significand};
        if (exponent >= 0) {
            coefficient <<= static_cast<unsigned>(exponent);
            append(coefficient.str());
            append(".0");
        } else {
            // Replace 2^-n with 5^n * 10^-n so schema comparisons use the
            // exact binary64 value, including when compared with an integer.
            coefficient *= boost::multiprecision::pow(cpp_int{5}, static_cast<unsigned>(-exponent));
            append(coefficient.str());
            append("e");
            append(std::to_string(exponent));
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
