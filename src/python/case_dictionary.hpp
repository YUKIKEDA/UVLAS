#pragma once

#include <pybind11/pybind11.h>

#include <string>
#include <uvlas/core/failure.hpp>

namespace uvlas::python {

[[nodiscard]] core::Result<std::string> copy_case_dictionary(pybind11::handle);

} // namespace uvlas::python
