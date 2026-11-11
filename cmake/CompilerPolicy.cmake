# SPDX-License-Identifier: Apache-2.0

if(MSVC)
    if(MSVC_VERSION LESS 1944)
        message(FATAL_ERROR "UVLAS requires MSVC 19.44 or newer with C++23 support.")
    endif()
    # CMake 3.30 and 4.1 map cxx_std_23 to /std:c++latest. Keep this
    # project on the explicit C++23 mode, without enabling later drafts.
    set(CMAKE_CXX23_STANDARD_COMPILE_OPTION "/std:c++23preview")
    set(CMAKE_CXX23_EXTENSION_COMPILE_OPTION "/std:c++23preview")
endif()

include(CheckCXXSourceCompiles)
check_cxx_source_compiles([[
    #include <expected>
    #if !defined(__cpp_lib_expected) || __cpp_lib_expected < 202211L
    #error std::expected with monadic operations is required
    #endif
    constexpr auto result = std::expected<int, int>{23}.transform(
        [](int value) { return value + 1; });
    static_assert(result.has_value() && result.value() == 24);
    int main() { return 0; }
]] UVLAS_HAS_REQUIRED_EXPECTED)

if(NOT UVLAS_HAS_REQUIRED_EXPECTED)
    message(FATAL_ERROR
        "UVLAS requires C++23 std::expected with monadic operations "
        "(__cpp_lib_expected >= 202211L). Select a supported compiler and "
        "standard library; GCC 11 is not supported.")
endif()
