#include "uvlas/core/build_info.hpp"
#include "validate_case.hpp"

#include <filesystem>
#include <iostream>
#include <string_view>

namespace {

constexpr std::string_view usage = "Usage: uvlas --help | --version | --build-info\n"
                                   "       uvlas validate case.json\n";

int run(int argc, std::filesystem::path::value_type* argv[]) {
    if (argc < 2) {
        std::cerr << "Expected exactly one option.\n" << usage;
        return 2;
    }

    const auto argument = std::filesystem::path{argv[1]}.u8string();
    const std::string_view option{reinterpret_cast<const char*>(argument.data()), argument.size()};
    if (option == "validate") {
        if (argc != 3) {
            std::cerr << "Expected one case path.\n" << usage;
            return 2;
        }
        return uvlas::app::validate_case(std::filesystem::path{argv[2]});
    }
    if (argc != 2) {
        std::cerr << "Expected exactly one option.\n" << usage;
        return 2;
    }
    if (option == "--help") {
        std::cout << usage;
        return 0;
    }
    if (option == "--version") {
        const auto info = uvlas::build_info();
        std::cout << "UVLAS " << info.version << '\n'
                  << "Compiler: " << info.compiler << '\n'
                  << "Platform: " << info.platform << '\n'
                  << "Configuration: " << info.configuration << '\n'
                  << "Instruction set: " << info.instruction_set << '\n';
        return 0;
    }

    if (option == "--build-info") {
        std::cout << uvlas::build_info_json();
        return 0;
    }

    std::cerr << "Unknown option: " << option << '\n' << usage;
    return 2;
}

} // namespace

// No exception leaves the process: an unexpected one ends it with exit code 1.
#if defined(_WIN32)
int wmain(int argc, wchar_t* argv[]) {
#else
int main(int argc, char* argv[]) {
#endif
    try {
        return run(argc, argv);
    } catch (...) {
        return 1;
    }
}
