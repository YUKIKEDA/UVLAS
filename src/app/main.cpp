#include "uvlas/core/build_info.hpp"

#include <iostream>
#include <string_view>

namespace {

constexpr std::string_view usage = "Usage: uvlas --help | --version | --build-info\n"
                                   "Development build: solver commands are not implemented.\n";

int run(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Expected exactly one option.\n" << usage;
        return 2;
    }

    const std::string_view option{argv[1]};
    if (option == "--help") {
        std::cout << usage;
        return 0;
    }
    if (option == "--version") {
        const auto info = uvlas::build_info();
        std::cout << "UVLAS " << info.version << '\n'
                  << "Compiler: " << info.compiler << '\n'
                  << "Platform: " << info.platform << '\n'
                  << "Configuration: " << info.configuration << '\n';
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
int main(int argc, char* argv[]) {
    try {
        return run(argc, argv);
    } catch (...) {
        return 1;
    }
}
