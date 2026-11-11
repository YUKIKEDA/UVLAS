#include "uvlas/core/build_info.hpp"

#include "uvlas/build_config.hpp"

namespace uvlas {

BuildInfo build_info() noexcept {
    return {UVLAS_VERSION, UVLAS_COMPILER, UVLAS_PLATFORM, UVLAS_BUILD_CONFIGURATION};
}

} // namespace uvlas
