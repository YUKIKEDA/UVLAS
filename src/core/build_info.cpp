#include "uvlas/core/build_info.hpp"

#include "uvlas/build_config.hpp"
#include "uvlas/git_revision.hpp"

#include <Eigen/Core>
#include <Spectra/Util/Version.h>
#include <boost/version.hpp>
#include <nlohmann/json.hpp>

#include <format>

// The widest instruction set the compiler may use for this build. It is part
// of the platform key: raising it is a change of numerics (QA-07).
#if defined(__AVX512F__)
#define UVLAS_INSTRUCTION_SET "avx512"
#elif defined(__AVX2__)
#define UVLAS_INSTRUCTION_SET "avx2"
#elif defined(__AVX__)
#define UVLAS_INSTRUCTION_SET "avx"
#else
#define UVLAS_INSTRUCTION_SET "sse2"
#endif

namespace uvlas {

BuildInfo build_info() noexcept {
    return {.version = UVLAS_VERSION,
            .commit = UVLAS_GIT_COMMIT,
            .uncommitted_changes = UVLAS_GIT_DIRTY != 0,
            .compiler = UVLAS_COMPILER,
            .platform = UVLAS_PLATFORM,
            .configuration = UVLAS_BUILD_CONFIGURATION,
            .floating_point = UVLAS_FLOATING_POINT,
            .instruction_set = UVLAS_INSTRUCTION_SET,
            .platform_key =
                UVLAS_OS_KEY "-" UVLAS_ARCH_KEY "-" UVLAS_COMPILER_KEY "-" UVLAS_INSTRUCTION_SET};
}

std::string build_info_json() {
    const BuildInfo info = build_info();
    return std::format("{{\n"
                       "  \"schema\": \"uvlas.build-info/1\",\n"
                       "  \"version\": \"{}\",\n"
                       "  \"commit\": \"{}\",\n"
                       "  \"uncommitted_changes\": {},\n"
                       "  \"compiler\": \"{}\",\n"
                       "  \"platform\": \"{}\",\n"
                       "  \"configuration\": \"{}\",\n"
                       "  \"floating_point\": \"{}\",\n"
                       "  \"instruction_set\": \"{}\",\n"
                       "  \"platform_key\": \"{}\",\n"
                       "  \"dependencies\": {{\n"
                       "    \"eigen\": \"{}.{}.{}\",\n"
                       "    \"spectra\": \"{}.{}.{}\",\n"
                       "    \"boost\": \"{}.{}.{}\",\n"
                       "    \"nlohmann_json\": \"{}.{}.{}\"\n"
                       "  }}\n"
                       "}}\n",
                       info.version, info.commit, info.uncommitted_changes, info.compiler,
                       info.platform, info.configuration, info.floating_point, info.instruction_set,
                       info.platform_key, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION,
                       EIGEN_PATCH_VERSION, SPECTRA_MAJOR_VERSION, SPECTRA_MINOR_VERSION,
                       SPECTRA_PATCH_VERSION, BOOST_VERSION / 100000, BOOST_VERSION / 100 % 1000,
                       BOOST_VERSION % 100, NLOHMANN_JSON_VERSION_MAJOR,
                       NLOHMANN_JSON_VERSION_MINOR, NLOHMANN_JSON_VERSION_PATCH);
}

} // namespace uvlas
