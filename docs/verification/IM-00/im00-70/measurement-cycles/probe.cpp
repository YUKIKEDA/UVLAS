// Windows-only measurement probe; no practical performance acceptance claim.
#include <array>
#include <bit>
#include <chrono>
#include <complex>
#include <cstdint>
#include <exception>
#include <format>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>
#include <uvlas/core/build_info.hpp>
#include <uvlas/numerics/eigen_arnoldi.hpp>
#include <vector>
#include <windows.h>

namespace {
using namespace uvlas;
namespace eig = numerics::eigen;
[[nodiscard]] eig::Context context() {
    return {{core::ModelFingerprint{}, core::ModelGeneration{1}, core::GeometryGeneration{2},
             core::ChartGeneration{3}, core::SparsityGeneration{4}, core::NumericalGeneration{5},
             core::ScaleGeneration{6}, core::PrecisionGeneration{7}, core::Precision::Binary64,
             core::WakePolicyGeneration{8}, core::TimeIntervalGeneration{9},
             core::BackendGeneration{10}},
            core::TrialId{11},
            core::ScopeId{12},
            core::ScopeId{13},
            0};
}
class Rotation final : public eig::ActionProvider<double> {
public:
    [[nodiscard]] core::Index order() const noexcept override {
        return 8;
    }
    [[nodiscard]] eig::ActionStamp stamp() const noexcept override {
        return {context(), core::NumericalGeneration{21}, {}};
    }
    [[nodiscard]] numerics::SupportedPrecision
    supported_precision(eig::ActionKind) const noexcept override {
        return numerics::SupportedPrecision::Binary64Only;
    }
    [[nodiscard]] core::Result<void> validate_frozen_state() const override {
        return {};
    }
    [[nodiscard]] core::Result<void> apply(eig::ActionKind kind, eig::WordView<double> x,
                                           std::span<double> high, std::span<double> low,
                                           std::span<double> error) const override {
        constexpr std::array<double, 8> diagonal{0, 0, .5, .25, .125, .0625, .03125, .015625};
        for (std::size_t i = 0; i < 8; ++i) {
            const auto source = i == 0 ? 1 : (i == 1 ? 0 : i);
            double coefficient = i == 0 ? -4 : (i == 1 ? 4 : diagonal[i]);
            if (i < 2 && kind == eig::ActionKind::Adjoint) {
                coefficient = -coefficient;
            }
            const auto product = numerics::two_product(coefficient, x.high()[source]);
            if (!product) {
                return std::unexpected(product.error());
            }
            high[i] = product->hi;
            low[i] = product->lo;
            error[i] = 0;
        }
        return {};
    }
};
class NoConsumer final : public eig::ArnoldiConsumerEvaluator {
public:
    [[nodiscard]] numerics::SupportedPrecision supported_precision() const noexcept override {
        return numerics::SupportedPrecision::Binary64Only;
    }
    [[nodiscard]] core::Result<core::ErrorCertificate>
    assess(const eig::ArnoldiCandidateView&, const eig::ArnoldiCandidateView&,
           const numerics::linear::ConsumerCheck&, core::CertificateContext) const override {
        return std::unexpected(
            core::Failure{core::FailureCode::InvalidInput, core::Operation::ExtractEigen});
    }
};

class Diagonal final : public eig::ActionProvider<double> {
public:
    core::Index order() const noexcept override {
        return 96;
    }
    eig::ActionStamp stamp() const noexcept override {
        return {context(), core::NumericalGeneration{22}, {}};
    }
    numerics::SupportedPrecision supported_precision(eig::ActionKind) const noexcept override {
        return numerics::SupportedPrecision::Binary64Only;
    }
    core::Result<void> validate_frozen_state() const override {
        return {};
    }
    core::Result<void> apply(eig::ActionKind, eig::WordView<double> x, std::span<double> hi,
                             std::span<double> lo, std::span<double> error) const override {
        for (std::size_t i = 0; i < 96; ++i) {
            const auto p = numerics::two_product(96 - static_cast<double>(i), x.high()[i]);
            if (!p)
                return std::unexpected(p.error());
            hi[i] = p->hi;
            lo[i] = p->lo;
            error[i] = 0;
        }
        return {};
    }
};

class ComplexDiagonal final : public eig::ActionProvider<std::complex<double>> {
public:
    core::Index order() const noexcept override {
        return 6;
    }
    eig::ActionStamp stamp() const noexcept override {
        return {context(), core::NumericalGeneration{23}, {}};
    }
    numerics::SupportedPrecision supported_precision(eig::ActionKind) const noexcept override {
        return numerics::SupportedPrecision::Binary64Only;
    }
    core::Result<void> validate_frozen_state() const override {
        return {};
    }
    core::Result<void> apply(eig::ActionKind kind, eig::WordView<std::complex<double>> x,
                             std::span<std::complex<double>> hi, std::span<std::complex<double>> lo,
                             std::span<double> error) const override {
        constexpr std::array<std::complex<double>, 6> d{
            {{7, 1}, {3, 0}, {2, 0}, {1, 0}, {.5, 0}, {.25, 0}}};
        for (std::size_t i = 0; i < 6; ++i) {
            const auto a = kind == eig::ActionKind::Adjoint ? std::conj(d[i]) : d[i];
            const auto value = x.high()[i];
            hi[i] = a * value;
            lo[i] = {};
            const double magnitude = (std::abs(a.real()) + std::abs(a.imag())) *
                                     (std::abs(value.real()) + std::abs(value.imag()));
            error[i] = std::nextafter(16 * std::numeric_limits<double>::epsilon() * magnitude +
                                          16 * std::numeric_limits<double>::denorm_min(),
                                      std::numeric_limits<double>::infinity());
        }
        return {};
    }
};
std::size_t compare_words(std::string_view name, const eig::ArnoldiSpectrum& result, bool cluster,
                          bool both) {
    std::vector<std::uint64_t> bits;
    const auto add = [&](std::span<const std::complex<double>> values) {
        for (auto value : values) {
            bits.push_back(std::bit_cast<std::uint64_t>(value.real()));
            bits.push_back(std::bit_cast<std::uint64_t>(value.imag()));
        }
    };
    const auto add_proof = [&](const core::ErrorCertificate& p) {
        const auto* achieved = std::get_if<core::Achieved>(&p.assessment());
        if (achieved == nullptr)
            throw std::runtime_error("unachieved probe proof");
        bits.push_back(std::bit_cast<std::uint64_t>(achieved->upper_bound));
    };
    const auto v = result.values();
    if (!v)
        throw std::runtime_error("missing probe values");
    add(v->eigenvalues.values);
    const auto add_view = [&](const auto& part) {
        if constexpr (requires { part.vectors; })
            add(part.vectors.values);
        else if constexpr (requires { part.basis; }) {
            add(part.basis.values);
            add(part.triangular.values);
            for (auto index : part.diagonal_columns)
                bits.push_back(static_cast<std::uint64_t>(index));
        }
    };
    std::visit(add_view, v->right);
    std::visit(add_view, v->left);
    for (auto side : {eig::ArnoldiSide::Right, eig::ArnoldiSide::Left}) {
        if (side == eig::ArnoldiSide::Left && !both)
            continue;
        if (cluster) {
            const auto m = result.subspace_metric(side);
            if (!m)
                throw std::runtime_error("probe basis metrics");
            add_proof(m->residual);
            add_proof(m->orthogonality);
        } else
            for (core::Index j = 0; j < static_cast<core::Index>(v->eigenvalues.values.size());
                 ++j) {
                const auto m = result.pair_metric(side, j);
                if (!m)
                    throw std::runtime_error("probe pair metrics");
                add_proof(m->residual);
                add_proof(m->normalization);
                bits.push_back(std::bit_cast<std::uint64_t>(m->vector_infinity_lower));
            }
    }
    for (const auto& s : result.sensitivities()) {
        bits.push_back(static_cast<std::uint64_t>(s.index()));
        if (const auto* estimate = std::get_if<eig::ArnoldiConditionEstimate>(&s))
            bits.push_back(std::bit_cast<std::uint64_t>(estimate->value));
    }
    static std::map<std::string, std::vector<std::uint64_t>> controls;
    const auto [position, inserted] = controls.try_emplace(std::string(name), bits);
    if (!inserted && position->second != bits)
        throw std::runtime_error("probe repeated numerical words differ");
    return bits.size();
}
void sample(std::string_view name, core::Index n, core::Index nev, core::Index ncv, bool cluster,
            bool both, unsigned round, unsigned repeat) {
    const Rotation rotation;
    const Diagonal diagonal;
    const eig::ActionProvider<double>& provider =
        n == 8 ? static_cast<const eig::ActionProvider<double>&>(rotation)
               : static_cast<const eig::ActionProvider<double>&>(diagonal);
    const NoConsumer evaluator;
    std::vector<double> seed(static_cast<std::size_t>(n)), d(static_cast<std::size_t>(n), 1);
    for (std::size_t i = 0; i < seed.size(); ++i)
        seed[i] = static_cast<double>(i + 1);
    const std::vector<numerics::linear::UnitDimensions> units(static_cast<std::size_t>(n));
    const auto initial = eig::WordView<double>::make(seed, {}, core::Precision::Binary64);
    if (!initial)
        throw std::runtime_error("seed");
    const eig::MetricRequirement residual{1e-12, core::ObservationId{31}, core::ScopeId{32}},
        norm{1e-12, core::ObservationId{33}, core::ScopeId{32}};
    eig::ArnoldiQualification qualification =
        eig::ArnoldiPairRequest{{residual, norm, eig::NotRequested{}}};
    if (cluster)
        qualification = eig::ArnoldiSubspaceRequest{{residual, norm, eig::NotRequested{}}};
    const eig::ArnoldiRequest<double> request{
        nev,
        ncv,
        *initial,
        1e-14,
        100,
        eig::ArnoldiLargestMagnitude{},
        qualification,
        both ? eig::ArnoldiSides::Both : eig::ArnoldiSides::RightOnly,
        {1e-12, 1e-12},
        numerics::linear::RetainedPrecision{core::NumericalGeneration{40}},
        {},
        {},
        evaluator,
        eig::ArnoldiTimingMode::Enabled};
    const eig::ArnoldiScales scales{{d, 1}, units, {}, context()};
    const eig::ArnoldiLimits limits{n,
                                    nev + 1,
                                    ncv,
                                    100,
                                    {10000, 10000, 0, both ? 2U : 1U},
                                    numerics::linear::EstimateAdmission{10000000}};
    core::ResourcePlan scratch(1000000), output(1000000);
    auto workspace = eig::ArnoldiWorkspace::create(n, nev + 1, scratch);
    if (!workspace)
        throw std::runtime_error("workspace");
    ULONG64 cycles_before = 0, cycles_after = 0;
    if (!QueryThreadCycleTime(GetCurrentThread(), &cycles_before))
        throw std::runtime_error("thread cycle start");
    const auto wall_start = std::chrono::steady_clock::now();
    auto result = eig::arnoldi_spectrum(provider, scales, request, limits, *workspace, output);
    const auto extraction_wall_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                        std::chrono::steady_clock::now() - wall_start)
                                        .count();
    if (!QueryThreadCycleTime(GetCurrentThread(), &cycles_after))
        throw std::runtime_error("thread cycle finish");
    const auto thread_cycles = cycles_after - cycles_before;
    if (!result)
        throw std::runtime_error(std::format("{} extraction failed {}", name,
                                             static_cast<unsigned>(result.error().cause.code)));
    auto measured = result->measurements();
    if (!measured)
        throw std::runtime_error("measurement");
    const auto timing = std::get<eig::ArnoldiTimings>(*measured);
    const auto memory = core::sample_process_memory(core::Operation::ExtractEigen);
    if (!memory)
        throw std::runtime_error("memory");
    const auto compared_words = compare_words(name, *result, cluster, both);
    if (repeat == 0)
        return;
    const auto values = result->values();
    if (!values)
        throw std::runtime_error("values");
    std::cout << "{\"case\":\"" << name << "\",\"compared_words\":" << compared_words
              << ",\"repeat_equal\":true,\"round\":" << round << ",\"repeat\":" << repeat
              << ",\"thread_cycles\":" << thread_cycles
              << ",\"extraction_wall_ns\":" << extraction_wall_ns
              << ",\"workspace_payload_bytes\":" << scratch.allocated_bytes()
              << ",\"result_payload_bytes\":" << output.allocated_bytes()
              << ",\"resident_bytes\":" << memory->resident_bytes
              << ",\"process_lifetime_peak_bytes\":" << memory->peak_resident_bytes
              << ",\"root_real_bits\":"
              << std::bit_cast<std::uint64_t>(values->eigenvalues.values[0].real())
              << ",\"root_imag_bits\":"
              << std::bit_cast<std::uint64_t>(values->eigenvalues.values[0].imag())
              << ",\"stages\":[";
    bool first = true;
    for (auto stage : {timing.preparation, timing.krylov, timing.qualification, timing.subspace}) {
        if (!first)
            std::cout << ',';
        first = false;
        std::cout << "{\"calls\":" << stage.calls
                  << ",\"nanoseconds\":" << stage.elapsed_nanoseconds << '}';
    }
    std::cout << "]}\n";
}

void sample_complex(unsigned round, unsigned repeat) {
    constexpr core::Index n = 6, nev = 2, ncv = 6;
    constexpr bool cluster = false, both = true;
    constexpr std::string_view name = "complex6-pairs-both";
    const ComplexDiagonal provider;
    const NoConsumer evaluator;
    std::vector<std::complex<double>> seed(static_cast<std::size_t>(n));
    const std::vector<double> d(static_cast<std::size_t>(n), 1);
    for (std::size_t i = 0; i < seed.size(); ++i)
        seed[i] = {static_cast<double>(i + 1), 1};
    const std::vector<numerics::linear::UnitDimensions> units(static_cast<std::size_t>(n));
    const auto initial =
        eig::WordView<std::complex<double>>::make(seed, {}, core::Precision::Binary64);
    if (!initial)
        throw std::runtime_error("seed");
    const eig::MetricRequirement residual{1e-12, core::ObservationId{31}, core::ScopeId{32}},
        norm{1e-12, core::ObservationId{33}, core::ScopeId{32}};
    eig::ArnoldiQualification qualification =
        eig::ArnoldiPairRequest{{residual, norm, eig::NotRequested{}}};
    if (cluster)
        qualification = eig::ArnoldiSubspaceRequest{{residual, norm, eig::NotRequested{}}};
    const eig::ArnoldiRequest<std::complex<double>> request{
        nev,
        ncv,
        *initial,
        1e-14,
        100,
        eig::ArnoldiLargestMagnitude{},
        qualification,
        both ? eig::ArnoldiSides::Both : eig::ArnoldiSides::RightOnly,
        {1e-12, 1e-12},
        numerics::linear::RetainedPrecision{core::NumericalGeneration{40}},
        {},
        {},
        evaluator,
        eig::ArnoldiTimingMode::Enabled};
    const eig::ArnoldiScales scales{{d, 1}, units, {}, context()};
    const eig::ArnoldiLimits limits{n,
                                    nev,
                                    ncv,
                                    100,
                                    {10000, 10000, 0, both ? 2U : 1U},
                                    numerics::linear::EstimateAdmission{10000000}};
    core::ResourcePlan scratch(1000000), output(1000000);
    auto workspace = eig::ArnoldiWorkspace::create(n, nev, scratch);
    if (!workspace)
        throw std::runtime_error("workspace");
    ULONG64 cycles_before = 0, cycles_after = 0;
    if (!QueryThreadCycleTime(GetCurrentThread(), &cycles_before))
        throw std::runtime_error("thread cycle start");
    const auto wall_start = std::chrono::steady_clock::now();
    auto result = eig::arnoldi_spectrum(provider, scales, request, limits, *workspace, output);
    const auto extraction_wall_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                        std::chrono::steady_clock::now() - wall_start)
                                        .count();
    if (!QueryThreadCycleTime(GetCurrentThread(), &cycles_after))
        throw std::runtime_error("thread cycle finish");
    const auto thread_cycles = cycles_after - cycles_before;
    if (!result)
        throw std::runtime_error(std::format("{} extraction failed {}", name,
                                             static_cast<unsigned>(result.error().cause.code)));
    auto measured = result->measurements();
    if (!measured)
        throw std::runtime_error("measurement");
    const auto timing = std::get<eig::ArnoldiTimings>(*measured);
    const auto memory = core::sample_process_memory(core::Operation::ExtractEigen);
    if (!memory)
        throw std::runtime_error("memory");
    const auto compared_words = compare_words(name, *result, cluster, both);
    if (repeat == 0)
        return;
    const auto values = result->values();
    if (!values)
        throw std::runtime_error("values");
    std::cout << "{\"case\":\"" << name << "\",\"compared_words\":" << compared_words
              << ",\"repeat_equal\":true,\"round\":" << round << ",\"repeat\":" << repeat
              << ",\"thread_cycles\":" << thread_cycles
              << ",\"extraction_wall_ns\":" << extraction_wall_ns
              << ",\"workspace_payload_bytes\":" << scratch.allocated_bytes()
              << ",\"result_payload_bytes\":" << output.allocated_bytes()
              << ",\"resident_bytes\":" << memory->resident_bytes
              << ",\"process_lifetime_peak_bytes\":" << memory->peak_resident_bytes
              << ",\"root_real_bits\":"
              << std::bit_cast<std::uint64_t>(values->eigenvalues.values[0].real())
              << ",\"root_imag_bits\":"
              << std::bit_cast<std::uint64_t>(values->eigenvalues.values[0].imag())
              << ",\"stages\":[";
    bool first = true;
    for (auto stage : {timing.preparation, timing.krylov, timing.qualification, timing.subspace}) {
        if (!first)
            std::cout << ',';
        first = false;
        std::cout << "{\"calls\":" << stage.calls
                  << ",\"nanoseconds\":" << stage.elapsed_nanoseconds << '}';
    }
    std::cout << "]}\n";
}

} // namespace
int main() {
    try {
        std::cout << uvlas::build_info_json() << '\n';
        DWORD_PTR process_mask = 0, system_mask = 0;
        if (!GetProcessAffinityMask(GetCurrentProcess(), &process_mask, &system_mask))
            throw std::runtime_error("cannot inspect probe affinity");
        std::cout << "{\"process_affinity_mask\":" << process_mask
                  << ",\"system_affinity_mask\":" << system_mask << "}\n";
        for (unsigned round = 0; round < 3; ++round) {
            for (unsigned repeat = 0; repeat <= 30; ++repeat)
                sample("real8-pairs-both", 8, 2, 8, false, true, round, repeat);
            for (unsigned repeat = 0; repeat <= 3; ++repeat)
                sample("real8-cluster-both", 8, 2, 8, true, true, round, repeat);
            for (unsigned repeat = 0; repeat <= 1; ++repeat)
                sample("real96-pairs", 96, 24, 64, false, false, round, repeat);
            for (unsigned repeat = 0; repeat <= 3; ++repeat)
                sample_complex(round, repeat);
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
