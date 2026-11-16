// Product complete-spectrum and derived-result baseline, with no rounded decimal values.
#include <array>
#include <vector>
#include <cstdint>
#include <exception>
#include <format>
#include <iostream>
#include <stdexcept>
#include <uvlas/core/build_info.hpp>
#include <uvlas/numerics/eigen_full.hpp>
#include <variant>
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
            0.};
}
[[nodiscard]] eig::FullRequest request(core::NumericalGeneration gen) {
    return {{{1e-12, core::ObservationId{1}, core::ScopeId{1}},
             {1e-12, core::ObservationId{2}, core::ScopeId{1}},
             eig::NotRequested{}},
            numerics::linear::RetainedPrecision{gen},
            {3, 3, {32, 32, 0, 2}, numerics::linear::EstimateAdmission{std::uint64_t{1024} * 1024}},
            {},
            nullptr,
            false};
}
[[nodiscard]] double bound(const core::ErrorCertificate& certificate) {
    const auto* a = std::get_if<core::Achieved>(&certificate.assessment());
    if (!a)
        throw std::runtime_error("fixed result not qualified");
    return a->upper_bound;
}
void emit(const char* name, double value) {
    std::cout << std::format("{} {:a}\n", name, value);
}
void report(const char* name,int repetition,const eig::FullTimings& t,std::uint64_t payload) {
 std::cout<<name<<" "<<repetition<<" "<<t.preparation.elapsed_nanoseconds<<" "<<t.backend.elapsed_nanoseconds<<" "<<t.qualification.elapsed_nanoseconds<<" "<<payload<<'\n';
}
void evaluate(){
 for(int repetition=0;repetition<3;++repetition)for(const std::size_t n:{std::size_t{6},std::size_t{24}}){
 core::ResourcePlan plan{std::uint64_t{64}*1024*1024};std::vector<double>a(n*n),d(n,1.);std::vector<numerics::linear::UnitDimensions>u(n),au(n*n);
 for(std::size_t i=0;i<n;++i){a[i+i*n]=static_cast<double>(i+1);if(i+1<n)a[i+(i+1)*n]=.125;}
 auto matrix=numerics::DenseReal::copy({static_cast<core::Index>(n),static_cast<core::Index>(n),static_cast<core::Index>(n)},a,plan,core::Operation::ExtractEigen);if(!matrix)throw std::runtime_error("timing input");
 auto input=numerics::PrecisionBuffer<double>::binary64(std::move(*matrix),core::NumericalGeneration{5});if(!input)throw std::runtime_error("timing precision");
 auto req=request(core::NumericalGeneration{40});req.limits.maximum_order=24;req.limits.maximum_results=24;req.limits.actions={200,200,0,2};req.timing=eig::FullTimingMode::Enabled;
 const auto before=plan.allocated_bytes();auto full=eig::full_spectrum(*input,context(),{},{d,1,u,{},au},req,plan);if(!full)throw std::runtime_error("timing spectrum");const auto t=full->timings();if(!t)throw std::runtime_error("no timing");report(n==6?"full6":"full24",repetition,*t,plan.allocated_bytes()-before);
 if(n==6){const std::array<core::Index,2>indices{0,2};auto selection=eig::FullSelection::make(*full,indices,plan);if(!selection)throw std::runtime_error("timing selection");
 auto output=req.output();output.returned=numerics::linear::RetainedPrecision{core::NumericalGeneration{41}};const eig::FullPairRequest pr{output,{{1e-12,core::ObservationId{3},core::ScopeId{1}},{1e-12,core::ObservationId{4},core::ScopeId{1}},eig::NotRequested{}},{0,1e-14},true,false};
 const auto pb=plan.allocated_bytes();auto pairs=eig::full_pairs(*full,*selection,pr,plan);if(!pairs)throw std::runtime_error("timing pairs");const auto pt=pairs->timings();if(!pt)throw std::runtime_error("no pair timing");report("pairs6x2",repetition,*pt,plan.allocated_bytes()-pb);
 auto cluster_req=req;cluster_req.returned=numerics::linear::RetainedPrecision{core::NumericalGeneration{42}};const eig::FullClusterRequest cr{cluster_req,{0,1e-14},true};const auto cb=plan.allocated_bytes();auto cluster=eig::full_cluster(*full,*selection,cr,plan);if(!cluster)throw std::runtime_error("timing cluster");const auto ct=cluster->timings();if(!ct)throw std::runtime_error("no cluster timing");report("cluster6x2",repetition,*ct,plan.allocated_bytes()-cb);
 }
 }
}

int run() {
    try {
        std::cout << "platform " << uvlas::build_info().platform_key << '\n';
        evaluate();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
} // namespace
int main() {
    try {
        return run();
    } catch (...) {
        return 1;
    }
}
