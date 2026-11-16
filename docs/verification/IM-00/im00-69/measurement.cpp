#include <uvlas/numerics/eigen_generalized.hpp>
#include <uvlas/core/build_info.hpp>
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace uvlas;
namespace eig=numerics::eigen;
eig::Context context() {
 return {{core::ModelFingerprint{},core::ModelGeneration{1},core::GeometryGeneration{2},core::ChartGeneration{3},core::SparsityGeneration{4},core::NumericalGeneration{5},core::ScaleGeneration{6},core::PrecisionGeneration{7},core::Precision::Binary64,core::WakePolicyGeneration{8},core::TimeIntervalGeneration{9},core::BackendGeneration{10}},core::TrialId{11},core::ScopeId{12},core::ScopeId{13},0.};
}
void sample(core::Index n, unsigned repetitions, unsigned round) {
 std::vector<double> k(static_cast<std::size_t>(n*n)),m(k.size()),d(static_cast<std::size_t>(n),1.);
 std::vector<numerics::linear::UnitDimensions> units(d.size());
 for (core::Index j=0;j<n;++j) for(core::Index i=0;i<n;++i) {
  k[static_cast<std::size_t>(j*n+i)]=i==j?3.+static_cast<double>(i):.01/static_cast<double>(i+j+1);
  m[static_cast<std::size_t>(j*n+i)]=i==j?2.:.05;
 }
 double ls=1,ms=1;
 if(n==2){k={-2.,-2.,-2.,10.};m={2.,2.,2.,5.};d={3.,.2};ls=7;ms=3;}
 core::ResourcePlan inputs(1000000);
 auto stiffness=numerics::DenseReal::copy({n,n,n},k,inputs,core::Operation::ExtractEigen);
 auto mass=numerics::DenseReal::copy({n,n,n},m,inputs,core::Operation::ExtractEigen);
 if(!stiffness||!mass)throw std::runtime_error("input allocation");
 const eig::GeneralizedRequest request{{{d,ls},ms,units,{},{},context()},
  {{1.e-12,core::ObservationId{20},core::ScopeId{21}},{1.e-12,core::ObservationId{22},core::ScopeId{23}},{1.e-12,core::ObservationId{24},core::ScopeId{25}},eig::NotRequested{}},
  {core::Precision::Binary64,numerics::linear::RetainedPrecision{core::NumericalGeneration{40}},{},nullptr},
  {n,n,{static_cast<std::uint64_t>(n),0,0,1},numerics::linear::EstimateAdmission{100000000}},{},false,eig::GeneralizedMeasurement::Enabled};
 std::array<core::StageTiming,5> total{};
 std::uint64_t payload=0,checksum=0;
 for(unsigned repeat=0;repeat<repetitions+1;++repeat) {
  core::ResourcePlan plan(10000000);
  auto result=eig::symmetric_generalized_spectrum({*stiffness,context()},{*mass,context()},request,plan);
  if(!result)throw std::runtime_error("extraction failed");
  if(repeat==0)continue;
  for(std::size_t s=0;s<5;++s){total[s].calls+=result->timings().stages[s].calls;total[s].elapsed_nanoseconds+=result->timings().stages[s].elapsed_nanoseconds;}
  payload=plan.allocated_bytes();
  const auto current=std::bit_cast<std::uint64_t>(result->eigenvalues()[0]);
  if(repeat>1&&current!=checksum)throw std::runtime_error("changed numerical bytes");checksum=current;
 }
 const auto memory=core::sample_process_memory(core::Operation::ExtractEigen);
 if(!memory)throw std::runtime_error("memory sampling");
 std::cout<<"{\"order\":"<<n<<",\"round\":"<<round<<",\"repetitions\":"<<repetitions<<",\"payload_bytes\":"<<payload<<",\"resident_bytes\":"<<memory->resident_bytes<<",\"process_peak_bytes\":"<<memory->peak_resident_bytes<<",\"first_root_bits\":"<<checksum<<",\"stages\":[";
 for(std::size_t s=0;s<5;++s){if(s)std::cout<<',';std::cout<<"{\"calls\":"<<total[s].calls<<",\"nanoseconds\":"<<total[s].elapsed_nanoseconds<<'}';}
 std::cout<<"]}\n";
}
int main() {
 try {
  std::cout<<build_info_json()<<'\n';
  for(unsigned round=0;round<3;++round){sample(2,100,round);sample(32,20,round);}
  return 0;
 }catch(...){return 1;}
}
