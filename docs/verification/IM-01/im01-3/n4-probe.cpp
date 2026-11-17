#include <uvlas/structure/pg_coefficients.hpp>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace uvlas;using namespace uvlas::structure;
template<class T>T take(core::Result<T> value){if(!value)throw std::runtime_error("Failure "+std::to_string(static_cast<unsigned>(value.error().code)));return std::move(*value);}
int main(){try{
 constexpr auto budget=std::uint64_t{64}*1024*1024*1024;
 const PgCoefficientLimits limits{{64,1000,4096,32,32,65536,2000000,ExactEstimateAdmission{std::uint64_t{16}*1024*1024*1024},ExactTimingMode::Enabled},1000,20000000,2000,{}};
 core::ResourcePlan input{budget};auto input_limits=limits.algebra;input_limits.timing=ExactTimingMode::Disabled;
 const auto number=[&](std::int64_t p,std::int64_t q){return take(ExactRational::from_integer_pair(std::to_string(p),std::to_string(q),input_limits,input));};
 const auto x0=number(2,7),w=number(4,5),a=number(3,20);
 std::array<ExactRational,15> coefficients;std::vector<ExactInputTerm> terms;terms.reserve(15);std::size_t ordinal=0;
 for(core::Index i=0;i<=4;++i)for(core::Index j=0;j<=4-i;++j){coefficients[ordinal]=number(i+2*j+1,i+j+1);terms.push_back({{i,j},coefficients[ordinal++]});}
 const auto f=take(ExactChartPolynomial::make(terms,input_limits,input));
 const auto graph=take(PgGraph::prepare({x0,w,a,f,4,PgGraphConvention::Version1},{input_limits,64,{}},input));
 core::Workspace workspace{core::PlannedBuffer<std::byte>{}};core::ResourcePlan output{budget};
 const auto before_memory=take(core::sample_process_memory(core::Operation::PrepareModel));
 const auto start=std::chrono::steady_clock::now();
 const auto value=take(PgCoefficientDictionary::prepare(graph,limits,workspace,output));
 const auto finish=std::chrono::steady_clock::now();
 const auto memory=take(core::sample_process_memory(core::Operation::PrepareModel));
 const auto& info=value.preparation();const auto& times=info.timings;if(!times)throw std::runtime_error("Missing product timings");
 const auto whole=std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count();
 const auto summed=times->coefficients.elapsed_nanoseconds+times->anchor.elapsed_nanoseconds+times->verification.elapsed_nanoseconds;
 if(whole<0 || static_cast<std::uint64_t>(whole)<summed || info.observed_columns!=933)throw std::runtime_error("Invalid observation");
 std::cout<<"{\"schema\":\"uvlas.pg-dictionary-stage-observation/1\",\"n\":4,\"samples\":1,\"whole_nanoseconds\":"<<whole
 <<",\"coefficients_nanoseconds\":"<<times->coefficients.elapsed_nanoseconds<<",\"anchor_nanoseconds\":"<<times->anchor.elapsed_nanoseconds
 <<",\"verification_nanoseconds\":"<<times->verification.elapsed_nanoseconds<<",\"outside_three_phases_nanoseconds\":"<<static_cast<std::uint64_t>(whole)-summed
 <<",\"logical_columns\":"<<info.observed_columns<<",\"stored_terms\":"<<info.observed_nonzero_terms<<",\"observed_maximum_value_bits\":"<<info.observed_maximum_value_bits
 <<",\"owned_cumulative_actual_bytes\":"<<info.allocated_owned_payload_bytes<<",\"owned_cumulative_forecast_bytes\":"<<info.planned_owned_payload_upper_bytes
 <<",\"external_payload_estimate_bytes\":"<<info.arithmetic.estimated_external_payload_bytes<<",\"resident_before_bytes\":"<<before_memory.resident_bytes
 <<",\"resident_after_bytes\":"<<memory.resident_bytes<<",\"lifetime_peak_bytes\":"<<memory.peak_resident_bytes<<"}\n";
 return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
