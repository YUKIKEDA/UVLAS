#include "p1_induction.hpp"
#include "aero/geometry_snapshot.hpp"
#include "aero/local_geometry_packet.hpp"
#include "aero/geometry_snapshot_fixture.hpp"
#include <chrono>
#include <bit>
#include <cstdio>
#include <cstdlib>
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
extern "C" void p1_profile_clear();
extern "C" const std::uint64_t* p1_profile_read();
using namespace uvlas;
using namespace uvlas::aero::detail;
template<class T> T get(T&& value) { return std::forward<T>(value); }
template<class T,class E> T take_result(std::expected<T,E> value) {
    if(!value) std::exit(100);
    return std::move(*value);
}
int main() {
    DWORD_PTR process_mask=0,system_mask=0;
    if(GetProcessAffinityMask(GetCurrentProcess(),&process_mask,&system_mask))
        SetProcessAffinityMask(GetCurrentProcess(),process_mask&(~process_mask+1));
    core::ResourcePlan setup{256ULL*1024*1024};
    auto geometry_workspace=take_result(GeometryWorkspace::create(setup));
    auto local_workspace=take_result(LocalPacketWorkspace::create(setup));
    auto kernel=take_result(P1Workspace::create(setup));
    const uvlas::test_geometry::TriangleFixture fixture;
    auto geometry=take_result(prepare_owned_geometry(fixture.input(),uvlas::test_geometry::limits(),*geometry_workspace,setup));
    for(unsigned which=0;which<2;++which) {
        auto query=take_result(resolve_query(geometry,{QueryId{44},FreePoint{FreeSingle{{.75,1.,which==0?2.:30.}}}},
            geometry.context(),uvlas::test_geometry::limits(),*geometry_workspace,setup));
        auto packet=take_result(prepare_local_packet(geometry,query,FacetId{70},{core::Precision::TwoWord,1e-10},
            uvlas::test_geometry::limits(),*local_workspace,setup));
        for(unsigned repeat=0;repeat<3;++repeat) {
            core::ResourcePlan output{P1Basis::payload_bytes()};
            p1_profile_clear();
            const auto start=std::chrono::steady_clock::now();
            auto result=take_result(evaluate_p1_basis(packet,{core::Precision::TwoWord,1e-10,1e-10},
                {numerics::linear::EstimateAdmission{28672},UINT64_C(1000000000)},*kernel,output));
            const auto total=std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now()-start).count();
            const auto* stages=p1_profile_read();
            std::printf("TIME %u %u %lld",which,repeat,static_cast<long long>(total));
            for(unsigned i=0;i<6;++i) std::printf(" %llu",static_cast<unsigned long long>(stages[i]));
            std::printf("\nBITS %u %u",which,repeat);
            for(std::size_t c=0;c<9;++c) for(std::size_t j=0;j<14;++j) {
                const auto& value=take_result(j==0?result.basis(c/3,c%3):result.derivative(c/3,c%3,j-1)).get();
                const auto words=value.center();
                double high=0,low=0;
                if(const auto* one=std::get_if<SingleWord>(&words)) high=one->value;
                else { const auto& pair=std::get<DoubleWord>(words); high=pair.high;low=pair.low; }
                std::printf(" %016llx %016llx %016llx %d",
                    static_cast<unsigned long long>(std::bit_cast<std::uint64_t>(high)),
                    static_cast<unsigned long long>(std::bit_cast<std::uint64_t>(low)),
                    static_cast<unsigned long long>(std::bit_cast<std::uint64_t>(value.radius())),value.exponent());
            }
            std::printf(" STEP %llu\n",static_cast<unsigned long long>(take_result(result.steps())));
        }
    }
}