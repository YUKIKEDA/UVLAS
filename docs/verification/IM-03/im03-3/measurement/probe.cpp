#include "aero/local_geometry_packet.hpp"
#include "geometry_snapshot_fixture.hpp"
#include <uvlas/core/build_info.hpp>
#include <bit>
#include <chrono>
#include <cmath>
#include <iostream>
#if defined(_WIN32)
#include <windows.h>
#endif
using namespace uvlas;
using namespace uvlas::test_geometry;
using namespace subject;
int main(){
#if defined(_WIN32)
    if(!SetProcessAffinityMask(GetCurrentProcess(),1))return 2;
#endif
    std::cout<<build_info_json()<<'\n';
    for(const int scenario:{0,1,2,3}){
        TriangleFixture fixture;
        if(scenario==1){fixture.reference={0,0,0,1,0,0,1,0x1p-60,0};fixture.current=fixture.reference;}
        if(scenario>=2){const int exponent=scenario==2?-900:900;for(auto& x:fixture.reference)x=std::scalbn(x,exponent);fixture.current=fixture.reference;fixture.radius=SingleWord{std::scalbn(.25,exponent)};}
        core::ResourcePlan setup{budget};auto source_workspace=workspace(setup);
        auto work=LocalPacketWorkspace::create(setup);if(!work)return 3;
        auto geometry=take(prepare_owned_geometry(fixture.input(),limits(),*source_workspace,setup));
        auto query=take(resolve_query(geometry,{QueryId{5},SourceFacePoint{FacetId{70},SingleWord{.25},SingleWord{.25}}},geometry.context(),limits(),*source_workspace,setup));
        {core::ResourcePlan warmup{budget};auto result=prepare_local_packet(geometry,query,FacetId{70},{core::Precision::TwoWord,1e-10,true},limits(),**work,warmup);if(!result)return 4;}
        for(int sample=0;sample<3;++sample){
            core::ResourcePlan output{budget};
            const auto begin=std::chrono::steady_clock::now();
            auto result=prepare_local_packet(geometry,query,FacetId{70},{core::Precision::TwoWord,1e-10,true},limits(),**work,output);
            const auto wall=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-begin).count();
            if(!result)return 5;
            const auto timings=take(result->timings());if(!timings)return 6;
            const auto memory=core::sample_process_memory(core::Operation::EvaluateAero);if(!memory)return 7;
            std::cout<<"{\"scenario\":"<<scenario<<",\"sample\":"<<sample<<",\"wall_ns\":"<<wall
                <<",\"classification_ns\":"<<timings->classification.elapsed_nanoseconds
                <<",\"qualification_ns\":"<<timings->qualification.elapsed_nanoseconds
                <<",\"enclosure_ns\":"<<timings->enclosure.elapsed_nanoseconds
                <<",\"root_ns\":"<<timings->root.elapsed_nanoseconds
                <<",\"width_ns\":"<<timings->width.elapsed_nanoseconds
                <<",\"workspace_payload\":"<<LocalPacketWorkspace::payload_bytes()
                <<",\"result_plan_payload\":"<<output.allocated_bytes()
                <<",\"observed_steps\":"<<take(result->steps())
                <<",\"rss\":"<<memory->resident_bytes<<",\"lifetime_peak_rss\":"<<memory->peak_resident_bytes
                <<",\"words\":[";
            for(unsigned i=0;i<12;++i){const auto& value=take(result->value(static_cast<LocalQuantity>(i))).get();const auto words=std::get<DoubleWord>(value.center());
                if(i!=0)std::cout<<',';
                std::cout<<"[\""<<std::hex<<std::bit_cast<std::uint64_t>(words.high)<<"\",\""<<std::bit_cast<std::uint64_t>(words.low)<<"\",\""<<std::bit_cast<std::uint64_t>(value.radius())<<"\","<<std::dec<<value.exponent()<<']';}
            std::cout<<"]}\n";
        }
    }
}
