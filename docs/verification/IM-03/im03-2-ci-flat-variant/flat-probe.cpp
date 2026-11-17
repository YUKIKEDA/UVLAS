#include <array>
#include <variant>
#include <cstdint>
#include <type_traits>
#include <cstdio>
struct NoMaterialSupport{};
using ParameterStorage = std::variant<double, std::array<double, 2>>;
using FreeSingleStorage = std::array<double, 3>;
using FreeDoubleStorage = std::array<double, 6>;
struct VertexStorage {
    std::uint64_t point = 0;
};
struct EdgeQueryStorage {
    std::uint64_t edge = 0;
    ParameterStorage t;
};
struct FaceQueryStorage {
    std::uint64_t facet = 0;
    ParameterStorage u, v;
};
struct PointSupportStorage {
    std::uint64_t id = 0;
};
struct EdgeSupportStorage {
    std::uint64_t id = 0;
};
struct FacetSupportStorage {
    std::uint64_t id = 0;
};
struct QueryStorage {
    std::uint64_t id = 0;
    std::variant<FreeSingleStorage, FreeDoubleStorage, VertexStorage, EdgeQueryStorage,
                 FaceQueryStorage>
        definition;
    std::variant<NoMaterialSupport, PointSupportStorage, EdgeSupportStorage, FacetSupportStorage>
        support;
};
// These selected alternatives construct/copy/move without a throwing inner
// variant conversion, so no outer strong-guarantee temporary is needed.
static_assert(std::is_nothrow_copy_constructible_v<FreeSingleStorage>);
static_assert(std::is_nothrow_copy_constructible_v<FreeDoubleStorage>);
static_assert(std::is_nothrow_move_constructible_v<EdgeQueryStorage>);
static_assert(std::is_nothrow_copy_constructible_v<EdgeQueryStorage>);
static_assert(std::is_nothrow_move_constructible_v<FaceQueryStorage>);
static_assert(std::is_nothrow_copy_constructible_v<FaceQueryStorage>);
static_assert(std::is_nothrow_copy_assignable_v<QueryStorage>);

extern void consume(const QueryStorage&);
void producer(int kind,const std::array<double,3>& high,const std::array<double,3>& low,std::uint64_t id,double t) {
 QueryStorage query;query.id=id;
 switch(kind){case 0: query.definition.emplace<FreeSingleStorage>(high); break;case 1: query.definition.emplace<FreeDoubleStorage>(std::array{high[0],high[1],high[2],low[0],low[1],low[2]}); break;
 case 2:query.definition=VertexStorage{id};break;
 case 3:query.definition=EdgeQueryStorage{id,ParameterStorage{t}};break;
 case 4:query.definition=FaceQueryStorage{id,ParameterStorage{t},ParameterStorage{std::array{t,t}}};break;
 default:return;}
 consume(query);
}
