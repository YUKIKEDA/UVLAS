#pragma once
#include "case_schema_test_adapter.hpp"
#include "reference_declarations.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace uvlas::test::declaration_fixture {
using namespace schema_adapter;
namespace rd = io::detail;
constexpr std::uint64_t budget = std::uint64_t{1} << 40U;
template <class T>
T take(core::Result<T> value) {
    if (!value)
        throw std::runtime_error("Declaration control failed: " +
                                 std::to_string(static_cast<unsigned>(value.error().code)));
    return std::move(*value);
}
inline rd::DeclarationLimits limits() {
    std::array<std::uint64_t, static_cast<std::size_t>(rd::DeclarationQuota::Count)> caps{};
    caps.fill(budget);
    return rd::DeclarationLimits{caps};
}
inline rd::DeclarationLimits with(rd::DeclarationQuota q, std::uint64_t value) {
    std::array<std::uint64_t, static_cast<std::size_t>(rd::DeclarationQuota::Count)> caps{};
    caps.fill(budget);
    caps[static_cast<std::size_t>(q)] = value;
    return rd::DeclarationLimits{caps};
}
inline io::InputDocument parsed(const JsonToken& value) {
    const auto wire = value.wire();
    const io::InputLimits syntax{
        10000000, 256,      1000000,  1000000,
        10000000, 10000000, 10000000, {io::ParserAllocationMode::EstimateAdmission, budget}};
    core::ResourcePlan plan{budget};
    std::optional<io::InputDiagnostic> diag;
    auto doc = take(io::read_input_document(std::as_bytes(std::span{wire.data(), wire.size()}),
                                            io::SourceId{12}, syntax, plan, diag));
    return doc;
}
inline io::SchemaCase checked(const JsonToken& value) {
    const auto doc = parsed(value);
    core::ResourcePlan plan{budget};
    std::array<std::uint64_t, static_cast<std::size_t>(io::SchemaQuota::Count)> caps{};
    caps.fill(budget);
    io::SchemaDiagnostic slot;
    auto result = io::validate_case_schema(doc, io::SchemaLimits{caps}, plan, &slot);
    if (!result)
        throw std::runtime_error("Schema control failed: " +
                                 std::string{slot.path().value_or("<no-path>")});
    return std::move(*result);
}
inline JsonToken json(std::string_view text) {
    return TokenReader{text}.read();
}
inline void reference(JsonToken& file, std::string_view id, std::string_view path, char digit) {
    file.set("source_id", json(nlohmann::json(id).dump()));
    file.set("path", json(nlohmann::json(path).dump()));
    file.set("sha256", json(nlohmann::json(std::string(64, digit)).dump()));
}
// Eight independently specified roles; ID order intentionally differs from declaration order.
inline JsonToken all_roles() {
    auto value = root("coupled_transient");
    auto gmsh = fragment("Geometry.gmsh");
    reference(gmsh.field("mesh"), "0300", R"(folder\./mesh.msh)", '3');
    reference(gmsh.field("sidecar"), "301", "folder/../sidecar.json", '4');
    reference(gmsh.field("reference").field("file"), "302", "参照/geometry.json", '5');
    auto nastran = fragment("Geometry.nastran");
    nastran.set("id", json("\"3\""));
    reference(nastran.field("deck"), "200", "model.bdf", '1');
    reference(nastran.field("sidecar"), "201", "model-sidecar.json", '2');
    auto& geometries = value.field("model").field("geometry").elements;
    geometries.insert(geometries.begin(), nastran);
    geometries.insert(geometries.begin(), gmsh);
    auto seed = fragment("BranchSeed.prepared");
    reference(seed.field("history"), "700", "seed.history", '6');
    value.field("analysis").field("branch").set("seed", seed);
    auto comparison = fragment("ComparisonReference.file");
    reference(comparison.field("source"), "900", "observation.dat", '8');
    comparison.set("columns", json("[\"x\",\"y\",\"z\"]"));
    comparison.set("uncertainty", json("[0,0,0]"));
    value.field("observations").elements[0].field("comparison").set("reference", comparison);
    auto initial = fragment("Initial.prepared_history");
    reference(initial.field("file"), "800", "initial.history", '7');
    value.set("initial", initial);
    return value;
}
inline void reverse_keys(JsonToken& value) {
    std::ranges::reverse(value.members);
    for (auto& member : value.members)
        reverse_keys(member.second);
    for (auto& child : value.elements)
        reverse_keys(child);
}
inline bool same_identity(const io::SchemaIdentity& a, const io::SchemaIdentity& b) {
    return a.source == b.source && a.contract_revision == b.contract_revision &&
           a.representation_revision == b.representation_revision;
}
inline std::string snapshot(const rd::DeclaredReferences& batch) {
    const auto identity = take(batch.identity());
    std::string result = std::to_string(identity.source.source.value()) + ":" +
                         std::to_string(identity.source.byte_count);
    for (auto byte : identity.source.digest.bytes)
        result += std::to_string(std::to_integer<unsigned>(byte)) + ",";
    result += std::to_string(identity.contract_revision) + ":" +
              std::to_string(identity.representation_revision);
    for (core::Index i = 0; i < static_cast<core::Index>(take(batch.size())); ++i) {
        const auto r = take(batch.entry(i));
        result += "|" + std::to_string(r.source_id.value()) + ":" +
                  std::to_string(static_cast<unsigned>(r.role)) + ":" +
                  std::to_string(r.declaration_node);
        for (auto b : r.declared_digest.bytes)
            result += std::to_string(std::to_integer<unsigned>(b)) + ",";
        for (auto s :
             {r.source_id_spelling, r.original_path, r.digest_spelling, r.declaration_pointer})
            result += std::to_string(s.size()) + ":" + std::string{s};
    }
    const auto m = take(batch.memory_report());
    for (const auto n :
         {m.owned_bytes, m.scratch_bytes, m.cumulative_bytes, m.work, m.planned_owned_bytes,
          m.planned_scratch_bytes, m.planned_work, m.admission_bytes})
        result += "/" + std::to_string(n);
    return result;
}
} // namespace uvlas::test::declaration_fixture
