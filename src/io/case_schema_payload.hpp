#pragma once

#include "case_schema_storage.hpp"

namespace uvlas::io::detail {

// The same field-lowering kernel serves the complete Case and its private contract
// tests. A scoped payload is a distinct internal result, never a public SchemaCase.
struct SchemaPayloadRequest {
    SchemaPayloadRequest(InputNodeView node, FieldType declared, std::uint64_t count,
                         SchemaDimension dimension, bool registry) noexcept
        : input(node),
          type(declared),
          components(count),
          dimensions(dimension),
          case_registry(registry) {}
    InputNodeView input;
    FieldType type;
    std::uint64_t components;
    SchemaDimension dimensions;
    bool case_registry;
};

struct SchemaPayload {
    SchemaPayload(core::PlannedBuffer<SchemaData> owner, std::uint64_t node,
                  FieldType declared) noexcept
        : data(std::move(owner)), root(node), type(declared) {}
    core::PlannedBuffer<SchemaData> data;
    std::uint64_t root;
    FieldType type;
};

[[nodiscard]] core::Result<SchemaPayload> validate_schema_payload(const InputDocument&,
                                                                  SchemaPayloadRequest,
                                                                  SchemaLimits, core::ResourcePlan&,
                                                                  SchemaDiagnostic*) noexcept;

} // namespace uvlas::io::detail
