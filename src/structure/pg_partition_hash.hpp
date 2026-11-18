#pragma once
#include "pg_partition_topology.hpp"

#include <boost/hash2/sha2.hpp>

namespace uvlas::structure::detail {
struct PgPartitionHash {
    PgPartitionAlgebra& op;
    boost::hash2::sha2_256 hash;
    void bytes(const void* data, std::size_t size) {
        op.visit(size);
        if (size != 0)
            hash.update(data, size);
    }
    void u64(std::uint64_t value) {
        std::array<unsigned char, 8> out{};
        for (std::size_t i = 0; i < 8; ++i) {
            out[7 - i] = static_cast<unsigned char>(value & 0xffU);
            value >>= 8U;
        }
        bytes(out.data(), out.size());
    }
    void integer(const IntegerBackend& value) {
        const auto bits =
            value.size() == 1 && value.limbs()[0] == 0
                ? std::uint64_t{0}
                : static_cast<std::uint64_t>(
                      (value.size() - 1) * IntegerBackend::limb_bits + IntegerBackend::limb_bits -
                      static_cast<std::size_t>(std::countl_zero(value.limbs()[value.size() - 1])));
        const auto count = (bits + 7) / 8;
        u64(count);
        for (auto i = count; i != 0; --i) {
            const auto bit = (i - 1) * 8;
            const auto byte = static_cast<unsigned char>(
                (value.limbs()[static_cast<std::size_t>(bit / IntegerBackend::limb_bits)] >>
                 static_cast<unsigned>(bit % IntegerBackend::limb_bits)) &
                0xffU);
            bytes(&byte, 1);
        }
    }
    void rational(const ExactRational& value) {
        const auto& b = ExactFactory::value(value).backend();
        const unsigned char sign = b.num().sign() ? 1 : 0;
        bytes(&sign, 1);
        integer(b.num());
        integer(b.denom());
    }
    void operator()(PgPartitionData& d) {
        constexpr std::string_view tag = "uvlas.pg.partition.v1";
        u64(tag.size());
        bytes(tag.data(), tag.size());
        u64(1);
        const auto graph = partition_graph(d);
        bytes(graph.identity.bytes.data(), graph.identity.bytes.size());
        u64(d.domain_count);
        u64(d.thickness_offset);
        u64(d.source_points.size());
        for (const auto& point : d.source_points.values()) {
            rational(point.s);
            rational(point.t);
        }
        u64(d.thickness_triangles.size());
        for (const auto& triangle : d.thickness_triangles.values())
            for (const auto i : triangle)
                u64(static_cast<std::uint64_t>(i));
        u64(d.thickness_values.size());
        for (const auto& value : d.thickness_values.values())
            rational(value);
        u64(d.layers.size());
        for (const auto& layer : d.layers.values()) {
            u64(layer.id.value());
            rational(layer.lower);
            rational(layer.upper);
        }
        u64(d.boundaries.size());
        for (const auto& boundary : d.boundaries.values()) {
            u64(static_cast<std::uint64_t>(boundary.first_point));
            u64(static_cast<std::uint64_t>(boundary.point_count));
            u64(static_cast<std::uint64_t>(boundary.text_count));
            const auto text =
                d.text.values().subspan(static_cast<std::size_t>(boundary.text_offset),
                                        static_cast<std::size_t>(boundary.text_count));
            bytes(text.data(), text.size());
        }
        u64(d.sources.size());
        for (const auto& source : d.sources.values()) {
            u64(source.descriptor.index());
            for (const auto i : source.endpoints)
                u64(static_cast<std::uint64_t>(i));
            if (const auto* r = std::get_if<PgPartitionRefinementSource>(&source.descriptor)) {
                u64(static_cast<std::uint64_t>(r->ordinal));
                u64(r->id);
                u64(r->version);
            } else if (const auto* t =
                           std::get_if<PgPartitionThicknessSource>(&source.descriptor)) {
                u64(static_cast<std::uint64_t>(t->triangle));
                u64(static_cast<std::uint64_t>(t->edge));
            } else if (const auto* n = std::get_if<PgPartitionNamedSource>(&source.descriptor)) {
                u64(static_cast<std::uint64_t>(n->boundary));
                u64(static_cast<std::uint64_t>(n->segment));
            } else
                u64(static_cast<std::uint64_t>(
                    std::get<PgPartitionDomainSource>(source.descriptor).edge));
        }
        const auto digest = hash.result();
        for (std::size_t i = 0; i < d.identity.bytes.size(); ++i)
            d.identity.bytes[i] = static_cast<std::byte>(digest[i]);
    }
};
} // namespace uvlas::structure::detail
