#pragma once

#include "pg_reference_arithmetic.hpp"

#include <boost/hash2/sha2.hpp>

#include <array>
#include <cstddef>
#include <string_view>

namespace uvlas::structure::detail {

/// Canonical signed-magnitude rational encoding; inputs are already canonical.
class ReferenceHash {
public:
    explicit ReferenceHash(ReferenceArithmetic& arithmetic) noexcept : math_(arithmetic) {}

    void bytes(const void* data, std::size_t size) {
        math_.cancel();
        hash_.update(data, size);
    }

    void rational(const ExactRational& value) {
        const auto& backend = ExactFactory::value(value).backend();
        const unsigned char sign = backend.num().sign() ? 1 : 0;
        bytes(&sign, 1);
        magnitude(backend.num(), value.numerator_bits());
        magnitude(backend.denom(), value.denominator_bits());
    }

    void word(std::uint64_t value) {
        std::array<unsigned char, 8> encoding{};
        for (auto i = encoding.rbegin(); i != encoding.rend(); ++i) {
            *i = static_cast<unsigned char>(value & 0xffU);
            value >>= 8U;
        }
        math_.charge();
        bytes(encoding.data(), encoding.size());
    }

    template <class Backend>
    void signed_magnitude(const Backend& value, std::uint64_t bits) {
        math_.charge();
        const unsigned char sign = value.sign() ? 1 : 0;
        bytes(&sign, 1);
        magnitude(value, bits);
    }

    [[nodiscard]] PgReferenceIdentity result() noexcept {
        PgReferenceIdentity result;
        const auto digest = hash_.result();
        for (std::size_t i = 0; i < result.bytes.size(); ++i)
            result.bytes[i] = static_cast<std::byte>(digest[i]);
        return result;
    }

private:
    template <class Backend>
    void magnitude(const Backend& value, std::uint64_t bits) {
        auto count = partition_add(bits, 7) / 8;
        std::array<unsigned char, 8> length{};
        auto remaining = count;
        for (auto i = length.rbegin(); i != length.rend(); ++i) {
            *i = static_cast<unsigned char>(remaining & 0xffU);
            remaining >>= 8U;
        }
        bytes(length.data(), length.size());
        for (; count > 0; --count) {
            math_.charge();
            const auto bit = (count - 1) * 8;
            constexpr auto limb_bits = Backend::limb_bits;
            const auto word = value.limbs()[static_cast<std::size_t>(bit / limb_bits)];
            const auto byte = static_cast<unsigned char>(
                (word >> static_cast<unsigned>(bit % limb_bits)) & 0xffU);
            bytes(&byte, 1);
        }
    }

    ReferenceArithmetic& math_;
    boost::hash2::sha2_256 hash_;
};

[[nodiscard]] inline PgReferenceIdentity make_reference_identity(ReferenceArithmetic& math,
                                                                 PgGraphContentIdentity graph,
                                                                 PgPartitionIdentity partition,
                                                                 const PgReferenceFrame& frame) {
    auto stream = ReferenceHash{math};
    constexpr std::string_view tag = "uvlas.pg.reference-geometry.v1";
    stream.bytes(tag.data(), tag.size());
    stream.bytes(graph.bytes.data(), graph.bytes.size());
    stream.bytes(partition.bytes.data(), partition.bytes.size());
    for (const auto& value : frame.rotation)
        stream.rational(value);
    for (const auto& value : frame.origin_metres)
        stream.rational(value);
    stream.rational(frame.length_metres);
    return stream.result();
}

} // namespace uvlas::structure::detail
