#include "case_schema_detail.hpp"
#include "case_schema_expected_fields.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <set>
#include <string>
#include <string_view>
#include <uvlas/io/case_schema.hpp>

TEST_CASE("Adopted field inventory independently fixes names types and required flags",
          "[verification][negative][IO4-01]") {
    using namespace uvlas;
    const auto records = io::detail::schema_records();
    const auto types = io::detail::schema_types();
    CHECK(records.size() == 194);
    std::size_t actual_members = 0;
    for (const auto& record : records)
        actual_members += record.fields.size();
    CHECK(actual_members == test::schema_field_inventory::fields.size());
    std::set<std::string> expected_records;
    for (const auto& expected : test::schema_field_inventory::fields) {
        INFO(expected.record << "/" << expected.name);
        expected_records.insert(std::string{expected.record});
        const auto record = std::ranges::find_if(
            records, [&](const auto& row) { return row.name == expected.record; });
        REQUIRE(record != records.end());
        const auto field = std::ranges::find_if(
            record->fields, [&](const auto& row) { return row.name == expected.name; });
        REQUIRE(field != record->fields.end());
        CHECK(field->optional == expected.optional);
        const auto& actual_type = types[static_cast<std::size_t>(field->type)];
        if (expected.record == "FileReference" && expected.name == "role") {
            // CS-02's role parameter is instantiated with all eight FD file purposes.
            CHECK(actual_type.kind == io::detail::FieldKind::Enumeration);
            constexpr std::array roles{"nastran-deck",       "nastran-sidecar",
                                       "gmsh-msh4.1-ascii",  "gmsh-sidecar",
                                       "reference-geometry", "preparation-history",
                                       "initial-history",    "observation-reference"};
            std::set<std::string> actual_roles;
            auto spelling = actual_type.spelling;
            REQUIRE(spelling.starts_with("Tag{"));
            REQUIRE(spelling.ends_with('}'));
            spelling = spelling.substr(4, spelling.size() - 5);
            for (;;) {
                const auto end = spelling.find(',');
                actual_roles.insert(std::string{spelling.substr(0, end)});
                if (end == std::string_view::npos)
                    break;
                spelling.remove_prefix(end + 1);
            }
            CHECK(actual_roles.size() == roles.size());
            for (const auto role : roles)
                CHECK(actual_roles.contains(role));
        } else
            CHECK(actual_type.spelling == expected.spelling);
    }
    CHECK(expected_records.size() == records.size());
}
