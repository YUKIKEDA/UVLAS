#pragma once

#include "gmsh_mesh.hpp"

#include <algorithm>

namespace uvlas::io::detail {

// Entity membership uses the original dimension and PhysicalTag, independent of names.
[[nodiscard]] inline bool gmsh_element_in_physical_group(const GmshMeshView& mesh,
                                                         const GmshElement& element, GmshKey group,
                                                         auto&& visit) {
    if (element.entity.dimension != group.dimension) {
        return false;
    }
    const auto compare = [&](const auto& lhs, const auto& rhs) {
        visit();
        return lhs < rhs;
    };
    const auto entity =
        std::ranges::lower_bound(mesh.entities, element.entity, compare, &GmshEntity::key);
    if (entity == mesh.entities.end() || entity->key != element.entity) {
        return false;
    }
    const auto links = mesh.entity_links.subspan(entity->first_physical, entity->physical_count);
    return std::ranges::any_of(links, [&](const auto& link) {
        visit();
        return link.key == group;
    });
}
} // namespace uvlas::io::detail
