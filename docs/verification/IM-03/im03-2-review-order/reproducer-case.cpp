TEST_CASE("Connection diagnostics select the lowest facet ID independently of source row order",
          "[UI-D05][UI-D06][N4-03][negative][verification]") {
    const TriangleFixture fixture;
    std::array rows{
        TriangleInput{
            FacetId{90}, PhysicalFaceId{999}, fixture.points, {EdgeId{10}, EdgeId{20}, EdgeId{30}}},
        TriangleInput{FacetId{10},
                      PhysicalFaceId{999},
                      fixture.points,
                      {EdgeId{10}, EdgeId{20}, EdgeId{30}}}};
    core::ResourcePlan prep{budget};
    auto w = workspace(prep);
    for (unsigned order = 0; order < 2; ++order) {
        auto input = fixture.input();
        input.triangles = rows;
        core::ResourcePlan plan{budget};
        const auto rejected = prepare_owned_geometry(input, limits(), *w, plan);
        failure(rejected, core::FailureCode::InvalidInput);
        if (rejected || !rejected.error().problem) {
            FAIL_CHECK("Missing typed facet failure");
            return;
        }
        const auto* problem = std::get_if<FacetId>(&*rejected.error().problem);
        if (problem == nullptr) {
            FAIL_CHECK("Wrong typed geometry problem");
            return;
        }
        CHECK(*problem == FacetId{10});
        std::swap(rows[0], rows[1]);
    }
}
