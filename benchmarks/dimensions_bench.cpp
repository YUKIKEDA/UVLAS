#include <array>
#include <benchmark/benchmark.h>
#include <uvlas/core/checked_dimensions.hpp>

namespace {
void checked_extent_product(benchmark::State& state) {
    const std::array<uvlas::core::Index, 3> dimensions{5000, 64, 3};
    for (auto iteration : state) {
        (void)iteration;
        auto result =
            uvlas::core::checked_extent_product(dimensions, uvlas::core::Operation::PrepareModel);
        benchmark::DoNotOptimize(result);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(checked_extent_product);
} // namespace
