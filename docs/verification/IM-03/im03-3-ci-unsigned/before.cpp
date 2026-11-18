#include <cstdint>
#include <initializer_list>
int main() {
    std::uint64_t total=0;
    for (const std::uint64_t cap : {16, 17}) total+=cap;
    return total==33 ? 0 : 1;
}
