#include "aero/p1_induction.hpp"
#include <iostream>
int main() {
    using namespace uvlas::aero::detail;
    std::cout << "workspace_payload_bytes=" << P1Workspace::payload_bytes() << '\n'
              << "result_payload_bytes=" << P1Basis::payload_bytes() << '\n';
}