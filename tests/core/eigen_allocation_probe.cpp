#include "../support/allocation_audit.hpp"

#include <string_view>
#if defined(_WIN32)
#include <crtdbg.h>
#include <cstdlib>
#endif

int main(int argc, char* argv[]) {
#if defined(_WIN32)
    _set_abort_behavior(0, _WRITE_ABORT_MSG | _CALL_REPORTFAULT);
    _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDERR);
#endif
    if (argc != 2) {
        return 2;
    }
    const std::string_view mode{argv[1]};
    Eigen::VectorXd prepared = Eigen::VectorXd::Ones(4);
    uvlas::test_support::AllocationAudit audit{
        uvlas::test_support::strict_kernel_backends,
        uvlas::test_support::AuditThreadDomain::CallingThread};
    if (mode == "clean") {
        prepared *= 2;
        return audit.finish().status == uvlas::test_support::AllocationAuditStatus::Clean ? 0 : 3;
    }
    if (mode == "allocate") {
        Eigen::VectorXd unexpected(8);
        unexpected.setZero();
        return 4;
    }
    if (mode == "reallocate") {
        prepared.conservativeResize(9);
        return 4;
    }
    return 2;
}
