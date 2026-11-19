#include "reference_path.hpp"

#include <algorithm>
#include <utility>

namespace uvlas::io::detail {
namespace {

constexpr auto operation = core::Operation::PrepareModel;

[[nodiscard]] bool is_separator(char value) noexcept {
    return value == '/' || value == '\\';
}

[[nodiscard]] bool is_drive_prefix(std::string_view text) noexcept {
    if (text.size() < 2 || text[1] != ':') {
        return false;
    }
    return (text.front() >= 'A' && text.front() <= 'Z') ||
           (text.front() >= 'a' && text.front() <= 'z');
}

[[nodiscard]] std::optional<ReferencePathDiagnostic> check_path(std::string_view input) noexcept {
    const auto nul = input.find('\0');
    if (nul != input.npos) {
        return ReferencePathDiagnostic{ReferencePathReason::InvalidPath, nul};
    }
    if (input.empty()) {
        return ReferencePathDiagnostic{ReferencePathReason::MissingFile, 0};
    }
    if (is_separator(input.front()) || is_drive_prefix(input)) {
        return ReferencePathDiagnostic{ReferencePathReason::RootedPath, 0};
    }
    const auto colon = input.find(':');
    if (colon != input.npos) {
        return ReferencePathDiagnostic{ReferencePathReason::StreamPath, colon};
    }
    if (is_separator(input.back())) {
        return ReferencePathDiagnostic{ReferencePathReason::MissingFile, input.size() - 1};
    }

    // Check the prefix depth in forward order so the first escape has a stable byte offset.
    std::size_t depth = 0;
    for (std::size_t begin = 0; begin < input.size();) {
        const auto end = std::min(input.find_first_of("/\\", begin), input.size());
        const auto component = input.substr(begin, end - begin);
        if (component == "..") {
            if (depth == 0) {
                return ReferencePathDiagnostic{ReferencePathReason::EscapeRoot, begin};
            }
            --depth;
        } else if (!component.empty() && component != ".") {
            ++depth;
        }
        if (end == input.size()) {
            break;
        }
        begin = end + 1;
    }
    if (depth == 0) {
        return ReferencePathDiagnostic{ReferencePathReason::MissingFile, std::nullopt};
    }
    return std::nullopt;
}

template <class Visit>
void visit_retained_components(std::string_view input, Visit visit) noexcept {
    // After the forward check, each ".." can cancel one preceding component. Walking backwards
    // identifies survivors in linear work without allocating a stack of removed components.
    std::size_t pending_parents = 0;
    for (auto end = input.size(); end != 0;) {
        auto begin = end;
        while (begin != 0 && !is_separator(input[begin - 1])) {
            --begin;
        }
        const auto component = input.substr(begin, end - begin);
        if (component == "..") {
            ++pending_parents;
        } else if (!component.empty() && component != ".") {
            if (pending_parents == 0) {
                visit(component);
            } else {
                --pending_parents;
            }
        }
        end = begin == 0 ? 0 : begin - 1;
    }
}

} // namespace

NormalizedReferencePath::NormalizedReferencePath(core::PlannedBuffer<char> bytes) noexcept
    : bytes_(std::move(bytes)) {}

std::string_view NormalizedReferencePath::text() const& noexcept {
    const auto bytes = bytes_.values();
    return bytes.empty() ? std::string_view{} : std::string_view{bytes.data(), bytes.size()};
}

core::Result<NormalizedReferencePath>
normalize_reference_path(std::string_view input, ReferencePathLimits limits,
                         core::ResourcePlan& plan, const std::stop_token& stop,
                         std::optional<ReferencePathDiagnostic>* diagnostic) noexcept {
    if (diagnostic != nullptr) {
        diagnostic->reset();
    }
    const auto reject = [&](core::FailureCode code, ReferencePathReason reason,
                            std::optional<std::size_t> offset = std::nullopt) noexcept {
        if (diagnostic != nullptr) {
            *diagnostic = ReferencePathDiagnostic{reason, offset};
        }
        return std::unexpected(core::Failure{code, operation});
    };
    if (stop.stop_requested()) {
        return std::unexpected(core::Failure{core::FailureCode::Cancelled, operation});
    }
    if (input.size() > limits.input_bytes) {
        return reject(core::FailureCode::ResourceExceeded, ReferencePathReason::InputBytes);
    }
    if (const auto problem = check_path(input)) {
        return reject(core::FailureCode::InvalidInput, problem->reason, problem->input_offset);
    }

    std::size_t depth = 0;
    std::size_t length = 0;
    visit_retained_components(input, [&](std::string_view component) noexcept {
        // Total survivor bytes plus their separating slashes cannot exceed input.size().
        length += component.size() + (depth != 0 ? 1U : 0U);
        ++depth;
    });
    if (length > limits.normalized_bytes || !std::in_range<core::Index>(length)) {
        return reject(core::FailureCode::ResourceExceeded, ReferencePathReason::NormalizedBytes);
    }
    if (depth > limits.normalized_components) {
        return reject(core::FailureCode::ResourceExceeded, ReferencePathReason::NormalizedDepth);
    }
    if (stop.stop_requested()) {
        return std::unexpected(core::Failure{core::FailureCode::Cancelled, operation});
    }
    auto bytes = plan.allocate<char>(static_cast<core::Index>(length), operation);
    if (!bytes) {
        if (diagnostic != nullptr) {
            *diagnostic = ReferencePathDiagnostic{ReferencePathReason::PlanBytes, std::nullopt};
        }
        return std::unexpected(bytes.error());
    }

    const auto output = bytes->values();
    auto cursor = length;
    visit_retained_components(input, [&](std::string_view component) noexcept {
        if (cursor != length) {
            output[--cursor] = '/';
        }
        cursor -= component.size();
        std::ranges::copy(component, output.subspan(cursor, component.size()).begin());
    });
    if (stop.stop_requested()) {
        return std::unexpected(core::Failure{core::FailureCode::Cancelled, operation});
    }
    return NormalizedReferencePath{std::move(*bytes)};
}

} // namespace uvlas::io::detail
