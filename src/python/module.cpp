#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <array>
#include <bit>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <memory>
#include <mutex>
#include <optional>
#include <span>
#include <type_traits>
#include <utility>
#include <uvlas/core/state.hpp>
#include <uvlas/numerics/storage.hpp>
#include <variant>

namespace py = pybind11;
namespace core = uvlas::core;
namespace numerics = uvlas::numerics;
static_assert(std::is_nothrow_move_constructible_v<core::AcceptedState>);
namespace {

struct Reply {
    std::variant<py::object, core::Failure> outcome;
    explicit Reply(py::object value) : outcome(std::move(value)) {}
    explicit Reply(core::Failure failure) : outcome(failure) {}
};

// Python allocation failures retain their native exceptions. Native processing
// exceptions become owned diagnostics without borrowing exception text.
template <class F>
Reply boundary(core::Operation operation, F&& action) {
    try {
        return std::forward<F>(action)();
    } catch (const py::error_already_set&) {
        throw;
    } catch (const std::bad_alloc&) {
        return Reply{core::Failure{core::FailureCode::ResourceExceeded, operation}};
    } catch (...) {
        return Reply{core::Failure{core::FailureCode::InvalidInput, operation}};
    }
}
Reply invalid(core::Operation operation) {
    return Reply{core::Failure{core::FailureCode::InvalidInput, operation}};
}
template <class T>
Reply reply(core::Result<T> result) {
    if (!result) {
        return Reply{result.error()};
    }
    if constexpr (std::is_void_v<T>) {
        return Reply{py::none()};
    } else {
        return Reply{py::cast(std::move(*result))};
    }
}
core::Result<std::uint64_t> integer(py::handle object, core::Operation operation) {
    if (!PyLong_CheckExact(object.ptr())) {
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    }
    const auto value = PyLong_AsUnsignedLongLong(object.ptr());
    if (PyErr_Occurred() != nullptr) {
        PyErr_Clear();
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    }
    return static_cast<std::uint64_t>(value);
}
template <class T>
void require_owner(const std::shared_ptr<T>& owner) {
    if (!owner) {
        throw py::type_error("A non-None owning instance is required");
    }
}

struct PlanOwner {
    explicit PlanOwner(std::uint64_t limit) : plan(limit) {}
    std::mutex mutex;
    core::ResourcePlan plan;
};
template <class T>
struct MutableOwner {
    explicit MutableOwner(T input) : value(std::move(input)) {}
    std::mutex mutex;
    T value;
};
using TrialOwner = MutableOwner<core::TrialState>;
using PacketOwner = MutableOwner<core::QualificationPacket>;

// Never block holding the GIL against a worker that must reacquire it. Every
// mutable getter uses the same nonblocking locks as processing operations.
class OwnerLocks {
public:
    explicit OwnerLocks(std::initializer_list<std::mutex*> mutexes) {
        std::array<std::mutex*, 3> ordered{};
        if (mutexes.size() > ordered.size()) {
            return;
        }
        // An explicitly bounded insertion sort avoids libstdc++ introsort's
        // large-range paths for this tiny prefix, without suppressing warnings.
        std::size_t count = 0;
        for (auto* mutex : mutexes) {
            if (mutex == nullptr) {
                return;
            }
            auto position = count;
            while (position > 0 && std::less<std::mutex*>{}(mutex, ordered[position - 1])) {
                ordered[position] = ordered[position - 1];
                --position;
            }
            ordered[position] = mutex;
            ++count;
        }
        for (std::size_t index = 0; index < count; ++index) {
            if (index > 0 && ordered[index] == ordered[index - 1]) {
                continue;
            }
            locks_[size_] = std::unique_lock<std::mutex>{*ordered[index], std::try_to_lock};
            if (!locks_[size_].owns_lock()) {
                // Release the successful prefix immediately on partial failure.
                // No caller waits with the GIL held, even for duplicate owners.
                while (size_ > 0) {
                    locks_[--size_].unlock();
                }
                return;
            }
            ++size_;
        }
        acquired_ = true;
    }
    [[nodiscard]] bool acquired() const noexcept {
        return acquired_;
    }

private:
    std::array<std::unique_lock<std::mutex>, 3> locks_;
    std::size_t size_ = 0;
    bool acquired_ = false;
};
template <class T>
struct DenseOwner {
    explicit DenseOwner(numerics::DenseMatrix<T> input) : matrix(std::move(input)) {}
    numerics::DenseMatrix<T> matrix;
    T empty{};
};
struct FrozenBytes {
    explicit FrozenBytes(core::PlannedBuffer<std::byte> input) : bytes(std::move(input)) {}
    core::PlannedBuffer<std::byte> bytes;
    std::byte empty{};
};
struct ModelOwner {
    explicit ModelOwner(const core::Model& input) : model(input) {}
    core::Model model;
    std::byte empty{};
};
struct AcceptedOwner {
    std::optional<core::AcceptedState> state;
    std::byte empty{};
};
template <class Owner>
py::array byte_view(const std::shared_ptr<Owner>& owner, std::span<const std::byte> bytes,
                    const std::byte* empty) {
    // NumPy treats nullptr as a request to allocate; retain an owner for empty
    // arrays too, using a stable owner-local sentinel that is never dereferenced.
    const auto* data = bytes.empty() ? empty : bytes.data();
    auto array = py::array{py::dtype::of<std::uint8_t>(),
                           {static_cast<py::ssize_t>(bytes.size())},
                           {py::ssize_t{1}},
                           data,
                           py::cast(owner)};
    array.attr("setflags")(false);
    return array;
}
core::Result<std::span<const std::byte>> byte_input(const py::array& array,
                                                    core::Operation operation) {
    if (!array.dtype().is(py::dtype::of<std::uint8_t>()) || array.ndim() != 1 ||
        (array.flags() & py::array::c_style) == 0) {
        return std::unexpected(core::Failure{core::FailureCode::InvalidInput, operation});
    }
    return std::span<const std::byte>{static_cast<const std::byte*>(array.data()),
                                      static_cast<std::size_t>(array.size())};
}
template <class T>
void bind_dense(py::module_& module, const char* name) {
    using Owner = DenseOwner<T>;
    py::class_<Owner, std::shared_ptr<Owner>>(
        module, name, "Immutable owning column-major matrix; no implicit dtype/layout conversions.")
        .def_static(
            "copy",
            [](const py::array& input, const std::shared_ptr<PlanOwner>& plan) {
                return boundary(core::Operation::PrepareModel, [&]() -> Reply {
                    if (!input.dtype().is(py::dtype::of<T>()) || input.ndim() != 2 ||
                        (input.flags() & py::array::f_style) == 0 ||
                        std::bit_cast<std::uintptr_t>(input.data()) % alignof(T) != 0) {
                        return invalid(core::Operation::PrepareModel);
                    }
                    OwnerLocks const locks{&plan->mutex};
                    if (!locks.acquired()) {
                        return invalid(core::Operation::PrepareModel);
                    }
                    // Snapshot while holding the GIL before exposing any NoGil work.
                    auto copied = numerics::DenseMatrix<T>::copy(
                        {input.shape(0), input.shape(1), input.shape(0)},
                        {static_cast<const T*>(input.data()),
                         static_cast<std::size_t>(input.size())},
                        plan->plan, core::Operation::PrepareModel);
                    if (!copied) {
                        return Reply{copied.error()};
                    }
                    return Reply{py::cast(std::make_shared<Owner>(std::move(*copied)))};
                });
            },
            py::arg("array").noconvert(), py::arg("plan").none(false),
            "Snapshot an aligned native-dtype Fortran-contiguous 2D array under the GIL. Returns "
            "Result.")
        .def(
            "view",
            [](const std::shared_ptr<Owner>& owner) {
                require_owner(owner);
                const auto shape = owner->matrix.shape();
                const auto values = std::as_const(owner->matrix).values();
                const auto* data = values.empty() ? &owner->empty : values.data();
                auto array = py::array{
                    py::dtype::of<T>(),
                    {static_cast<py::ssize_t>(shape.rows), static_cast<py::ssize_t>(shape.columns)},
                    {static_cast<py::ssize_t>(sizeof(T)),
                     static_cast<py::ssize_t>(sizeof(T)) * shape.leading_dimension},
                    data,
                    py::cast(owner)};
                array.attr("setflags")(false);
                return array;
            },
            "Return a read-only zero-copy NumPy view retaining this owner, including empty arrays.")
        .def(
            "validate",
            [](const Owner& owner) {
                core::Result<void> result;
                {
                    py::gil_scoped_release const release;
                    result = owner.matrix.validate(core::Operation::PrepareModel);
                }
                return reply(std::move(result));
            },
            "Check stored components for finiteness with the GIL released. Returns Result.");
}

void bind_state(py::module_& module) {
    py::class_<core::ObservationRequirement>(module, "ObservationRequirement",
                                             "Observation, scope, error budget and arithmetic "
                                             "representation, in a caller-defined common scale.")
        .def_static(
            "create",
            [](py::handle observation, py::handle scope, double budget, core::Precision precision) {
                return boundary(core::Operation::PrepareModel, [&]() -> Reply {
                    auto obs = integer(observation, core::Operation::PrepareModel);
                    auto domain = integer(scope, core::Operation::PrepareModel);
                    if (!obs) {
                        return Reply{obs.error()};
                    }
                    if (!domain) {
                        return Reply{domain.error()};
                    }
                    auto certificate = core::ErrorCertificate::unevaluated(
                        {core::ObservationId{*obs}, core::ScopeId{*domain},
                         core::NumericalGeneration{0}, precision},
                        budget);
                    if (!certificate) {
                        return Reply{certificate.error()};
                    }
                    return Reply{py::cast(core::ObservationRequirement{
                        core::ObservationId{*obs}, core::ScopeId{*domain}, budget, precision})};
                });
            },
            "Validate identities and a nonnegative finite error budget. Returns Result.")
        .def_property_readonly(
            "observation",
            [](const core::ObservationRequirement& value) { return value.observation.value(); },
            "Identity of the observed quantity, preserved as a full uint64 integer.")
        .def_property_readonly(
            "scope", [](const core::ObservationRequirement& value) { return value.scope.value(); },
            "Identity of the evaluation domain, preserved without numeric conversion.")
        .def_readonly("budget", &core::ObservationRequirement::budget,
                      "Nonnegative finite error budget in the evaluator-defined common scale.")
        .def_readonly("representation", &core::ObservationRequirement::representation,
                      "Arithmetic representation required of the evaluated field.");
    py::class_<ModelOwner, std::shared_ptr<ModelOwner>>(
        module, "Model",
        "Immutable owning canonical model definition; this protocol does not define physical model "
        "serialization.")
        .def_static(
            "create",
            [](const py::array& definition, py::handle generation,
               const core::AcceptancePlan& requirements, const std::shared_ptr<PlanOwner>& plan) {
                return boundary(core::Operation::PrepareModel, [&]() -> Reply {
                    auto bytes = byte_input(definition, core::Operation::PrepareModel);
                    auto number = integer(generation, core::Operation::PrepareModel);
                    if (!bytes) {
                        return Reply{bytes.error()};
                    }
                    if (!number) {
                        return Reply{number.error()};
                    }
                    OwnerLocks const locks{&plan->mutex};
                    if (!locks.acquired()) {
                        return invalid(core::Operation::PrepareModel);
                    }
                    // Canonical bytes are copied with the GIL held. They are opaque
                    // uint8 protocol bytes, never a reinterpretation of float fields.
                    auto model = core::Model::create(*bytes, core::ModelGeneration{*number},
                                                     requirements, plan->plan);
                    if (!model) {
                        return Reply{model.error()};
                    }
                    return Reply{py::cast(std::make_shared<ModelOwner>(std::move(*model)))};
                });
            },
            py::arg("definition").noconvert(), py::arg("generation"), py::arg("requirements"),
            py::arg("plan").none(false),
            "Snapshot contiguous uint8 definition and four residual/geometry/time/precision "
            "requirements. Returns Result.")
        .def(
            "definition",
            [](const std::shared_ptr<ModelOwner>& owner) {
                require_owner(owner);
                return byte_view(owner, owner->model.definition(), &owner->empty);
            },
            "Return a read-only uint8 view retaining the immutable model owner.")
        .def_property_readonly(
            "generation", [](const ModelOwner& owner) { return owner.model.generation().value(); },
            "Immutable model generation supplied at preparation.")
        .def_property_readonly(
            "fingerprint",
            [](const ModelOwner& owner) {
                const auto fingerprint = owner.model.fingerprint();
                const auto bytes = std::bit_cast<std::array<char, 32>>(fingerprint.bytes);
                return py::bytes{bytes.data(), bytes.size()};
            },
            "32-byte SHA-256 identity of actual definition bytes and acceptance requirements.");
    py::class_<core::EvaluationTarget>(module, "EvaluationTarget",
                                       "Immutable candidate stamp; capture BEFORE evaluation and "
                                       "carry it with the computed bounds.")
        .def_property_readonly(
            "trial", [](const core::EvaluationTarget& target) { return target.trial.value(); },
            "Identity of the candidate captured for evaluation.")
        .def_property_readonly(
            "numerical_generation",
            [](const core::EvaluationTarget& target) {
                return target.numerical_generation.value();
            },
            "Numerical generation captured before evaluation.")
        .def_readonly(
            "revision", &core::EvaluationTarget::revision,
            "Revision of the candidate payload, incremented on every successful replacement.")
        .def_property_readonly(
            "parent",
            [](const core::EvaluationTarget& target) -> py::object {
                if (target.parent) {
                    return py::int_(target.parent->id.value());
                }
                return py::none();
            },
            "Accepted parent identity, or None for an initial candidate.");
    py::class_<TrialOwner, std::shared_ptr<TrialOwner>>(
        module, "TrialState",
        "Owned mutable unaccepted candidate. Concurrent operations/getters fail immediately.")
        .def_static(
            "initial",
            [](const std::shared_ptr<ModelOwner>& model, py::handle id, py::handle generation,
               const py::array& payload, const std::shared_ptr<PlanOwner>& plan) {
                return boundary(core::Operation::Initialize, [&]() -> Reply {
                    auto number = integer(id, core::Operation::Initialize);
                    auto gen = integer(generation, core::Operation::Initialize);
                    auto bytes = byte_input(payload, core::Operation::Initialize);
                    if (!number) {
                        return Reply{number.error()};
                    }
                    if (!gen) {
                        return Reply{gen.error()};
                    }
                    if (!bytes) {
                        return Reply{bytes.error()};
                    }
                    OwnerLocks const locks{&plan->mutex};
                    if (!locks.acquired()) {
                        return invalid(core::Operation::Initialize);
                    }
                    auto trial = core::TrialState::initial(model->model, core::TrialId{*number},
                                                           core::NumericalGeneration{*gen}, *bytes,
                                                           plan->plan);
                    if (!trial) {
                        return Reply{trial.error()};
                    }
                    return Reply{py::cast(std::make_shared<TrialOwner>(std::move(*trial)))};
                });
            },
            py::arg("model").none(false), py::arg("id"), py::arg("generation"),
            py::arg("payload").noconvert(), py::arg("plan").none(false),
            "Snapshot an initial canonical uint8 candidate under the GIL. Returns Result; no "
            "acceptance is implied.")
        .def_static(
            "begin",
            [](const std::shared_ptr<AcceptedOwner>& parent, py::handle id, py::handle generation,
               const py::array& payload, const std::shared_ptr<PlanOwner>& plan) {
                return boundary(core::Operation::AdvanceTrial, [&]() -> Reply {
                    auto number = integer(id, core::Operation::AdvanceTrial);
                    auto gen = integer(generation, core::Operation::AdvanceTrial);
                    auto bytes = byte_input(payload, core::Operation::AdvanceTrial);
                    if (!number) {
                        return Reply{number.error()};
                    }
                    if (!gen) {
                        return Reply{gen.error()};
                    }
                    if (!bytes) {
                        return Reply{bytes.error()};
                    }
                    OwnerLocks const locks{&plan->mutex};
                    if (!locks.acquired()) {
                        return invalid(core::Operation::AdvanceTrial);
                    }
                    if (!parent->state) {
                        return Reply{core::Failure{core::FailureCode::StaleGeneration,
                                                   core::Operation::AdvanceTrial}};
                    }
                    auto trial = core::TrialState::begin(*parent->state, core::TrialId{*number},
                                                         core::NumericalGeneration{*gen}, *bytes,
                                                         plan->plan);
                    if (!trial) {
                        return Reply{trial.error()};
                    }
                    return Reply{py::cast(std::make_shared<TrialOwner>(std::move(*trial)))};
                });
            },
            py::arg("parent").none(false), py::arg("id"), py::arg("generation"),
            py::arg("payload").noconvert(), py::arg("plan").none(false),
            "Snapshot a newer candidate without mutating its accepted parent. Returns Result.")
        .def(
            "target",
            [](const std::shared_ptr<TrialOwner>& owner) {
                require_owner(owner);
                return boundary(core::Operation::EvaluateStep, [&]() -> Reply {
                    OwnerLocks const locks{&owner->mutex};
                    if (!locks.acquired()) {
                        return invalid(core::Operation::EvaluateStep);
                    }
                    if (!owner->value.active()) {
                        return Reply{core::Failure{core::FailureCode::StaleGeneration,
                                                   core::Operation::EvaluateStep}};
                    }
                    const auto target = [&]() {
                        py::gil_scoped_release const release;
                        return owner->value.target();
                    }();
                    return Reply{py::cast(target)};
                });
            },
            "Hash the owned candidate with the GIL released. Capture this Result's target before "
            "computing bounds.")
        .def(
            "active",
            [](const std::shared_ptr<TrialOwner>& owner) {
                require_owner(owner);
                OwnerLocks const locks{&owner->mutex};
                if (!locks.acquired()) {
                    return invalid(core::Operation::EvaluateStep);
                }
                return Reply{py::bool_(owner->value.active())};
            },
            "Return Result containing whether this candidate remains usable; concurrent access is "
            "rejected.")
        .def(
            "payload_copy",
            [](const std::shared_ptr<TrialOwner>& owner, const std::shared_ptr<PlanOwner>& plan) {
                require_owner(owner);
                return boundary(core::Operation::EvaluateStep, [&]() -> Reply {
                    OwnerLocks const locks{&owner->mutex, &plan->mutex};
                    if (!locks.acquired()) {
                        return invalid(core::Operation::EvaluateStep);
                    }
                    const auto source = owner->value.payload();
                    auto storage = plan->plan.allocate<std::byte>(
                        static_cast<core::Index>(source.size()), core::Operation::EvaluateStep);
                    if (!storage) {
                        return Reply{storage.error()};
                    }
                    std::copy(source.begin(), source.end(), storage->values().begin());
                    return Reply{py::cast(std::make_shared<FrozenBytes>(std::move(*storage)))};
                });
            },
            py::arg("plan").none(false),
            "Return Result containing an independent ByteSnapshot charged to plan; later "
            "replacement/cancel/commit "
            "cannot invalidate it.")
        .def(
            "replace_payload",
            [](const std::shared_ptr<TrialOwner>& owner, const py::array& payload,
               const std::shared_ptr<PlanOwner>& plan) -> py::object {
                // Allocate the Python success carrier BEFORE modifying native state.
                if (!owner) {
                    throw py::type_error("TrialState is required");
                }
                auto prepared = py::cast(Reply{py::none()});
                auto& destination = prepared.cast<Reply&>();
                auto result = boundary(core::Operation::AdvanceTrial, [&]() -> Reply {
                    auto bytes = byte_input(payload, core::Operation::AdvanceTrial);
                    if (!bytes) {
                        return Reply{bytes.error()};
                    }
                    OwnerLocks const locks{&owner->mutex, &plan->mutex};
                    if (!locks.acquired()) {
                        return invalid(core::Operation::AdvanceTrial);
                    }
                    return reply(owner->value.replace_payload(*bytes, plan->plan));
                });
                destination.outcome = std::move(result.outcome);
                return prepared;
            },
            py::arg("payload").noconvert(), py::arg("plan").none(false),
            "Transactionally snapshot contiguous uint8 bytes under the GIL. Success invalidates "
            "earlier reports/packets.")
        .def(
            "cancel",
            [](const std::shared_ptr<TrialOwner>& owner) -> py::object {
                if (!owner) {
                    throw py::type_error("TrialState is required");
                }
                auto prepared = py::cast(Reply{
                    core::Failure{core::FailureCode::Cancelled, core::Operation::AdvanceTrial}});
                auto& destination = prepared.cast<Reply&>();
                OwnerLocks const locks{&owner->mutex};
                destination.outcome = locks.acquired()
                                          ? owner->value.cancel()
                                          : core::Failure{core::FailureCode::InvalidInput,
                                                          core::Operation::AdvanceTrial};
                return prepared;
            },
            "Invalidate this candidate and return Cancelled, preserving its accepted parent. Reuse "
            "returns StaleGeneration.");
    py::class_<core::EvaluationReport>(
        module, "EvaluationReport",
        "Validated reported bounds; does not compute equations or establish physical correctness.")
        .def_static(
            "create",
            [](const std::shared_ptr<ModelOwner>& model,
               const core::EvaluationTarget& evaluated_target,
               const std::array<std::optional<double>, 4>& bounds) {
                return boundary(core::Operation::CertifyError, [&]() -> Reply {
                    const auto assess = [&](std::size_t index) {
                        const auto requirement = model->model.requirements()[index];
                        const core::CertificateContext context{
                            requirement.observation, requirement.scope,
                            evaluated_target.numerical_generation, requirement.representation};
                        return bounds[index]
                                   ? core::ErrorCertificate::assess(context, requirement.budget,
                                                                    *bounds[index])
                                   : core::ErrorCertificate::unevaluated(context,
                                                                         requirement.budget);
                    };
                    auto a = assess(0);
                    auto b = assess(1);
                    auto c = assess(2);
                    auto d = assess(3);
                    if (!a) {
                        return Reply{a.error()};
                    }
                    if (!b) {
                        return Reply{b.error()};
                    }
                    if (!c) {
                        return Reply{c.error()};
                    }
                    if (!d) {
                        return Reply{d.error()};
                    }
                    return reply(core::EvaluationReport::create(model->model, evaluated_target,
                                                                {*a, *b, *c, *d}));
                });
            },
            py::arg("model").none(false), py::arg("evaluated_target"), py::arg("bounds"),
            "Record four caller-computed conservative bounds (None means unevaluated) against the "
            "target captured BEFORE evaluation.");
    py::class_<PacketOwner, std::shared_ptr<PacketOwner>>(
        module, "QualificationPacket",
        "Single-use qualification bound to one candidate revision and accepted parent.")
        .def(
            "consumed",
            [](const std::shared_ptr<PacketOwner>& owner) {
                require_owner(owner);
                OwnerLocks const locks{&owner->mutex};
                if (!locks.acquired()) {
                    return invalid(core::Operation::Commit);
                }
                return Reply{py::bool_(owner->value.consumed())};
            },
            "Return Result containing consumption status; concurrent access is rejected.");
    py::class_<AcceptedOwner, std::shared_ptr<AcceptedOwner>>(
        module, "AcceptedState",
        "Immutable accepted protocol state owning its payload independently of trial and plan.")
        .def_property_readonly(
            "id",
            [](const AcceptedOwner& owner) {
                if (!owner.state) {
                    throw py::value_error("Accepted state is not initialized");
                }
                return owner.state->id().value();
            },
            "Identity of the immutable accepted state.")
        .def_property_readonly(
            "numerical_generation",
            [](const AcceptedOwner& owner) {
                if (!owner.state) {
                    throw py::value_error("Accepted state is not initialized");
                }
                return owner.state->numerical_generation().value();
            },
            "Numerical generation of the accepted payload.")
        .def(
            "payload",
            [](const std::shared_ptr<AcceptedOwner>& owner) {
                require_owner(owner);
                if (!owner->state) {
                    throw py::value_error("Accepted state is not initialized");
                }
                return byte_view(owner, owner->state->payload(), &owner->empty);
            },
            "Return a read-only uint8 view retaining this immutable accepted state.");
    module.def(
        "qualify",
        [](const std::shared_ptr<TrialOwner>& trial, const core::EvaluationReport& report) {
            return boundary(core::Operation::Commit, [&]() -> Reply {
                OwnerLocks const locks{&trial->mutex};
                if (!locks.acquired()) {
                    return invalid(core::Operation::Commit);
                }
                auto packet = [&]() {
                    py::gil_scoped_release const release;
                    return core::qualify(trial->value, report);
                }();
                if (!packet) {
                    return Reply{packet.error()};
                }
                return Reply{py::cast(std::make_shared<PacketOwner>(std::move(*packet)))};
            });
        },
        py::arg("trial").none(false), py::arg("report"),
        "Qualify only a matching current candidate with all four reported bounds achieved. Returns "
        "Result.");
    const auto commit =
        [](const std::shared_ptr<AcceptedOwner>& parent, const std::shared_ptr<TrialOwner>& trial,
           const std::shared_ptr<PacketOwner>& packet, py::handle initial_id) -> py::object {
        // Prepare ALL owning/Python carriers before consuming trial or packet.
        // AcceptedState's move construction and optional emplacement are noexcept.
        auto accepted = std::make_shared<AcceptedOwner>();
        auto prepared = py::cast(Reply{py::cast(accepted)});
        auto& destination = prepared.cast<Reply&>();
        auto number = integer(initial_id, core::Operation::Commit);
        if (!number) {
            destination.outcome = number.error();
            return prepared;
        }
        OwnerLocks const locks{&trial->mutex, &packet->mutex};
        if (!locks.acquired()) {
            destination.outcome =
                core::Failure{core::FailureCode::InvalidInput, core::Operation::Commit};
            return prepared;
        }
        auto result = [&]() {
            py::gil_scoped_release const release;
            if (parent) {
                if (!parent->state) {
                    return core::Result<core::AcceptedState>{std::unexpected(core::Failure{
                        core::FailureCode::StaleGeneration, core::Operation::Commit})};
                }
                return core::commit(*parent->state, trial->value, packet->value);
            }
            return core::commit_initial(trial->value, packet->value,
                                        core::AcceptedStateId{*number});
        }();
        if (result) {
            accepted->state.emplace(std::move(*result));
        } else {
            destination.outcome = result.error();
        }
        return prepared;
    };
    module.def(
        "commit_initial",
        [commit](const std::shared_ptr<TrialOwner>& trial,
                 const std::shared_ptr<PacketOwner>& packet,
                 py::handle id) { return commit(nullptr, trial, packet, id); },
        py::arg("trial").none(false), py::arg("packet").none(false), py::arg("id"),
        "Accept a qualified initial candidate exactly once; preallocated result publication "
        "preserves ownership on failure.");
    module.def(
        "commit",
        [commit](const std::shared_ptr<AcceptedOwner>& parent,
                 const std::shared_ptr<TrialOwner>& trial,
                 const std::shared_ptr<PacketOwner>& packet) {
            return commit(parent, trial, packet, py::int_(0));
        },
        py::arg("parent").none(false), py::arg("trial").none(false), py::arg("packet").none(false),
        "Commit a current qualified child while preserving its immutable parent; consumes "
        "trial/packet only on success.");
}
} // namespace

PYBIND11_MODULE(_uvlas, module) {
    module.doc() = "Low-level owning product storage/state protocols; physical analysis and file "
                   "APIs are not implemented.";
    py::class_<Reply>(module, "Result",
                      "An owned success value or Failure. Inspect ok before value.")
        .def_property_readonly(
            "ok",
            [](const Reply& value) { return std::holds_alternative<py::object>(value.outcome); },
            "True when value contains an owned success; false when failure contains a diagnostic.")
        .def_property_readonly(
            "value",
            [](const Reply& value) {
                if (const auto* success = std::get_if<py::object>(&value.outcome)) {
                    return *success;
                }
                throw py::value_error("Result has no success value; inspect failure");
            },
            "Owned success object; raises ValueError when ok is false.")
        .def_property_readonly(
            "failure",
            [](const Reply& value) -> py::object {
                if (const auto* failure = std::get_if<core::Failure>(&value.outcome)) {
                    return py::cast(*failure);
                }
                return py::none();
            },
            "Owned Failure when ok is false, otherwise None.");
    py::enum_<core::FailureCode>(
        module, "FailureCode",
        "Distinct product failure categories; names are stable, ordinals are not.")
        .value("InvalidInput", core::FailureCode::InvalidInput)
        .value("Unsupported", core::FailureCode::Unsupported)
        .value("GeometryInvalid", core::FailureCode::GeometryInvalid)
        .value("RankUnresolved", core::FailureCode::RankUnresolved)
        .value("ChartSingular", core::FailureCode::ChartSingular)
        .value("Nonfinite", core::FailureCode::Nonfinite)
        .value("FactorizationFailed", core::FailureCode::FactorizationFailed)
        .value("LinearUnconverged", core::FailureCode::LinearUnconverged)
        .value("NonlinearUnconverged", core::FailureCode::NonlinearUnconverged)
        .value("EigenCoverageUnresolved", core::FailureCode::EigenCoverageUnresolved)
        .value("PrecisionUnsupported", core::FailureCode::PrecisionUnsupported)
        .value("PrecisionUnresolved", core::FailureCode::PrecisionUnresolved)
        .value("ResourceExceeded", core::FailureCode::ResourceExceeded)
        .value("Cancelled", core::FailureCode::Cancelled)
        .value("StaleGeneration", core::FailureCode::StaleGeneration)
        .value("CheckpointCorrupt", core::FailureCode::CheckpointCorrupt)
        .value("CheckpointIncompatible", core::FailureCode::CheckpointIncompatible);
    py::enum_<core::Operation>(module, "Operation",
                               "Product entry operation responsible for a failure.")
        .value("PrepareModel", core::Operation::PrepareModel)
        .value("EvaluateStructure", core::Operation::EvaluateStructure)
        .value("EvaluateAero", core::Operation::EvaluateAero)
        .value("Initialize", core::Operation::Initialize)
        .value("EvaluateStep", core::Operation::EvaluateStep)
        .value("SolveLinear", core::Operation::SolveLinear)
        .value("AdvanceTrial", core::Operation::AdvanceTrial)
        .value("Commit", core::Operation::Commit)
        .value("LinearizeEquilibrium", core::Operation::LinearizeEquilibrium)
        .value("SearchFlutter", core::Operation::SearchFlutter)
        .value("Recover", core::Operation::Recover)
        .value("CertifyError", core::Operation::CertifyError)
        .value("ExtractEigen", core::Operation::ExtractEigen)
        .value("SaveCheckpoint", core::Operation::SaveCheckpoint)
        .value("WriteResults", core::Operation::WriteResults);
    py::enum_<core::ResourceKind>(module, "ResourceKind",
                                  "Unit shared by required and observed resource counts.")
        .value("Bytes", core::ResourceKind::Bytes)
        .value("Iterations", core::ResourceKind::Iterations)
        .value("OperatorApplications", core::ResourceKind::OperatorApplications)
        .value("Nanoseconds", core::ResourceKind::Nanoseconds);
    py::class_<core::Failure>(module, "Failure",
                              "Owned failure with category, operation and diagnostic context.")
        .def_readonly("code", &core::Failure::code,
                      "Named failure category; inspect names rather than enum ordinals.")
        .def_readonly("operation", &core::Failure::operation,
                      "Product entry responsible for the failure.")
        .def_property_readonly(
            "details",
            [](const core::Failure& failure) {
                py::dict result;
                const auto& details = failure.details;
                const auto put = [&result](const char* name, const auto& value) {
                    if (value) {
                        result[name] = py::int_(value->value());
                    } else {
                        result[name] = py::none();
                    }
                };
                put("model", details.model);
                put("parent", details.parent);
                put("trial", details.trial);
                put("numerical_generation", details.numerical_generation);
                put("component", details.component);
                put("last_accepted", details.last_accepted);
                if (details.resources) {
                    py::dict resource;
                    resource["kind"] = py::cast(details.resources->kind);
                    resource["required"] = details.resources->required;
                    resource["observed"] = details.resources->observed;
                    result["resources"] = resource;
                } else {
                    result["resources"] = py::none();
                }
                if (details.problem) {
                    result["problem"] = std::visit(
                        [](const auto& id) { return py::int_(id.value()); }, *details.problem);
                    constexpr std::array domains{"dof",   "panel",       "wake_vertex", "wake_face",
                                                 "layer", "observation", "scope"};
                    result["problem_domain"] = domains[details.problem->index()];
                } else {
                    result["problem"] = py::none();
                    result["problem_domain"] = py::none();
                }
                return result;
            },
            "Owned identity/resource dictionary. Unknown fields are None; resource kind and "
            "problem domain retain their units/type.");
    py::enum_<core::Precision>(module, "Precision")
        .value("Binary64", core::Precision::Binary64)
        .value("TwoWord", core::Precision::TwoWord);
    py::class_<PlanOwner, std::shared_ptr<PlanOwner>>(
        module, "ResourcePlan",
        "Cumulative native payload budget; concurrent use fails immediately. Python/allocator "
        "overhead is excluded.")
        .def_static(
            "create",
            [](py::handle limit) {
                return boundary(core::Operation::PrepareModel, [&]() -> Reply {
                    auto parsed = integer(limit, core::Operation::PrepareModel);
                    if (!parsed) {
                        return Reply{parsed.error()};
                    }
                    return Reply{py::cast(std::make_shared<PlanOwner>(*parsed))};
                });
            },
            "Create a budget from an exact unsigned 64-bit Python integer. Returns Result.")
        .def(
            "usage",
            [](const std::shared_ptr<PlanOwner>& owner) {
                require_owner(owner);
                OwnerLocks const locks{&owner->mutex};
                if (!locks.acquired()) {
                    return invalid(core::Operation::PrepareModel);
                }
                return Reply{py::make_tuple(owner->plan.byte_limit(), owner->plan.allocated_bytes(),
                                            owner->plan.remaining_bytes())};
            },
            "Return Result containing (limit, allocated, remaining) bytes; concurrent access is "
            "rejected.");
    bind_dense<double>(module, "DenseReal");
    bind_dense<std::complex<double>>(module, "DenseComplex");
    py::class_<FrozenBytes, std::shared_ptr<FrozenBytes>>(
        module, "ByteSnapshot",
        "Independent immutable canonical protocol bytes, without physical field interpretation.")
        .def(
            "view",
            [](const std::shared_ptr<FrozenBytes>& owner) {
                require_owner(owner);
                return byte_view(owner, owner->bytes.values(), &owner->empty);
            },
            "Return a read-only uint8 view retaining the independent snapshot owner.");
    bind_state(module);
}
