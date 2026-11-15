"""Exercise ownership and failure boundaries against the actual product extension."""

from __future__ import annotations

import gc
import importlib
import subprocess
import sys
import threading
import time
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(sys.argv.pop(1)).resolve()))
uvlas = importlib.import_module("_uvlas")


class OwnershipTests(unittest.TestCase):
    """IO4-02/N4-03: real storage/state processing, never a synthetic solver."""

    def success(self, result):
        self.assertTrue(result.ok, None if result.ok else result.failure.code)
        self.assertIsNone(result.failure)
        return result.value

    def failure(self, result, code):
        self.assertFalse(result.ok)
        self.assertEqual(result.failure.code, code)
        with self.assertRaises(ValueError):
            _ = result.value
        return result.failure

    def plan(self, limit=1_000_000):
        return self.success(uvlas.ResourcePlan.create(limit))

    def model(self, plan):
        requirements = [
            self.success(uvlas.ObservationRequirement.create(i, 10, 1.0, uvlas.Precision.Binary64))
            for i in range(4)
        ]
        return self.success(uvlas.Model.create(np.arange(5, dtype=np.uint8), 7, requirements, plan))

    def trial(self, model, plan, number=1):
        return self.success(
            uvlas.TrialState.initial(model, number, 20, np.arange(6, dtype=np.uint8), plan)
        )

    def packet(self, model, trial):
        target = self.success(trial.target())
        report = self.success(uvlas.EvaluationReport.create(model, target, [0.5] * 4))
        return self.success(uvlas.qualify(trial, report))

    def readonly(self, array):
        self.assertFalse(array.flags.writeable)
        self.assertIsNotNone(array.base)
        for derived in (array, array.view(), array.reshape(-1, order="A"), array[...]):
            with self.assertRaises(ValueError):
                derived.setflags(write=True)
            with self.assertRaises(ValueError):
                derived.flags.writeable = True
        self.assertTrue(memoryview(array).readonly)

    def test_dense_ownership_and_empty(self):
        for scalar, owner_type in (
            (np.float64, uvlas.DenseReal),
            (np.complex128, uvlas.DenseComplex),
        ):
            for shape in ((3, 4), (0, 3), (3, 0), (0, 0)):
                with self.subTest(scalar=scalar, shape=shape):
                    source = np.ones(shape, dtype=scalar, order="F")
                    owner = self.success(owner_type.copy(source, self.plan()))
                    source.fill(9)
                    view = owner.view()
                    self.readonly(view)
                    self.assertEqual(view.shape, shape)
                    self.assertTrue(owner.validate().ok)
                    del owner
                    gc.collect()
                    np.testing.assert_array_equal(view, np.ones(shape, dtype=scalar))

    def test_empty_canonical_owners(self):
        plan = self.plan(0)
        requirements = [
            self.success(uvlas.ObservationRequirement.create(i, 0, 0.0, uvlas.Precision.Binary64))
            for i in range(4)
        ]
        model = self.success(uvlas.Model.create(np.empty(0, dtype=np.uint8), 0, requirements, plan))
        trial = self.success(
            uvlas.TrialState.initial(model, 0, 0, np.empty(0, dtype=np.uint8), plan)
        )
        definition = model.definition()
        snapshot = self.success(trial.payload_copy(plan)).view()
        target = self.success(trial.target())
        report = self.success(uvlas.EvaluationReport.create(model, target, [0.0] * 4))
        packet = self.success(uvlas.qualify(trial, report))
        accepted = self.success(uvlas.commit_initial(trial, packet, 0))
        payload = accepted.payload()
        del accepted, model, trial, packet, report, plan
        gc.collect()
        for array in (definition, snapshot, payload):
            self.assertEqual(array.size, 0)
            self.readonly(array)

    def test_strict_numpy_admission(self):
        base = np.ones((4, 6), dtype=np.float64, order="F")
        unaligned = np.ndarray((2, 3), dtype=np.float64, buffer=bytearray(49), offset=1, order="F")
        for array in (
            base.astype(np.float32),
            np.ones((4, 6)),
            base[::2],
            base[:, ::2],
            base[:, ::-1],
            base.ravel(order="F"),
            base.astype(">f8"),
            unaligned,
        ):
            with self.subTest(dtype=array.dtype, shape=array.shape, strides=array.strides):
                self.failure(
                    uvlas.DenseReal.copy(array, self.plan()), uvlas.FailureCode.InvalidInput
                )
        with self.assertRaises(TypeError):
            uvlas.DenseReal.copy([[1.0]], self.plan())
        for value in (np.inf, -np.inf, np.nan):
            self.failure(
                uvlas.DenseReal.copy(np.array([[value]]), self.plan()), uvlas.FailureCode.Nonfinite
            )
        self.failure(
            uvlas.DenseComplex.copy(np.array([[complex(1, np.inf)]]), self.plan()),
            uvlas.FailureCode.Nonfinite,
        )

    def test_failure_and_operation_names(self):
        self.assertEqual(
            set(uvlas.FailureCode.__members__),
            {
                "InvalidInput",
                "Unsupported",
                "GeometryInvalid",
                "RankUnresolved",
                "ChartSingular",
                "Nonfinite",
                "FactorizationFailed",
                "LinearUnconverged",
                "NonlinearUnconverged",
                "EigenCoverageUnresolved",
                "PrecisionUnsupported",
                "PrecisionUnresolved",
                "ResourceExceeded",
                "Cancelled",
                "StaleGeneration",
                "CheckpointCorrupt",
                "CheckpointIncompatible",
            },
        )
        self.assertEqual(
            set(uvlas.Operation.__members__),
            {
                "PrepareModel",
                "EvaluateStructure",
                "EvaluateAero",
                "Initialize",
                "EvaluateStep",
                "SolveLinear",
                "AdvanceTrial",
                "Commit",
                "LinearizeEquilibrium",
                "SearchFlutter",
                "Recover",
                "CertifyError",
            },
        )

    def test_integer_and_resource_diagnostics(self):
        largest = 2**64 - 1
        plan = self.plan(largest)
        self.assertEqual(self.success(plan.usage()), (largest, 0, largest))
        for value in (-1, 2**64, 1.5, "1", True, np.uint64(1)):
            self.failure(uvlas.ResourcePlan.create(value), uvlas.FailureCode.InvalidInput)
        failure = self.failure(
            uvlas.DenseReal.copy(np.ones((2, 3), order="F"), self.plan(47)),
            uvlas.FailureCode.ResourceExceeded,
        )
        self.assertEqual(failure.details["resources"]["required"], 48)
        self.assertEqual(failure.details["resources"]["observed"], 0)
        self.assertIsNone(failure.details["parent"])

    def test_state_lifecycle_parent_and_snapshot(self):
        plan = self.plan()
        model = self.model(plan)
        definition = model.definition()
        self.readonly(definition)
        trial = self.trial(model, plan, 2**64 - 1)
        self.assertEqual(self.success(trial.target()).trial, 2**64 - 1)
        snapshot = self.success(trial.payload_copy(self.plan())).view()
        self.readonly(snapshot)
        packet = self.packet(model, trial)
        parent = self.success(uvlas.commit_initial(trial, packet, 2**53 + 1))
        self.assertEqual(parent.id, 2**53 + 1)
        self.assertTrue(self.success(packet.consumed()))
        self.failure(uvlas.commit_initial(trial, packet, 0), uvlas.FailureCode.StaleGeneration)
        parent_view = parent.payload()
        self.readonly(parent_view)
        child = self.success(
            uvlas.TrialState.begin(parent, 2, 21, np.array([99], dtype=np.uint8), plan)
        )
        self.success(child.replace_payload(np.array([88, 77], dtype=np.uint8), plan))
        cancelled = self.success(
            uvlas.TrialState.begin(parent, 3, 21, np.array([55], dtype=np.uint8), plan)
        )
        self.failure(cancelled.cancel(), uvlas.FailureCode.Cancelled)
        np.testing.assert_array_equal(parent_view, np.arange(6))
        child_packet = self.packet(model, child)
        accepted_child = self.success(uvlas.commit(parent, child, child_packet))
        self.assertEqual(accepted_child.id, parent.id + 1)
        np.testing.assert_array_equal(accepted_child.payload(), [88, 77])
        del model, trial, packet, parent, child, child_packet, accepted_child, plan
        gc.collect()
        np.testing.assert_array_equal(definition, np.arange(5))
        np.testing.assert_array_equal(parent_view, np.arange(6))
        np.testing.assert_array_equal(snapshot, np.arange(6))

    def test_stale_unevaluated_exceeded_cancel(self):
        plan = self.plan()
        model = self.model(plan)
        trial = self.trial(model, plan)
        target = self.success(trial.target())
        other = self.trial(model, plan, 2)
        report = self.success(uvlas.EvaluationReport.create(model, target, [0.5] * 4))
        self.failure(uvlas.qualify(other, report), uvlas.FailureCode.StaleGeneration)
        old_packet = self.success(uvlas.qualify(trial, report))
        self.success(trial.replace_payload(np.array([3], dtype=np.uint8), plan))
        self.failure(uvlas.qualify(trial, report), uvlas.FailureCode.StaleGeneration)
        self.failure(uvlas.commit_initial(trial, old_packet, 0), uvlas.FailureCode.StaleGeneration)
        current = self.success(trial.target())
        for category, code in enumerate(
            (
                uvlas.FailureCode.NonlinearUnconverged,
                uvlas.FailureCode.GeometryInvalid,
                uvlas.FailureCode.PrecisionUnresolved,
                uvlas.FailureCode.PrecisionUnresolved,
            )
        ):
            for missing in (False, True):
                bounds = [0.0] * 4
                bounds[category] = None if missing else 2.0
                bound_report = self.success(uvlas.EvaluationReport.create(model, current, bounds))
                expected = uvlas.FailureCode.PrecisionUnresolved if missing else code
                self.failure(uvlas.qualify(trial, bound_report), expected)
        self.failure(trial.payload_copy(self.plan(0)), uvlas.FailureCode.ResourceExceeded)

        for bounds, code in (
            ([None] * 4, uvlas.FailureCode.PrecisionUnresolved),
            ([2.0, 0.0, 0.0, 0.0], uvlas.FailureCode.NonlinearUnconverged),
        ):
            failed_report = self.success(uvlas.EvaluationReport.create(model, current, bounds))
            self.failure(uvlas.qualify(trial, failed_report), code)
        before = self.success(trial.payload_copy(self.plan())).view()
        self.failure(
            trial.replace_payload(np.arange(9, dtype=np.uint8), self.plan(0)),
            uvlas.FailureCode.ResourceExceeded,
        )
        self.assertEqual(self.success(trial.target()).revision, current.revision)
        np.testing.assert_array_equal(self.success(trial.payload_copy(self.plan())).view(), before)
        self.failure(trial.cancel(), uvlas.FailureCode.Cancelled)
        self.failure(trial.cancel(), uvlas.FailureCode.StaleGeneration)
        self.assertFalse(self.success(trial.active()))

    def test_none_is_rejected(self):
        plan = self.plan()
        model = self.model(plan)
        trial = self.trial(model, plan)
        packet = self.packet(model, trial)
        report = self.success(
            uvlas.EvaluationReport.create(model, self.success(trial.target()), [0.0] * 4)
        )
        array = np.arange(3, dtype=np.uint8)
        calls = (
            lambda: uvlas.DenseReal.copy(np.ones((1, 1)), None),
            lambda: uvlas.Model.create(array, 0, [], None),
            lambda: uvlas.TrialState.initial(None, 0, 0, array, plan),
            lambda: uvlas.TrialState.initial(model, 0, 0, array, None),
            lambda: uvlas.TrialState.begin(None, 0, 0, array, plan),
            lambda: trial.replace_payload(array, None),
            lambda: trial.payload_copy(None),
            lambda: uvlas.EvaluationReport.create(None, self.success(trial.target()), [0.0] * 4),
            lambda: uvlas.qualify(None, report),
            lambda: uvlas.commit_initial(None, packet, 0),
            lambda: uvlas.commit_initial(trial, None, 0),
            lambda: uvlas.commit(None, trial, packet),
        )
        for call in calls:
            with self.assertRaises(TypeError):
                call()

    def test_none_self_is_rejected_in_subprocess(self):
        # An unchecked shared_ptr self crashes the process, so isolate this
        # boundary negative control from the remainder of the ownership suite.
        script = """
import sys
sys.path.insert(0, sys.argv[1])
import _uvlas as u
methods = [u.DenseReal.view, u.DenseComplex.view, u.Model.definition,
           u.TrialState.target, u.TrialState.active, u.TrialState.cancel,
           u.QualificationPacket.consumed, u.AcceptedState.payload,
           u.ResourcePlan.usage, u.ByteSnapshot.view]
for method in methods:
    try:
        method(None)
    except TypeError:
        continue
    raise AssertionError(method)
for method, rest in [(u.TrialState.payload_copy, [None]),
                     (u.TrialState.replace_payload, [None, None])]:
    try:
        method(None, *rest)
    except TypeError:
        continue
    raise AssertionError(method)
"""
        completed = subprocess.run(
            [sys.executable, "-c", script, str(Path(uvlas.__file__).parent)],
            capture_output=True,
            text=True,
            check=False,
            timeout=30,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_owned_native_processing_releases_gil_and_rejects_overlap(self):
        # Hash enough real owned payload to overlap another Python thread. This
        # is a synchronization assertion, with no wall-time/performance threshold.
        plan = self.plan(40_000_010)
        model = self.model(plan)
        source = np.zeros(40_000_000, dtype=np.uint8)
        trial = self.success(uvlas.TrialState.initial(model, 1, 1, source, plan))
        packet = self.packet(model, trial)
        started = threading.Event()
        completed = threading.Event()
        results = []
        errors = []

        def hash_candidate():
            try:
                started.set()
                results.append(trial.target())
            except Exception as error:
                errors.append(error)
            finally:
                completed.set()

        worker = threading.Thread(target=hash_candidate)
        worker.start()
        self.assertTrue(started.wait(timeout=10))
        overlap = False
        deadline = time.monotonic() + 30
        while not completed.is_set() and time.monotonic() < deadline:
            result = uvlas.commit_initial(trial, packet, 0)
            if not result.ok:
                self.failure(result, uvlas.FailureCode.InvalidInput)
                overlap = True
                # A rejected multi-owner call must release any acquired packet lock.
                self.assertFalse(self.success(packet.consumed()))
                source.fill(9)
                break
        worker.join(timeout=30)
        self.assertFalse(worker.is_alive())
        if errors:
            raise errors[0]
        self.assertTrue(
            overlap, "Python did not execute while native target hashing held the owner lock"
        )
        self.success(results[0])
        np.testing.assert_array_equal(
            self.success(trial.payload_copy(self.plan(40_000_000))).view(), 0
        )
        accepted = self.success(uvlas.commit_initial(trial, packet, 0))
        self.assertTrue(self.success(packet.consumed()))
        np.testing.assert_array_equal(accepted.payload(), 0)


if __name__ == "__main__":
    unittest.main()
