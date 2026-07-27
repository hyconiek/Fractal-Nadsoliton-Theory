#!/usr/bin/env python3
"""Regression tests for the FIN P240/P241/P242 laboratory transfer package."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path
import csv
import hashlib
import json
import tempfile
import unittest

import numpy as np
from scipy.linalg import expm

import fin_lab_p240_optimal_tomography as p240
import fin_lab_p241_validator as p241
import fin_lab_p242_pipeline as p242


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def integer_counts(probabilities: np.ndarray, shots: int) -> np.ndarray:
    counts = np.zeros_like(probabilities, dtype=np.int64)
    for column in range(probabilities.shape[1]):
        raw = np.floor(probabilities[:, column] * shots).astype(np.int64)
        residual = shots - int(raw.sum())
        fractions = probabilities[:, column] * shots - raw
        for row in np.argsort(fractions)[::-1][:residual]:
            raw[row] += 1
        counts[:, column] = raw
    return counts


class TestProgram240(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.a, cls.w, cls.s = p240.strict_operator()
        cls.eigenvalues = np.linalg.eigvalsh(cls.a)

    def test_row_sum_reference(self):
        self.assertAlmostEqual(cls := self.s, p240.ROW_SUM_REFERENCE, places=14)
        self.assertGreater(cls, 0)

    def test_all_off_diagonal_weights_positive(self):
        mask = ~np.eye(12, dtype=bool)
        self.assertTrue(np.all(self.w[mask] > 0))

    def test_dirichlet_generator(self):
        self.assertLess(abs(self.eigenvalues[0]), 1e-12)
        self.assertGreater(self.eigenvalues[1], 0)
        self.assertLess(np.max(np.abs(self.a.sum(axis=1))), 1e-12)

    def test_exact_time_optimum(self):
        self.assertGreater(
            p240.scalar_inverse_log_amplification(0.8),
            p240.scalar_inverse_log_amplification(1.0),
        )
        self.assertGreater(
            p240.scalar_inverse_log_amplification(1.2),
            p240.scalar_inverse_log_amplification(1.0),
        )

    def test_transition_minimum_is_exp_minus_one(self):
        tau = 1.0 / self.eigenvalues[-1]
        transition = expm(-tau * self.a)
        self.assertAlmostEqual(
            float(np.linalg.eigvalsh(transition)[0]), np.exp(-1.0), places=13
        )

    def test_bernstein_radius_decreases(self):
        small = p240.bernstein_operator_radius(12, 10_000, 0.05)
        large = p240.bernstein_operator_radius(12, 50_000, 0.05)
        self.assertGreater(small, large)

    def test_design_lock_is_canonical(self):
        record = json.loads(p240.LOCK_PATH.read_text(encoding="utf-8"))
        self.assertEqual(
            record["canonical_core_sha256"],
            p240.canonical_digest(record["core"]),
        )
        self.assertFalse(record["core"]["external_validation_claimed"])


class TestProgram241(unittest.TestCase):
    def make_fixture(self, root: Path, resource_type: str) -> Path:
        bundle = root / f"{resource_type}_bundle"
        bundle.mkdir()
        now = datetime.now(timezone.utc)
        registered = now - timedelta(days=2)
        acquired = now - timedelta(days=1)

        if resource_type == "heat_process":
            fieldnames = sorted(p241.PROCESS_REQUIRED_COLUMNS)
            rows = []
            event_id = 0
            for subset, multiple in (("calibration", 1), ("holdout", 2)):
                for preparation in range(12):
                    for outcome in range(12):
                        rows.append(
                            {
                                "event_id": str(event_id),
                                "timestamp_utc": acquired.isoformat(),
                                "run_id": f"{subset}-{preparation}",
                                "subset": subset,
                                "preparation_id": str(preparation),
                                "outcome_id": str(outcome),
                                "evolution_multiple": str(multiple),
                                "intervention": "none",
                            }
                        )
                        event_id += 1
        else:
            fieldnames = sorted(p241.DOUBLE_SLIT_REQUIRED_COLUMNS)
            rows = []
            for event_id, configuration in enumerate(sorted(p241.DOUBLE_SLIT_CONFIGS)):
                rows.append(
                    {
                        "event_id": str(event_id),
                        "timestamp_utc": acquired.isoformat(),
                        "run_id": configuration,
                        "subset": "holdout" if event_id % 2 else "calibration",
                        "configuration": configuration,
                        "x_detector": str(event_id / 10),
                        "y_detector": "0.0",
                        "intervention": f"shutter:{configuration}",
                    }
                )
        with (bundle / "events.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        objects = {
            "calibration.json": {
                "clock_model": "test",
                "clock_unit": "s",
                "clock_uncertainty": "1e-6 s",
                "time_origin": "trigger",
                "detector_geometry": "fixture",
                "efficiency_calibration": "fixture",
                "dark_count_calibration": "fixture",
                "blur_or_psf_calibration": "fixture",
            },
            "controls.json": {
                "negative_controls": (
                    ["both_closed", "left_only", "right_only"]
                    if resource_type == "double_slit"
                    else [
                        "time_order_permutation",
                        "preparation_label_permutation",
                        "nearest_neighbour_C12",
                    ]
                )
            },
            "environment.json": {
                "background_record": "fixture",
                "temperature_or_not_applicable": "not applicable",
                "source_stability_record": "fixture",
                "apparatus_state_record": "fixture",
            },
            "preregistration.json": {
                "registered_at_utc": registered.isoformat(),
                "held_out_target": "P2=P1@P1",
                "frozen_analysis_hash": "fixture",
            },
        }
        for name, value in objects.items():
            (bundle / name).write_text(
                json.dumps(value, indent=2) + "\n", encoding="utf-8"
            )
        files = {
            name: file_hash(bundle / name)
            for name in sorted(p241.REQUIRED_FILES)
        }
        manifest = {
            "schema_version": p241.SCHEMA_VERSION,
            "resource_type": resource_type,
            "evidence_status": "external",
            "acquisition_started_at_utc": acquired.isoformat(),
            "provider": {
                "name": "Fixture Provider",
                "institution": "Fixture Lab",
                "contact": "fixture@example.invalid",
                "license": "CC BY 4.0",
            },
            "custody": {
                "provider_id": "provider",
                "registrar_id": "registrar",
                "analyst_id": "analyst",
            },
            "declaration": {
                "synthetic": False,
                "unit_test_fixture": True,
                "rendered_image_only": False,
            },
            "holdout": {
                "sealed_before_analysis": True,
                "committed_before_analysis": True,
            },
            "files": files,
        }
        (bundle / "bundle_manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        return bundle

    def test_template_is_not_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            target = Path(temporary) / "template"
            p241.emit_templates(target)
            self.assertTrue((target / "README.txt").is_file())
            manifest = json.loads(
                (target / "bundle_manifest.template.json").read_text()
            )
            self.assertEqual(manifest["evidence_status"], "EXTERNAL_REQUIRED")

    def test_process_fixture_cannot_pass_external_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = self.make_fixture(Path(temporary), "heat_process")
            result = p241.validate_bundle(bundle)
            self.assertFalse(result["core"]["program_241_gate_passed"])
            self.assertEqual(result["core"]["passed_fields"], 9)
            self.assertFalse(result["core"]["registrar_signature_verified"])

    def test_double_slit_fixture_cannot_unlock_semigroup(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = self.make_fixture(Path(temporary), "double_slit")
            result = p241.validate_bundle(bundle)
            self.assertFalse(result["core"]["program_242_semigroup_ready"])

    def test_hash_tampering_is_detected(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = self.make_fixture(Path(temporary), "heat_process")
            with (bundle / "events.csv").open("a", encoding="utf-8") as handle:
                handle.write("\n")
            result = p241.validate_bundle(bundle)
            hash_check = result["core"]["checks"][1]
            self.assertFalse(hash_check["passed"])

    def test_rendered_image_cannot_substitute(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = self.make_fixture(Path(temporary), "double_slit")
            manifest_path = bundle / "bundle_manifest.json"
            manifest = json.loads(manifest_path.read_text())
            manifest["declaration"]["rendered_image_only"] = True
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
            result = p241.validate_bundle(bundle)
            self.assertFalse(result["core"]["checks"][0]["passed"])


class TestProgram242(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        a, _, _ = p240.strict_operator()
        eigenvalues = np.linalg.eigvalsh(a)
        cls.tau = 1.0 / eigenvalues[-1]
        cls.p1 = expm(-cls.tau * a)
        cls.p2 = cls.p1 @ cls.p1

    def test_analysis_lock_canonical(self):
        record = p242.verify_analysis_lock()
        self.assertEqual(
            record["canonical_core_sha256"],
            p242.canonical_digest(record["core"]),
        )

    def test_exact_semigroup_not_falsified(self):
        counts1 = integer_counts(self.p1, 1_000_000)
        counts2 = integer_counts(self.p2, 1_000_000)
        result = p242.analyze_counts(counts1, counts2)
        self.assertEqual(
            result["primary"]["decision"],
            "NOT_FALSIFIED_AT_REGISTERED_LEVEL",
        )
        self.assertLess(result["primary"]["observed"], 1e-4)

    def test_wrong_holdout_is_falsified(self):
        counts1 = integer_counts(self.p1, 1_000_000)
        counts2 = integer_counts(np.eye(12), 1_000_000)
        result = p242.analyze_counts(counts1, counts2)
        self.assertEqual(
            result["primary"]["decision"],
            "FALSIFIED_AT_REGISTERED_LEVEL",
        )

    def test_primary_radius_decreases_with_shots(self):
        small = p242.registered_primary_radius(10_000, 10_000)
        large = p242.registered_primary_radius(50_000, 50_000)
        self.assertGreater(
            small["registered_tv_radius"], large["registered_tv_radius"]
        )

    def test_fingerprint_target(self):
        signature, diagnostic = p242.projective_fingerprint(self.p1)
        self.assertIsNotNone(signature)
        target = p240.projective_signature_from_generator(p240.strict_operator()[0])
        self.assertLess(np.max(np.abs(signature - target)), 1e-12)
        self.assertGreater(diagnostic["minimum_transition_eigenvalue"], 0)

    def test_execution_token_is_atomic(self):
        with tempfile.TemporaryDirectory() as temporary:
            ledger = Path(temporary)
            first = p242.consume_execution_token(ledger, "abc")
            self.assertTrue(first.is_file())
            with self.assertRaises(FileExistsError):
                p242.consume_execution_token(ledger, "abc")


if __name__ == "__main__":
    unittest.main()
