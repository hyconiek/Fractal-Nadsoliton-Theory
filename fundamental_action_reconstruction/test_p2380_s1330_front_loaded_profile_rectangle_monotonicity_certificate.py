from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2380_s1330_front_loaded_profile_rectangle_monotonicity_certificate.py"
OUT = ROOT / "generated" / "p2380_s1330_front_loaded_profile_rectangle_monotonicity_certificate.json"
MD = ROOT / "generated" / "p2380_s1330_front_loaded_profile_rectangle_monotonicity_certificate.md"

PREREQ_SCRIPTS = [
    ROOT / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.py",
    ROOT / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.py",
    ROOT / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.py",
    ROOT / "p2366_s1316_selector_candidate_phase_reference_chi11_audit_probe.py",
    ROOT / "p2367_s1317_selector_phase_origin_admissibility_no_go_probe.py",
    ROOT / "p2368_s1318_self_recorded_endpoint_anchor_selector_candidate_probe.py",
    ROOT / "p2369_s1319_self_recorded_ledger_closed_form_uniqueness_theorem.py",
    ROOT / "p2370_s1320_d5_bandpass_support_closed_form_theorem.py",
    ROOT / "p2371_s1321_aut_invariant_unit_bandpass_obstruction_theorem.py",
    ROOT / "p2372_s1322_bridge_kernel_direct_band_polarity_audit.py",
    ROOT / "p2373_s1323_bridge_kernel_polarity_correction_cone_theorem.py",
    ROOT / "p2374_s1324_damping_compression_band_polarity_candidate_theorem.py",
    ROOT / "p2375_s1325_damping_compression_polarity_interval_robustness_theorem.py",
    ROOT / "p2376_s1326_damping_compression_eta_beta_rectangle_robustness_theorem.py",
    ROOT / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.py",
    ROOT / "p2378_s1328_unit_normalized_transport_coupling_insufficiency_theorem.py",
    ROOT / "p2379_s1329_front_loaded_normalized_transport_density_existence_probe.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2380FrontLoadedProfileRectangleMonotonicityCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["front_loaded_profile_rectangle_monotonicity_certificate"]
        cls.certificate = cls.probe["front_loaded_profile_rectangle_monotonicity_certificate"]
        cls.interval = cls.certificate["interval_monotonicity_certificate"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2380_s1330_v1")
        self.assertEqual(self.payload["packet_id"], "P2380")
        self.assertEqual(self.payload["stage_id"], "S1330")
        self.assertEqual(self.payload["result_kind"], "FRONT_LOADED_PROFILE_RECTANGLE_MONOTONICITY_CERTIFICATE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_interval_monotonicity_certificate(self) -> None:
        self.assertEqual(self.interval["interval_grid_per_axis"], 128)
        self.assertEqual(self.interval["interval_box_count"], 128 * 128)
        self.assertTrue(self.interval["eta_derivative_negative_on_all_boxes"])
        self.assertTrue(self.interval["beta_tors_derivative_positive_on_all_boxes"])
        self.assertTrue(self.interval["denominator_positive_on_all_boxes"])
        self.assertTrue(self.interval["rectangle_uniform_sufficiency_certified_by_monotonicity"])
        self.assertLess(self.interval["worst_eta_derivative_hi"], 0.0)
        self.assertGreater(self.interval["worst_beta_tors_derivative_lo"], 0.0)
        self.assertEqual(self.interval["failure_counts"]["eta_derivative_sign"], 0)
        self.assertEqual(self.interval["failure_counts"]["beta_tors_derivative_sign"], 0)
        self.assertEqual(self.interval["failure_counts"]["denominator_positive"], 0)

    def test_corner_maximum_and_support_replay(self) -> None:
        corner = self.interval["corner_maximum_candidate"]
        self.assertEqual(corner["eta"], 1.8)
        self.assertEqual(corner["beta_tors"], 0.1)
        self.assertAlmostEqual(corner["lambda_threshold_gt"], 0.7916644842269429)
        self.assertLess(corner["lambda_threshold_gt"], 0.8)
        rows = self.certificate["sample_support_replay"]
        self.assertEqual(len(rows), 9)
        for row in rows:
            self.assertGreater(row["lambda_test_margin"], 0.0)
            self.assertTrue(row["front_loaded_score_audit"]["d5_chamber"])
            self.assertEqual(row["front_loaded_score_audit"]["maximizer_h1_h5_pair_distribution"], {"0,4": 12})

    def test_gatekeepers_fingerprint_limits_and_next_step(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("Interval arithmetic over the P2376 rectangle certifies dT/deta<0", theorem["claim"])
        self.assertIn("strict variational source theorem deriving rho_lambda or lambda=0.8", theorem["not_licensed"])
        self.assertIn("derive rho_lambda and lambda>=0.8 from strict dynamics", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
