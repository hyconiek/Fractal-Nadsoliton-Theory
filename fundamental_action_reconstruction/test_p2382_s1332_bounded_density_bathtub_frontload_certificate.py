from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2382_s1332_bounded_density_bathtub_frontload_certificate.py"
OUT = ROOT / "generated" / "p2382_s1332_bounded_density_bathtub_frontload_certificate.json"
MD = ROOT / "generated" / "p2382_s1332_bounded_density_bathtub_frontload_certificate.md"

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
    ROOT / "p2380_s1330_front_loaded_profile_rectangle_monotonicity_certificate.py",
    ROOT / "p2381_s1331_affine_frontload_source_obligation_theorem.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2382BoundedDensityBathtubFrontloadCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["bounded_density_bathtub_frontload_theorem"]
        cls.certificate = cls.probe["bounded_density_bathtub_frontload_certificate"]
        cls.obligation = cls.certificate["rectangle_worst_cap_source_obligation"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2382_s1332_v1")
        self.assertEqual(self.payload["packet_id"], "P2382")
        self.assertEqual(self.payload["stage_id"], "S1332")
        self.assertEqual(self.payload["result_kind"], "BOUNDED_DENSITY_BATHTUB_FRONTLOAD_CERTIFICATE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_bathtub_monotonicity_and_rectangle_acceptance(self) -> None:
        qprime = self.certificate["qprime_monotonicity_certificate"]
        self.assertTrue(qprime["qprime_strictly_negative_on_all_boxes"])
        self.assertLess(qprime["worst_qprime_box"]["qprime_interval"]["hi"], 0.0)
        acceptance = self.certificate["rectangle_cap_acceptance_certificate"]
        self.assertTrue(acceptance["strictly_positive_margin_on_all_boxes"])
        self.assertGreater(acceptance["worst_margin_box"]["margin_interval"]["lo"], 0.0)
        self.assertTrue(acceptance["worst_corner_replay"]["d5_chamber"])

    def test_cap_threshold_negative_control_and_positive_replay(self) -> None:
        self.assertAlmostEqual(self.obligation["cap_threshold_gt"], 1.574821357435363)
        self.assertLess(self.obligation["cap_threshold_gt"], 1.6)
        self.assertFalse(self.obligation["below_threshold_negative_control"]["d5_chamber"])
        self.assertGreater(self.obligation["below_threshold_support_score"]["a_over_b"], 1.0 / 3.0)
        self.assertTrue(self.obligation["cap_test_positive_replay"]["d5_chamber"])
        self.assertEqual(self.obligation["cap_test_support_score"]["maximizer_h1_h5_pair_distribution"], {"0,4": 12})

    def test_sample_table_dense_audit_midpoint_and_limits(self) -> None:
        rows = self.certificate["sample_local_cap_table"]
        self.assertEqual(len(rows), 9)
        for row in rows:
            self.assertLess(row["cap_threshold_gt"], 1.6)
            self.assertGreater(row["cap_test_margin"], 0.0)
            self.assertEqual(row["cap_test_support_score"]["maximizer_h1_h5_pair_distribution"], {"0,4": 12})
        dense = self.certificate["dense_cap_threshold_audit"]
        self.assertTrue(dense["cap_test_exceeds_all_grid_thresholds"])
        self.assertLess(dense["max_cap_threshold_row"]["cap_threshold_gt"], 1.6)
        midpoint = self.certificate["midpoint_replay_at_worst_corner"]
        self.assertLess(midpoint["absolute_error"], 2.0e-4)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("bathtub problem", theorem["claim"])
        self.assertIn("strict variational source theorem deriving an M-capped bang-bang density from nadsoliton dynamics", theorem["not_licensed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))


if __name__ == "__main__":
    unittest.main()
