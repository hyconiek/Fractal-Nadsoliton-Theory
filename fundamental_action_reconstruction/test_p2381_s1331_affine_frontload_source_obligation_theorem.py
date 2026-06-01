from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2381_s1331_affine_frontload_source_obligation_theorem.py"
OUT = ROOT / "generated" / "p2381_s1331_affine_frontload_source_obligation_theorem.json"
MD = ROOT / "generated" / "p2381_s1331_affine_frontload_source_obligation_theorem.md"

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
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2381AffineFrontloadSourceObligationTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["affine_frontload_source_obligation_theorem"]
        cls.certificate = cls.probe["affine_frontload_source_obligation_certificate"]
        cls.obligation = cls.certificate["rectangle_worst_case_source_obligation"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2381_s1331_v1")
        self.assertEqual(self.payload["packet_id"], "P2381")
        self.assertEqual(self.payload["stage_id"], "S1331")
        self.assertEqual(self.payload["result_kind"], "AFFINE_FRONTLOAD_SOURCE_OBLIGATION_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_rectangle_worst_case_obligations(self) -> None:
        self.assertAlmostEqual(self.obligation["lambda_must_exceed"], 0.7916644842269429)
        obligations = self.obligation["necessary_profile_obligations_for_any_rectangle_uniform_affine_source"]
        self.assertGreater(obligations["early_half_mass_int_0_to_1_2"], 0.6979)
        self.assertLess(obligations["transport_time_barycenter_int_s_rho"], 0.3681)
        self.assertGreater(obligations["endpoint_contrast_rho0_over_rho1"], 8.59)
        self.assertGreater(obligations["l1_distance_from_uniform"], 0.395)
        self.assertGreater(obligations["l2_squared_distance_from_uniform"], 0.208)

    def test_negative_control_positive_replay_and_sample_table(self) -> None:
        self.assertFalse(self.obligation["below_threshold_negative_control_score"]["d5_chamber"])
        self.assertGreater(self.obligation["below_threshold_negative_control_score"]["a_over_b"], 1.0 / 3.0)
        self.assertTrue(self.obligation["lambda_test_corner_score"]["d5_chamber"])
        self.assertEqual(self.obligation["lambda_test_corner_score"]["maximizer_h1_h5_pair_distribution"], {"0,4": 12})
        rows = self.certificate["sample_local_obligation_table"]
        self.assertEqual(len(rows), 9)
        for row in rows:
            self.assertLess(row["lambda_threshold_gt"], 0.8)
            self.assertGreater(row["lambda_test_margin"], 0.0)
            self.assertIn("transport_time_barycenter_int_s_rho", row["necessary_obligations_at_local_threshold"])

    def test_gatekeepers_fingerprint_limits_and_next_step(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("necessary source burden", theorem["claim"])
        self.assertIn("strict variational source theorem deriving the required affine asymmetry", theorem["not_licensed"])
        self.assertIn("profile remains an explicit non-strict premise", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
