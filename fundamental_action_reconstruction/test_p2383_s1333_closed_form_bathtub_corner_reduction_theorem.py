from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2383_s1333_closed_form_bathtub_corner_reduction_theorem.py"
OUT = ROOT / "generated" / "p2383_s1333_closed_form_bathtub_corner_reduction_theorem.json"
MD = ROOT / "generated" / "p2383_s1333_closed_form_bathtub_corner_reduction_theorem.md"

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
    ROOT / "p2382_s1332_bounded_density_bathtub_frontload_certificate.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2383ClosedFormBathtubCornerReductionTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["closed_form_bathtub_corner_reduction_theorem"]
        cls.ratio = cls.probe["closed_form_qprime_ratio_certificate"]
        cls.corner = cls.probe["cap_threshold_corner_reduction_certificate"]
        cls.burden = cls.probe["source_burden_translation"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2383_s1333_v1")
        self.assertEqual(self.payload["packet_id"], "P2383")
        self.assertEqual(self.payload["stage_id"], "S1333")
        self.assertEqual(self.payload["result_kind"], "CLOSED_FORM_BATHTUB_CORNER_REDUCTION_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_closed_form_ratio_beats_sqrt3_and_replays_qprime(self) -> None:
        self.assertGreater(self.ratio["R_min_closed_form_corner"], self.ratio["sqrt3"])
        self.assertGreater(self.ratio["R_min_square_minus_3"], 0.0)
        self.assertLess(self.ratio["qprime_factor_upper_bound_3_minus_Rmin_square"], 0.0)
        dense = self.ratio["dense_replay_worst_R"]
        self.assertAlmostEqual(dense["R"], self.ratio["R_min_closed_form_corner"], places=14)
        self.assertLess(dense["qprime"], 0.0)
        self.assertEqual(dense["eta"], 1.8)
        self.assertEqual(dense["beta_tors"], 0.0)
        self.assertEqual(dense["s"], 1.0)

    def test_cap_corner_reduction_derivative_signs_and_threshold(self) -> None:
        derivative_audits = self.corner["derivative_audits"]
        self.assertGreater(derivative_audits["margin_d_eta_on_cap_band"]["min"]["value"], 0.0)
        self.assertLess(derivative_audits["margin_d_beta_on_cap_band"]["max"]["value"], 0.0)
        self.assertGreater(derivative_audits["margin_d_cap_on_cap_band"]["min"]["value"], 0.0)
        self.assertEqual(self.corner["worst_corner"], {"eta": 1.8, "beta_tors": 0.1})
        self.assertAlmostEqual(self.corner["worst_corner_threshold"], 1.574821357435363)
        self.assertTrue(self.corner["positive_at_M_1_6"])
        self.assertTrue(self.corner["negative_at_M_1_5"])

    def test_source_burden_translation_and_gatekeepers(self) -> None:
        self.assertAlmostEqual(self.burden["minimal_rectangle_cap_gt"], 1.574821357435363)
        self.assertAlmostEqual(self.burden["M_1_6_early_interval_length"], 0.625)
        self.assertAlmostEqual(self.burden["M_1_6_early_mass_on_first_half"], 0.8)
        self.assertAlmostEqual(self.burden["M_1_6_barycenter"], 0.3125)
        self.assertAlmostEqual(self.burden["M_1_6_barycenter_shift_from_uniform"], 0.1875)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("closed ratio reduction", theorem["claim"])
        self.assertIn("strict source theorem deriving the cap M or bang-bang density", theorem["not_licensed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))


if __name__ == "__main__":
    unittest.main()
