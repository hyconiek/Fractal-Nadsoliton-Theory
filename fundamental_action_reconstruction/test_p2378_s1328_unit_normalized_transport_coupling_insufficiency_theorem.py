from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2378_s1328_unit_normalized_transport_coupling_insufficiency_theorem.py"
OUT = ROOT / "generated" / "p2378_s1328_unit_normalized_transport_coupling_insufficiency_theorem.json"
MD = ROOT / "generated" / "p2378_s1328_unit_normalized_transport_coupling_insufficiency_theorem.md"

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
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2378UnitNormalizedTransportCouplingInsufficiencyTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["unit_normalized_transport_coupling_insufficiency_theorem"]
        cls.certificate = cls.probe["unit_normalized_transport_coupling_insufficiency_certificate"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2378_s1328_v1")
        self.assertEqual(self.payload["packet_id"], "P2378")
        self.assertEqual(self.payload["stage_id"], "S1328")
        self.assertEqual(self.payload["result_kind"], "UNIT_NORMALIZED_TRANSPORT_COUPLING_INSUFFICIENCY_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_rectangle_mass_budget_proof(self) -> None:
        proof = self.certificate["rectangle_proof"]
        self.assertGreater(proof["denominator_min_value"], 0)
        self.assertLess(proof["denominator_max_value"], proof["numerator_3K1_minus_K5"])
        self.assertTrue(proof["unit_mass_insufficient_on_rectangle"])
        self.assertGreater(proof["tau_threshold_range"]["minimum_tau_gt"], 1.0)
        self.assertGreater(proof["tau_threshold_range"]["maximum_tau_gt"], proof["tau_threshold_range"]["minimum_tau_gt"])
        self.assertEqual(proof["denominator_min_corner"], {"eta": 1.8, "beta_tors": 0.1})
        self.assertEqual(proof["denominator_max_corner"], {"eta": 2.0, "beta_tors": 0.0})

    def test_unit_failure_and_threshold_success_scans(self) -> None:
        rows = self.certificate["sample_mass_budget_support_audits"]
        self.assertEqual(len(rows), 9)
        for row in rows:
            self.assertGreater(row["tau_threshold"], 1.0)
            self.assertFalse(row["unit_mass_blend_score_audit"]["d5_chamber"])
            self.assertGreater(row["unit_mass_blend_ratio"], 1.0 / 3.0)
            self.assertEqual(row["unit_mass_blend_score_audit"]["maximizer_h1_h5_pair_distribution"], {"3,3": 24})
            self.assertTrue(row["just_above_threshold_score_audit"]["d5_chamber"])
            self.assertEqual(row["just_above_threshold_score_audit"]["maximizer_h1_h5_pair_distribution"], {"0,4": 12})

    def test_gatekeepers_fingerprint_limits_and_next_step(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("unit-normalized transport primitive is insufficient", theorem["claim"])
        self.assertIn("strict variational source theorem fixing M above threshold", theorem["not_licensed"])
        self.assertIn("Do not treat P2377 transport provenance as selector closure", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
