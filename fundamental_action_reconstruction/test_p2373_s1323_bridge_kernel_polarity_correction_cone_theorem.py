from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2373_s1323_bridge_kernel_polarity_correction_cone_theorem.py"
OUT = ROOT / "generated" / "p2373_s1323_bridge_kernel_polarity_correction_cone_theorem.json"
MD = ROOT / "generated" / "p2373_s1323_bridge_kernel_polarity_correction_cone_theorem.md"

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
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2373BridgeKernelPolarityCorrectionConeTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["bridge_kernel_polarity_correction_cone_theorem"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2373_s1323_v1")
        self.assertEqual(self.payload["packet_id"], "P2373")
        self.assertEqual(self.payload["stage_id"], "S1323")
        self.assertEqual(self.payload["result_kind"], "BRIDGE_KERNEL_POLARITY_CORRECTION_CONE_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_correction_thresholds(self) -> None:
        cone = self.probe["correction_cone_certificate"]
        self.assertGreater(cone["baseline_weights"]["a0_over_b0"], 1 / 3)
        self.assertGreater(cone["minimal_corrections"]["pure_h5_boost_lambda_gt"], 0)
        self.assertGreater(cone["minimal_corrections"]["pure_h1_suppression_mu_gt"], 0)
        self.assertGreater(cone["minimal_corrections"]["antisymmetric_gamma_for_plus_h5_minus_h1_gt"], 0)
        self.assertGreater(cone["minimal_corrections"]["pure_h5_boost_lambda_gt_in_units_of_K5"], 50)

    def test_score_audits_enter_d5_chamber(self) -> None:
        audits = self.probe["correction_cone_certificate"]["score_audits"]
        self.assertEqual(audits["baseline_direct_kernel"]["maximizer_h1_h5_pair_distribution"], {"4,0": 12})
        self.assertEqual(audits["pure_h5_boost_at_open_threshold_plus_epsilon"]["maximizer_h1_h5_pair_distribution"], {"0,4": 12})
        self.assertEqual(audits["pure_h1_suppression_at_open_threshold_plus_epsilon"]["maximizer_h1_h5_pair_distribution"], {"0,4": 12})
        self.assertEqual(audits["antisymmetric_h5_minus_h1_polarity_at_open_threshold_plus_epsilon"]["maximizer_h1_h5_pair_distribution"], {"0,4": 12})

    def test_gatekeepers_fingerprint_and_limits(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("necessary correction-cone bounds", theorem["claim"])
        self.assertIn("QW-2191 discharge", theorem["not_licensed"])
        self.assertIn("actual bridge-completed dynamical source", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
