from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2372_s1322_bridge_kernel_direct_band_polarity_audit.py"
OUT = ROOT / "generated" / "p2372_s1322_bridge_kernel_direct_band_polarity_audit.json"
MD = ROOT / "generated" / "p2372_s1322_bridge_kernel_direct_band_polarity_audit.md"

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
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2372BridgeKernelDirectBandPolarityAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["bridge_kernel_direct_band_polarity_audit"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2372_s1322_v1")
        self.assertEqual(self.payload["packet_id"], "P2372")
        self.assertEqual(self.payload["stage_id"], "S1322")
        self.assertEqual(self.payload["result_kind"], "BRIDGE_KERNEL_DIRECT_BAND_POLARITY_AUDIT")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_strict_direct_kernel_fails_d5_chamber(self) -> None:
        audit = self.probe["polarity_audit"]
        self.assertGreater(audit["strict_kernel_values"]["K1_over_K5"], 1 / 3)
        strict_direct = audit["candidate_score_audits"]["strict_direct_kernel_pair_weight"]
        self.assertFalse(strict_direct["d5_chamber_a_over_b_lt_1_over_3"])
        self.assertEqual(strict_direct["maximizer_h1_h5_pair_distribution"], {"4,0": 12})
        self.assertEqual(audit["candidate_score_audits"]["apd_completed_direct_pair_weight_equals_strict"]["maximizer_h1_h5_pair_distribution"], {"4,0": 12})

    def test_legacy_controls_and_gatekeepers(self) -> None:
        audit = self.probe["polarity_audit"]
        self.assertLess(audit["legacy_amplitude_normalized_kernel_values"]["K5"], 0)
        legacy_abs = audit["candidate_score_audits"]["legacy_amplitude_normalized_absolute_pair_weight"]
        self.assertFalse(legacy_abs["d5_chamber_a_over_b_lt_1_over_3"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_fingerprint_and_limits(self) -> None:
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("direct distance-1/distance-5 pair weights", theorem["claim"])
        self.assertIn("QW-2191 discharge", theorem["not_licensed"])
        self.assertIn("Do not use direct K(d) pair weights", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
