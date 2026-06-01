from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2371_s1321_aut_invariant_unit_bandpass_obstruction_theorem.py"
OUT = ROOT / "generated" / "p2371_s1321_aut_invariant_unit_bandpass_obstruction_theorem.json"
MD = ROOT / "generated" / "p2371_s1321_aut_invariant_unit_bandpass_obstruction_theorem.md"

PREREQ_SCRIPTS = [
    ROOT / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.py",
    ROOT / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.py",
    ROOT / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.py",
    ROOT / "p2366_s1316_selector_candidate_phase_reference_chi11_audit_probe.py",
    ROOT / "p2367_s1317_selector_phase_origin_admissibility_no_go_probe.py",
    ROOT / "p2368_s1318_self_recorded_endpoint_anchor_selector_candidate_probe.py",
    ROOT / "p2369_s1319_self_recorded_ledger_closed_form_uniqueness_theorem.py",
    ROOT / "p2370_s1320_d5_bandpass_support_closed_form_theorem.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2371AutInvariantUnitBandpassObstructionTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["aut_invariant_unit_bandpass_obstruction_theorem"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2371_s1321_v1")
        self.assertEqual(self.payload["packet_id"], "P2371")
        self.assertEqual(self.payload["stage_id"], "S1321")
        self.assertEqual(self.payload["result_kind"], "AUT_INVARIANT_UNIT_BANDPASS_OBSTRUCTION_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_chamber_certificate(self) -> None:
        chamber = self.probe["chamber_certificate"]
        self.assertEqual(chamber["support_count"], 792)
        self.assertEqual(chamber["d5_path_orbit_count"], 12)
        self.assertEqual(chamber["d5_path_h1_h5_pairs"], [[0, 4]])
        self.assertEqual(chamber["d1_path_h1_h5_pairs"], [[4, 0]])
        self.assertEqual(chamber["mixed_h1_eq_h5_eq_3_count"], 24)
        self.assertFalse(chamber["aut_invariant_a_eq_b_selects_d5_path_orbit"])

    def test_score_chambers(self) -> None:
        chamber = self.probe["chamber_certificate"]
        self.assertEqual(chamber["full_aut_invariant_a_eq_b_maximizers"]["maximizer_count"], 24)
        self.assertEqual(chamber["full_aut_invariant_a_eq_b_maximizers"]["maximizer_h1_h5_pair_distribution"], {"3,3": 24})
        self.assertEqual(chamber["pure_h5_bandpass_maximizers"]["maximizer_h1_h5_pair_distribution"], {"0,4": 12})
        self.assertEqual(chamber["threshold_tie_a_over_b_eq_1_over_3"]["maximizer_h1_h5_pair_distribution"], {"0,4": 12, "3,3": 24})
        self.assertEqual(chamber["d5_selecting_example_a_over_b_eq_1_over_4"]["maximizer_h1_h5_pair_distribution"], {"0,4": 12})
        self.assertIn("a/b < 1/3", chamber["closed_form_chamber_statement"])

    def test_gatekeepers_fingerprint_and_limits(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("full-Aut unit-band", theorem["claim"])
        self.assertIn("QW-2191 discharge", theorem["not_licensed"])
        self.assertIn("band-polarity inequality", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
