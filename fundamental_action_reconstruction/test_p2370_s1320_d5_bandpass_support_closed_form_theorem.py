from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2370_s1320_d5_bandpass_support_closed_form_theorem.py"
OUT = ROOT / "generated" / "p2370_s1320_d5_bandpass_support_closed_form_theorem.json"
MD = ROOT / "generated" / "p2370_s1320_d5_bandpass_support_closed_form_theorem.md"

PREREQ_SCRIPTS = [
    ROOT / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.py",
    ROOT / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.py",
    ROOT / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.py",
    ROOT / "p2366_s1316_selector_candidate_phase_reference_chi11_audit_probe.py",
    ROOT / "p2367_s1317_selector_phase_origin_admissibility_no_go_probe.py",
    ROOT / "p2368_s1318_self_recorded_endpoint_anchor_selector_candidate_probe.py",
    ROOT / "p2369_s1319_self_recorded_ledger_closed_form_uniqueness_theorem.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2370D5BandpassSupportClosedFormTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["d5_bandpass_support_closed_form_theorem"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2370_s1320_v1")
        self.assertEqual(self.payload["packet_id"], "P2370")
        self.assertEqual(self.payload["stage_id"], "S1320")
        self.assertEqual(self.payload["result_kind"], "D5_BANDPASS_SUPPORT_CLOSED_FORM_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_bandpass_closed_form_certificate(self) -> None:
        cert = self.probe["bandpass_closed_form_certificate"]
        self.assertEqual(cert["gcd_5_12"], 1)
        self.assertEqual(cert["support_count_cross_check"], 792)
        self.assertEqual(cert["closed_form_max_h5"], 4)
        self.assertEqual(cert["brute_force_max_h5"], 4)
        self.assertEqual(cert["maximizer_count"], 12)
        self.assertTrue(cert["maximizers_equal_step5_path_orbit"])
        self.assertTrue(cert["all_maximizers_connected_5_paths"])
        self.assertEqual(cert["h5_count_distribution"]["4"], 12)

    def test_unit_orbit_chain_and_scratch_replay(self) -> None:
        unit = self.probe["unit_step_orbit_certificate"]
        self.assertEqual(unit["aut_z12_units"], [1, 5, 7, 11])
        self.assertEqual(unit["oriented_images_of_step5"], [1, 5, 7, 11])
        self.assertFalse(unit["step5_unoriented_pair_fixed_by_aut_z12"])
        self.assertTrue(unit["step5_unoriented_pair_collapses_with_step1_under_full_aut_z12"])

        chain = self.probe["conditional_chain_certificate"]
        self.assertEqual(chain["p2369_loaded_packet"], "P2369")
        self.assertTrue(chain["p2369_arrow_tiebreak_unique"])
        self.assertTrue(chain["p2369_endpoint_anchor_distinct"])
        self.assertTrue(chain["bandpass_maximizers_equal_step5_orbit"])
        self.assertIn("band-pass action", chain["remaining_missing_premise"])

        replay = self.probe["scratch_bandpass_replay"]
        self.assertEqual(replay["max_h5"], 4)
        self.assertEqual(replay["maximizer_count"], 12)
        self.assertTrue(replay["maximizers_equal_fifth_step_orbit"])

    def test_gatekeepers_fingerprint_and_limits(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("distance-5 graph", theorem["claim"])
        self.assertIn("QW-2191 discharge", theorem["not_licensed"])
        self.assertIn("distance-5 band-pass action", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
