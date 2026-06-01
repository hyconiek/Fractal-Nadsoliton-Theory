from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2369_s1319_self_recorded_ledger_closed_form_uniqueness_theorem.py"
OUT = ROOT / "generated" / "p2369_s1319_self_recorded_ledger_closed_form_uniqueness_theorem.json"
MD = ROOT / "generated" / "p2369_s1319_self_recorded_ledger_closed_form_uniqueness_theorem.md"

PREREQ_SCRIPTS = [
    ROOT / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.py",
    ROOT / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.py",
    ROOT / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.py",
    ROOT / "p2366_s1316_selector_candidate_phase_reference_chi11_audit_probe.py",
    ROOT / "p2367_s1317_selector_phase_origin_admissibility_no_go_probe.py",
    ROOT / "p2368_s1318_self_recorded_endpoint_anchor_selector_candidate_probe.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2369SelfRecordedLedgerClosedFormUniquenessTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["self_recorded_ledger_closed_form_uniqueness_theorem"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2369_s1319_v1")
        self.assertEqual(self.payload["packet_id"], "P2369")
        self.assertEqual(self.payload["stage_id"], "S1319")
        self.assertEqual(
            self.payload["result_kind"],
            "SELF_RECORDED_LEDGER_CLOSED_FORM_UNIQUENESS_THEOREM",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_closed_form_ripple_certificate(self) -> None:
        cert = self.probe["closed_form_ripple_certificate"]
        self.assertEqual(cert["base_value"], 1)
        self.assertEqual(cert["remainder"], 3)
        self.assertEqual(cert["closed_form_lower_bound"], 14)
        self.assertEqual(cert["permutation_minimizer_count"], 10)
        self.assertEqual(cert["brute_force_composition_count_cross_check"], 35)
        self.assertTrue(cert["closed_form_matches_bruteforce"])
        self.assertTrue(all(delta > 0 for delta in cert["observed_positive_smoothing_deltas"]))

    def test_arrow_tiebreak_and_endpoint_anchor(self) -> None:
        arrow = self.probe["arrow_tiebreak_certificate"]
        self.assertTrue(arrow["arrow_zero_iff_nonincreasing_on_minimizer_permutations"])
        self.assertEqual(arrow["minimum_arrow"], 0)
        self.assertEqual(arrow["winner_count"], 1)
        self.assertEqual(arrow["winners"], [[2, 2, 2, 1, 1]])
        self.assertTrue(arrow["unique_winner_is_balanced_ledger"])

        endpoint = self.probe["endpoint_anchor_closed_form_certificate"]
        self.assertEqual(endpoint["endpoint_values"], [2, 1])
        self.assertTrue(endpoint["endpoint_values_are_distinct"])
        self.assertIn("unordered value multiset", endpoint["why_internal_order_is_required"])

    def test_eta_identity_gatekeepers_and_limits(self) -> None:
        eta = self.probe["eta_identity_certificate"]
        self.assertEqual(eta["q_power_product"], "256/243")
        self.assertEqual(eta["lhs_power_integer"], eta["rhs_power_integer"])
        self.assertTrue(eta["exact_eta_identity_holds"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("convex smoothing", theorem["claim"])
        self.assertIn("QW-2191 discharge", theorem["not_licensed"])
        self.assertIn("ordered d5 support", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
