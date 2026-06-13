from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2677_s1627_s6_o3_typed_seed_route_no_go_audit.py"
OUT = ROOT / "generated" / "p2677_s1627_s6_o3_typed_seed_route_no_go_audit.json"
MD = ROOT / "generated" / "p2677_s1627_s6_o3_typed_seed_route_no_go_audit.md"


class P2677S6O3TypedSeedRouteNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_covers_orientation_source_escape_hatch(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "xor_xnor_reversal_content",
            "internal_selector_source_content",
            "sigma_f301_arrow_content",
            "typed_seed_route_content",
            "collapse_forbidden_content",
            "fiat_convention_blocker_content",
            "closure_guard_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_upstream_p2676_consistency(self) -> None:
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2676_formal_xor_xnor_pair_exists"])
        self.assertTrue(upstream["p2676_export_gate_zero"])
        self.assertTrue(upstream["p2676_no_sigma_to_f301_arrow"])
        self.assertTrue(upstream["p2676_s6_not_passed"])
        self.assertTrue(upstream["p2676_o3_not_passed"])

    def test_orientation_source_candidates_all_blocked(self) -> None:
        rows = self.payload["orientation_source_candidate_table"]
        self.assertEqual(len(rows), 8)
        self.assertTrue(all(row["mentioned_or_search_visible"] for row in rows))
        self.assertFalse(any(row["passes_orientation_source_gate"] for row in rows))
        blockers = " ".join(row["blocker"] for row in rows)
        for phrase in ("legacy-role", "theta-like", "Q_basis", "projector", "fiat", "symmetry", "observer"):
            self.assertIn(phrase, blockers)

    def test_finite_no_go_lattice(self) -> None:
        lattice = self.payload["finite_no_go_lattice"]
        self.assertEqual(lattice["total_states"], 128)
        self.assertEqual(lattice["passing_states"], 1)
        self.assertEqual(lattice["hamming_distance_to_pass"], 5)
        self.assertTrue(lattice["current_state"]["formal_sigma_f301_reversal_pair_identified"])
        self.assertTrue(lattice["current_state"]["route_closure_guards_preserved"])
        self.assertIn("strict_internal_orientation_source_exported", lattice["missing_current_obligations"])
        self.assertIn("nonfiat_orientation_choice_proof", lattice["missing_current_obligations"])

    def test_bounded_no_go_not_global_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["passing_orientation_source_candidates"], 0)
        self.assertFalse(decision["global_no_go_claimed"])
        self.assertFalse(decision["s6_current_route_passed_now"])
        self.assertFalse(decision["o3_current_route_passed_now"])
        self.assertFalse(decision["o4_o5_allowed_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2677/S1627", MD.read_text(encoding="utf-8"))
        self.assertIn("P2677/S1627", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2677/S1627", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
