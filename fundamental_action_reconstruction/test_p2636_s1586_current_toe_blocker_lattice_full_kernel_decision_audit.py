from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2636_s1586_current_toe_blocker_lattice_full_kernel_decision_audit.py"
OUT = ROOT / "generated" / "p2636_s1586_current_toe_blocker_lattice_full_kernel_decision_audit.json"
MD = ROOT / "generated" / "p2636_s1586_current_toe_blocker_lattice_full_kernel_decision_audit.md"


class P2636CurrentToeBlockerLatticeFullKernelDecisionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present_and_nonempty(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        self.assertIn("full_kernel_finality_content", audit["patterns"])
        self.assertIn("phase_node_inverse_hierarchy_content", audit["patterns"])

    def test_current_blocker_lattice_lists_mandatory_open_blockers(self) -> None:
        blockers = self.payload["current_blocker_lattice"]
        self.assertEqual(len(blockers), 7)
        names = {row["blocker"] for row in blockers}
        self.assertIn("legacy_to_strict_bridge_completion", names)
        self.assertIn("strict_core_selector_source_qw2191", names)
        self.assertIn("blind_empirical_frozen_kernel_confirmation", names)
        self.assertTrue(all(row["is_closed"] is False for row in blockers))

    def test_full_kernel_decision_rejects_finality_despite_toe_symptoms(self) -> None:
        decision = self.payload["full_kernel_decision"]
        self.assertTrue(decision["gates"]["toe_symptoms_present"])
        self.assertTrue(decision["gates"]["strict_kernel_stability_positive"])
        self.assertFalse(decision["gates"]["all_mandatory_blockers_closed"])
        self.assertFalse(decision["is_full_kernel_now"])
        self.assertEqual(decision["classification"], "ROBUST_TOE_LIKE_WORKING_KERNEL_NOT_FULL_KERNEL")

    def test_minimal_blocker_cut_prioritizes_core_obstructions(self) -> None:
        cut = self.payload["minimal_blocker_cut"]
        self.assertEqual(len(cut), 4)
        blockers = [row["blocker"] for row in cut]
        self.assertIn("legacy_to_strict_bridge_completion", blockers)
        self.assertIn("strict_core_selector_source_qw2191", blockers)
        self.assertGreaterEqual(cut[0]["severity_0_to_1"], cut[-1]["severity_0_to_1"])

    def test_negative_exports_and_recommendation(self) -> None:
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        rec = self.payload["recommended_next_honest_step"]
        self.assertIn("phase-frequency/node quotient-exhaustion certificate", rec)
        self.assertIn("no post-hoc retuning", rec)

    def test_docs_updated(self) -> None:
        self.assertIn("P2636/S1586", MD.read_text(encoding="utf-8"))
        self.assertIn("P2636/S1586", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2636/S1586", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
