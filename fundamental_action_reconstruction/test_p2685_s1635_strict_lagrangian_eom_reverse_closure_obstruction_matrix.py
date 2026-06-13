from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2685_s1635_strict_lagrangian_eom_reverse_closure_obstruction_matrix.py"
OUT = ROOT / "generated" / "p2685_s1635_strict_lagrangian_eom_reverse_closure_obstruction_matrix.json"
MD = ROOT / "generated" / "p2685_s1635_strict_lagrangian_eom_reverse_closure_obstruction_matrix.md"


class P2685StrictLagrangianEomReverseClosureObstructionMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_reverse_closure_frontier(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("reverse-closure", audit["mode"])
        for key in (
            "selector_independent_lagrangian_eom",
            "full_tensor_nonproxy_obligations",
            "reverse_closure_or_helmholtz",
            "open_obstruction_witnesses",
            "forbidden_imports",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_artifact_state_confirms_real_forward_progress_but_open_full_closure(self) -> None:
        state = self.payload["current_artifact_state"]
        self.assertTrue(state["p2684_reverse_closure_pivot_present"])
        self.assertTrue(state["p2316_best_working_ltotal_identified"])
        self.assertTrue(state["p2316_full_task3_theorem_not_claimed"])
        self.assertTrue(state["p2329_all_terms_selector_independent"])
        self.assertTrue(state["p2362_no_selector_prerequisite_for_eom_export"])
        self.assertTrue(state["p1795_metric_full_tensor_closure_open"])
        self.assertTrue(state["p1809_tg1_locked_by_unified_nonproxy_residual"])

    def test_anisotropic_residual_is_symbolic_obstruction(self) -> None:
        residual = self.payload["anisotropic_residual_rank"]
        self.assertTrue(residual["obstruction_detected"])
        self.assertGreater(residual["symbolic_jacobian_rank"], 0)
        self.assertTrue(residual["isotropic_limit_zero"])
        self.assertGreaterEqual(residual["numeric_nonzero_samples"], 3)

    def test_reverse_closure_matrix_fails_without_false_pass(self) -> None:
        rows = self.payload["reverse_closure_rows"]
        self.assertEqual(len(rows), 4)
        self.assertTrue(any(row["available_forward_export"] for row in rows))
        self.assertFalse(all(row["satisfied_now"] for row in rows))
        lattice = self.payload["policy_lattice"]
        self.assertFalse(lattice["current_reverse_closure_success"])
        self.assertIn("R3_anisotropic_background_transport", self.payload["decision"]["failed_rows"])

    def test_docs_updated_and_no_closure_exports(self) -> None:
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        decision = self.payload["decision"]
        self.assertFalse(decision["full_tensor_resolved_eom_closed_now"])
        self.assertFalse(decision["role_bearing_ltotal_exported_now"])
        self.assertFalse(decision["toe_closed_now"])
        self.assertIn("P2685/S1635", MD.read_text(encoding="utf-8"))
        self.assertIn("P2685/S1635", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2685/S1635", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2685/S1635", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
