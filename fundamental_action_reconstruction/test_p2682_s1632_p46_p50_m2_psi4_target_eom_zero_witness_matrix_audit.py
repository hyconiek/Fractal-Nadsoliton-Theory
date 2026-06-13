from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2682_s1632_p46_p50_m2_psi4_target_eom_zero_witness_matrix_audit.py"
OUT = ROOT / "generated" / "p2682_s1632_p46_p50_m2_psi4_target_eom_zero_witness_matrix_audit.json"
MD = ROOT / "generated" / "p2682_s1632_p46_p50_m2_psi4_target_eom_zero_witness_matrix_audit.md"


class P2682P46P50M2Psi4TargetEomZeroWitnessMatrixAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_action_eom_and_forbidden_imports(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "p46_target_action_frontier",
            "ax12_action_closure_boundary",
            "p50_target_eom_frontier",
            "forbidden_imports_and_closure",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_reads_current_state_action_closed_eom_not_closed(self) -> None:
        read = self.payload["action_eom_artifact_read"]
        self.assertTrue(read["ax12_external_action_zero_witness_present"])
        self.assertFalse(read["ax12_strict_core_promotion"])
        self.assertTrue(read["r39_target_eom_defect_packet_present"])
        self.assertFalse(read["r39_target_eom_zero_witness_present"])
        self.assertEqual(read["r39_defect_polynomial"], "(m2_psi4) - (mu_m2_plus3_segment_psi1_psi4)")
        self.assertEqual(read["r39_common_support"], "psi4(x)")

    def test_symbolic_matrix_blocks_import_without_transport_or_identity(self) -> None:
        matrix = self.payload["symbolic_defect_matrix"]
        self.assertEqual(matrix["total_states"], 64)
        self.assertGreater(matrix["passing_states"], 0)
        self.assertTrue(matrix["bounded_no_go_now"])
        self.assertFalse(matrix["current_state"]["coefficient_identity_m2_equals_mu_exported_for_target_eom"])
        self.assertFalse(matrix["current_state"]["action_to_eom_transport_theorem_exported"])
        self.assertFalse(matrix["current_state"]["zero_witness_exported_without_field_division"])

    def test_upstream_p2681_and_no_closure(self) -> None:
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2681_exists"])
        self.assertTrue(upstream["p2681_selected_finite_p46_lane"])
        self.assertTrue(upstream["p2681_toe_not_closed"])
        decision = self.payload["closure_decision"]
        self.assertFalse(decision["p46_target_action_reopened"])
        self.assertFalse(decision["p50_target_eom_zero_witness_exported_now"])
        self.assertFalse(decision["strict_core_promoted"])
        self.assertFalse(decision["toe_closed_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_docs_updated(self) -> None:
        self.assertIn("P2682/S1632", MD.read_text(encoding="utf-8"))
        self.assertIn("P2682/S1632", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2682/S1632", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
