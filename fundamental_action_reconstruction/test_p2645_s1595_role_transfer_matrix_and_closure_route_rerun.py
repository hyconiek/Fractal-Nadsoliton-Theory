from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2645_s1595_role_transfer_matrix_and_closure_route_rerun.py"
OUT = ROOT / "generated" / "p2645_s1595_role_transfer_matrix_and_closure_route_rerun.json"
MD = ROOT / "generated" / "p2645_s1595_role_transfer_matrix_and_closure_route_rerun.md"


class P2645RoleTransferMatrixRerunTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("role_transfer_matrix_content", audit["patterns"])
        self.assertIn("beta_alpha_role_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))

    def test_matrix_has_exactly_one_current_pass(self) -> None:
        rows = self.payload["role_transfer_matrix"]
        passing = [row for row in rows if row["passes_current_gate"]]
        self.assertEqual([row["role"] for row in passing], ["modified_compressed_successor_role"])
        rejected = self.payload["closure_decision"]["rejected_current_entries"]
        self.assertIn("legacy_integer_node_gauge_role", rejected)
        self.assertIn("unchanged_inverse_hierarchy_role", rejected)

    def test_alpha_and_beta_roles_remain_blocked(self) -> None:
        by_role = {row["role"]: row for row in self.payload["role_transfer_matrix"]}
        self.assertIn("target_independent_beta_identity", by_role["strict_beta_source_role"]["missing_atoms"])
        self.assertIn("alpha_geo_survives_completion", by_role["legacy_ew_angle_alpha_geo_role"]["missing_atoms"])
        self.assertIn("beta_tors_maps_to_strict_beta", by_role["legacy_alpha_em_beta_tors_role"]["missing_atoms"])
        self.assertFalse(self.payload["closure_decision"]["beta_source_route_ready"])

    def test_closure_guards_and_docs_are_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertFalse(decision["full_legacy_role_transfer_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["q_w_2191_discharged_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2645/S1595", MD.read_text(encoding="utf-8"))
        self.assertIn("P2645/S1595", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2645/S1595", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
