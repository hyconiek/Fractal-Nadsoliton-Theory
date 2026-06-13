from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2686_s1636_shared_background_nonproxy_component_residual_table.py"
OUT = ROOT / "generated" / "p2686_s1636_shared_background_nonproxy_component_residual_table.json"
MD = ROOT / "generated" / "p2686_s1636_shared_background_nonproxy_component_residual_table.md"


class P2686SharedBackgroundNonproxyComponentResidualTableTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_component_residual_frontier(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("component residual table", audit["mode"])
        for key in (
            "shared_background_nonproxy_runpack",
            "component_residual_table",
            "bianchi_anisotropic_obstruction",
            "gate_locks",
            "forbidden_closure_claims",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_current_state_uses_p2685_and_unified_runpack(self) -> None:
        state = self.payload["current_state"]
        self.assertTrue(state["p2685_next_was_component_table"])
        self.assertEqual(state["shared_background_family_id"], "bg_family_v1")
        self.assertTrue(state["p1806_tg1_locked"])
        self.assertTrue(state["p1806_ea_zero"])
        self.assertFalse(state["p1806_eh_zero"])
        self.assertFalse(state["p1806_elg_zero"])

    def test_component_table_has_zero_and_nonzero_rows(self) -> None:
        table = self.payload["component_residual_table"]
        self.assertGreaterEqual(len(table), 7)
        self.assertTrue(any(row["row_id"] == "runpack_EA" and row["zero"] for row in table))
        self.assertTrue(any(row["row_id"] == "runpack_EH" and not row["zero"] for row in table))
        self.assertTrue(any(row["row_id"].startswith("bianchi_I_") and not row["zero"] for row in table))

    def test_bianchi_source_routes_are_bounded_no_go(self) -> None:
        bianchi = self.payload["bianchi_rows"]
        self.assertTrue(bianchi["minimal_source_cancels_if_admitted"])
        self.assertFalse(bianchi["strict_source_derivation_exported"])
        self.assertTrue(bianchi["positive_energy_provider_bounded_no_go"])
        self.assertTrue(bianchi["energy_neutral_transport_obstructed"])

    def test_decision_is_no_go_boundary_without_closure_exports(self) -> None:
        decision = self.payload["decision"]
        self.assertFalse(decision["zero_table_exported"])
        self.assertTrue(decision["bounded_no_go_boundary_active"])
        self.assertFalse(decision["full_nonproxy_closure_exported_now"])
        self.assertFalse(decision["role_bearing_ltotal_exported_now"])
        self.assertFalse(decision["toe_closed_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2686/S1636", MD.read_text(encoding="utf-8"))
        self.assertIn("P2686/S1636", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2686/S1636", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2686/S1636", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
