from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2687_s1637_one_new_strict_anisotropic_source_class_audit.py"
OUT = ROOT / "generated" / "p2687_s1637_one_new_strict_anisotropic_source_class_audit.json"
MD = ROOT / "generated" / "p2687_s1637_one_new_strict_anisotropic_source_class_audit.md"


class P2687OneNewStrictAnisotropicSourceClassAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_source_class_frontier(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("one-new strict anisotropic source-class audit", audit["mode"])
        for key in (
            "p2686_continuation_rule",
            "derived_lapse_energy_source",
            "non_energy_neutral_tensor_transport",
            "provider_nonexport_registry",
            "forbidden_imports",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_current_state_uses_p2686_and_blockers(self) -> None:
        state = self.payload["current_state"]
        self.assertTrue(state["p2686_requested_one_new_source_class"])
        self.assertTrue(state["p1975_minimal_source_cancels_if_admitted"])
        self.assertFalse(state["p1975_strict_source_derivation_exported"])
        self.assertTrue(state["p1977_bounded_no_go_passed"])
        self.assertTrue(state["p1978_bounded_obstruction_passed"])
        self.assertTrue(state["p1976_nonexport_audit_passed"])

    def test_source_class_rows_are_near_miss_not_exports(self) -> None:
        rows = self.payload["source_class_rows"]
        ids = {row["source_class_id"] for row in rows}
        self.assertEqual(ids, {"derived_lapse_energy_source", "non_energy_neutral_tensorial_shear_transport"})
        for row in rows:
            self.assertFalse(row["typed_source_exported_now"])
            self.assertFalse(row["evades_blocker_now"])
            self.assertEqual(row["verdict"], "NEAR_MISS_NOT_EXPORTED")
        self.assertTrue(any("u0 =" in row["symbolic_requirement"] for row in rows))

    def test_obligation_matrix_blocks_reverse_closure_reopen(self) -> None:
        decision = self.payload["decision"]
        self.assertFalse(decision["new_source_class_exported"])
        self.assertTrue(decision["all_current_evasion_obligations_failed"])
        self.assertTrue(decision["freeze_lagrangian_eom_reverse_closure_lane"])
        self.assertFalse(decision["ltotal_promoted_now"])
        self.assertFalse(decision["toe_closed_now"])
        self.assertFalse(decision["selector_bridge_role_imported_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_documents_updated_without_closure_exports(self) -> None:
        self.assertIn("P2687/S1637", MD.read_text(encoding="utf-8"))
        self.assertIn("P2687/S1637", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2687/S1637", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2687/S1637", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))
        self.assertIn("P46/N49", self.payload["decision"]["next_honest_step"])


if __name__ == "__main__":
    unittest.main()
