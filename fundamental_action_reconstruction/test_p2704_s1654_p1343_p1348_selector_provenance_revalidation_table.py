from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2704_s1654_p1343_p1348_selector_provenance_revalidation_table.py"
OUT = ROOT / "generated" / "p2704_s1654_p1343_p1348_selector_provenance_revalidation_table.json"
MD = ROOT / "generated" / "p2704_s1654_p1343_p1348_selector_provenance_revalidation_table.md"


class P2704P1343P1348SelectorProvenanceRevalidationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_packet_basis_extracts_exact_selector_object(self) -> None:
        basis = self.payload["packet_basis"]
        self.assertTrue(basis["selector_object_declared"])
        self.assertTrue(basis["operator_basis_declared"])
        self.assertEqual(basis["coefficients_raw"], "1.0,0.35,0.20,0.15")
        self.assertTrue(basis["positive_weight_condition_declared"])
        self.assertTrue(basis["p1348_depends_on_full_validation_chain"])

    def test_generated_report_statuses_match_expected(self) -> None:
        rows = self.payload["report_status_table"]
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(row["passes"] for row in rows))
        self.assertIn("S_strict_internal_v1", rows[0]["key_objects"])

    def test_p1344_csv_computation_is_finite_and_passes(self) -> None:
        comp = self.payload["p1344_computation"]
        self.assertTrue(comp["passes"])
        self.assertEqual(comp["csv_recomputed"]["total_rows"], 12480)
        self.assertEqual(comp["csv_recomputed"]["admissible_rows"], 3216)
        self.assertEqual(comp["csv_recomputed"]["sign_flips"], 0)
        self.assertGreater(comp["csv_recomputed"]["minimum_selector_margin_abs"], comp["summary"]["ambiguity_alarm_threshold"])

    def test_revalidation_is_declared_scope_only_no_false_pass(self) -> None:
        matrix = self.payload["revalidation_matrix"]
        self.assertEqual(len(matrix), 5)
        self.assertTrue(all(row["passes"] for row in matrix))
        decision = self.payload["decision"]
        self.assertTrue(decision["p1343_p1348_declared_scope_provenance_revalidated"])
        self.assertIn("declared_scope", decision["current_qw2191_status_reading"])
        self.assertIn("P2705", decision["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_docs_updated(self) -> None:
        self.assertIn("P2704/S1654", MD.read_text(encoding="utf-8"))
        self.assertIn("P2704/S1654", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2704/S1654", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2704/S1654", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
