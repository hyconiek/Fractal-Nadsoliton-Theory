from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2633_s1583_diagram_grounded_strict_kernel_characteristic_preservation_audit.py"
OUT = ROOT / "generated" / "p2633_s1583_diagram_grounded_strict_kernel_characteristic_preservation_audit.json"
MD = ROOT / "generated" / "p2633_s1583_diagram_grounded_strict_kernel_characteristic_preservation_audit.md"


class P2633DiagramGroundedStrictKernelCharacteristicPreservationAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_and_diagram_source_are_present(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        source = self.payload["diagram_source"]
        self.assertIn("DIAGRAMS_KERNEL_TRANSFORMATION.md", source["path"])
        self.assertEqual(len(source["sha256"]), 64)

    def test_diagram_claims_are_extracted_by_content(self) -> None:
        claims = self.payload["diagram_claims"]["claims"]
        self.assertTrue(claims["has_four_mechanism_product"])
        self.assertTrue(claims["claims_inverse_hierarchy"])
        self.assertTrue(claims["claims_hyperbolic_not_exponential"])
        self.assertTrue(claims["claims_integer_node_pattern"])

    def test_inverse_hierarchy_is_not_preserved_on_audited_strict_grid(self) -> None:
        ratios = self.payload["finite_witness"]["inverse_hierarchy_ratio_abs_k7_over_abs_k1"]
        self.assertGreater(ratios["legacy_amplitude_normalized"], 1.0)
        self.assertLess(ratios["strict"], 1.0)
        self.assertFalse(ratios["strict_preserves_legacy_ratio_above_one"])

    def test_declared_integer_nodes_are_not_formula_exact(self) -> None:
        audit = self.payload["finite_witness"]["declared_integer_node_audit"]
        self.assertFalse(audit["all_declared_integer_nodes_are_formula_zeros"])
        self.assertTrue(all(not row["is_formula_zero_with_tolerance_1e_12"] for row in audit["rows"]))
        self.assertAlmostEqual(audit["true_continuous_zero_lattice_first_values"][0], 4.0 / 3.0)

    def test_preservation_matrix_and_acceptance_remain_nonpromoting(self) -> None:
        matrix = self.payload["preservation_matrix"]
        self.assertEqual(len(matrix), 7)
        self.assertTrue(any(row["preservation_verdict"] == "diagram_claim_not_formula_exact" for row in matrix))
        self.assertTrue(any(row["preservation_verdict"] == "not_preserved_numerically_without_reinterpretation" for row in matrix))
        acceptance = self.payload["source_acceptance"]
        self.assertFalse(acceptance["accepts_all_diagram_characteristics_preserved_by_strict"])
        self.assertIn("legacy_inverse_hierarchy_numerically_retained", acceptance["failed_gates"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_recommendation_and_docs_updated(self) -> None:
        self.assertIn("phase-frequency certificate", self.payload["recommended_next_honest_step"])
        self.assertIn("P2633/S1583", MD.read_text(encoding="utf-8"))
        self.assertIn("P2633/S1583", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2633/S1583", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
