import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2836_s1786_combined_graph_functional_units_normalization_obligation_audit.py"
JSON_PATH = ROOT / "generated" / "p2836_s1786_combined_graph_functional_units_normalization_obligation_audit.json"
MD_PATH = ROOT / "generated" / "p2836_s1786_combined_graph_functional_units_normalization_obligation_audit.md"


class P2836CombinedGraphFunctionalUnitsNormalizationObligationAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status(self):
        self.assertEqual(
            self.payload["status"],
            "P2836_UNITS_NORMALIZATION_OBLIGATION_NO_GO_NO_COUPLING_NO_CLOSURE",
        )

    def test_separator_rechecked(self):
        audit = self.payload["combined_graph_functional_units_normalization_obligation_audit"]
        sep = audit["combined_separator_rechecked"]
        self.assertEqual(sep["full_carrier_graph_count"], 16828)
        self.assertEqual(sep["combined_class_count"], 16828)
        self.assertEqual(sep["combined_collision_class_count"], 0)
        self.assertEqual(sep["p2834_patch_graph_count"], 138)
        self.assertEqual(len(sep["combined_separator_rows_sha256"]), 64)

    def test_finite_normalization_candidates(self):
        audit = self.payload["combined_graph_functional_units_normalization_obligation_audit"]
        candidates = {row["candidate"]: row for row in audit["finite_combinatorial_normalization_candidates"]}
        self.assertEqual(set(candidates), {
            "per_vertex",
            "per_edge",
            "per_unordered_pair",
            "per_two_edge_toggle_pair",
            "per_full_carrier_graph",
            "per_patched_residual_graph",
        })
        self.assertEqual(candidates["per_vertex"]["normalization_factor"]["denominator"], 16)
        self.assertEqual(candidates["per_edge"]["normalization_factor"]["denominator"], 32)
        self.assertEqual(candidates["per_unordered_pair"]["normalization_factor"]["denominator"], 120)
        self.assertEqual(candidates["per_two_edge_toggle_pair"]["normalization_factor"]["denominator"], 7140)
        self.assertTrue(all(row["status"] == "finite_dimensionless_available" for row in candidates.values()))
        self.assertFalse(any(row["exports_physical_units"] for row in candidates.values()))

    def test_scale_orbit_and_missing_obligations(self):
        audit = self.payload["combined_graph_functional_units_normalization_obligation_audit"]
        scale = audit["scale_orbit_witness"]
        self.assertTrue(scale["separation_invariant_under_positive_rescaling"])
        self.assertTrue(scale["finite_witness_order_invariant_under_positive_rescaling"])
        self.assertFalse(scale["canonical_representative_exported_by_current_artifacts"])
        self.assertEqual(
            audit["missing_blocking_obligations"],
            [
                "target_independent_physical_units",
                "canonical_scale_orbit_quotient",
                "coupling_coefficient_with_units",
            ],
        )

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        facts = matrix["facts"]
        self.assertTrue(facts["p2835_combined_separator_rechecked"])
        self.assertTrue(facts["finite_dimensionless_normalizations_available"])
        self.assertFalse(facts["target_independent_physical_units_exported"])
        self.assertFalse(facts["canonical_scale_orbit_quotient_exported"])
        self.assertFalse(facts["coupling_coefficient_with_units_exported"])
        self.assertFalse(facts["selector_bridge_or_role_transfer_imported"])
        self.assertTrue(matrix["accepted_as_finite_dimensionless_normalization_audit"])
        self.assertFalse(matrix["accepted_as_target_independent_units_source"])
        self.assertTrue(matrix["accepted_as_units_normalization_no_go"])

    def test_documents_updated(self):
        self.assertIn("P2836/S1786", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2836/S1786", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2836/S1786", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2836", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
