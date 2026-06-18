import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2835_s1785_combined_witness_source_law_theorem_obligation_audit.py"
JSON_PATH = ROOT / "generated" / "p2835_s1785_combined_witness_source_law_theorem_obligation_audit.json"
MD_PATH = ROOT / "generated" / "p2835_s1785_combined_witness_source_law_theorem_obligation_audit.md"
MANIFEST_PATH = ROOT / "generated" / "p2835_s1785_combined_witness_separator_manifest.json"


class P2835CombinedWitnessSourceLawTheoremObligationAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

    def test_status(self):
        self.assertEqual(
            self.payload["status"],
            "P2835_COMBINED_FINITE_SEPARATOR_THEOREM_OBLIGATION_NO_GO_NO_COUPLING_NO_CLOSURE",
        )

    def test_combined_separator_counts(self):
        audit = self.payload["combined_witness_source_law_theorem_obligation_audit"]
        separator = audit["combined_separator"]
        self.assertEqual(audit["decoded_full_carrier_graph_count"], 16828)
        self.assertEqual(separator["p2833_class_count"], 16757)
        self.assertEqual(separator["p2833_residual_digest_count"], 67)
        self.assertEqual(separator["p2834_patch_graph_count"], 138)
        self.assertEqual(separator["combined_class_count"], 16828)
        self.assertEqual(separator["combined_collision_class_count"], 0)
        self.assertEqual(separator["combined_collision_graph_count"], 0)
        self.assertEqual(separator["combined_class_size_histogram"]["1"], 16828)

    def test_manifest_shape(self):
        self.assertEqual(self.manifest["status"], "P2835_COMBINED_FINITE_SEPARATOR_MANIFEST")
        self.assertEqual(self.manifest["full_carrier_graph_count"], 16828)
        self.assertEqual(self.manifest["p2833_class_count"], 16757)
        self.assertEqual(self.manifest["p2834_patch_graph_count"], 138)
        self.assertEqual(self.manifest["combined_class_count"], 16828)
        self.assertEqual(self.manifest["combined_separator_row_count"], 16828)
        self.assertEqual(len(self.manifest["patched_graph_indices_0_based"]), 138)
        self.assertEqual(len(self.manifest["combined_separator_rows_sha256"]), 64)

    def test_theorem_obligation_boundaries(self):
        audit = self.payload["combined_witness_source_law_theorem_obligation_audit"]
        missing = audit["missing_blocking_obligations"]
        self.assertEqual(
            missing,
            [
                "target_independent_units_and_normalization",
                "typed_domain_codomain_into_kernel_or_lagrangian_density",
                "proved_variational_derivative_into_K_or_Ltotal",
                "coupling_source_law_with_units_and_coefficient",
            ],
        )
        statuses = {row["obligation"]: row["status"] for row in audit["theorem_obligation_matrix"]}
        self.assertEqual(statuses["finite_full_carrier_separation"], "satisfied")
        self.assertEqual(statuses["non_label_graph_formula"], "satisfied_as_finite_piecewise_formula")
        self.assertEqual(statuses["selector_bridge_role_transfer_independence"], "satisfied_by_exclusion")

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        self.assertTrue(matrix["facts"]["full_carrier_count_matches"])
        self.assertTrue(matrix["facts"]["p2834_patch_count_matches"])
        self.assertTrue(matrix["facts"]["combined_separator_is_full_carrier_injective"])
        self.assertTrue(matrix["facts"]["blocking_theorem_obligations_remain"])
        self.assertFalse(matrix["facts"]["selector_bridge_or_role_transfer_imported"])
        self.assertTrue(matrix["accepted_as_combined_finite_separator"])
        self.assertFalse(matrix["accepted_as_source_law_coupling"])
        self.assertTrue(matrix["accepted_as_theorem_obligation_no_go"])

    def test_documents_updated(self):
        self.assertIn("P2835/S1785", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2835/S1785", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2835/S1785", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2835", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
