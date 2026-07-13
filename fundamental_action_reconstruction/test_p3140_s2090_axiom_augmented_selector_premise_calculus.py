import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3140_s2090_axiom_augmented_selector_premise_calculus.py"
OUT = ROOT / "generated" / "p3140_s2090_axiom_augmented_selector_premise_calculus.json"
MD = ROOT / "generated" / "p3140_s2090_axiom_augmented_selector_premise_calculus.md"


class P3140AxiomAugmentedSelectorPremiseCalculusTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_inputs_and_counts(self):
        self.assertEqual(self.payload["status"], "P3140_AXIOM_AUGMENTED_SELECTOR_PREMISE_CALCULUS_NON_STRICT")
        self.assertTrue(all(self.payload["input_hashes"].values()))
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["pair_space_size"], 24)
        self.assertEqual(counts["no_axiom_orbit_size_under_translation_aut_lambda_flip"], 24)
        self.assertEqual(counts["axiom_packages_tested"], 8)
        self.assertEqual(counts["minimal_non_strict_closing_packages"], 1)

    def test_axiom_package_lattice(self):
        rows = {tuple(row["axioms"]): row for row in self.payload["axiom_package_rows"]}
        self.assertEqual(rows[()]["surviving_pairs"], 24)
        self.assertEqual(rows[("A_origin",)]["surviving_pairs"], 2)
        self.assertEqual(rows[("A_lambda",)]["surviving_pairs"], 12)
        full = rows[("A_origin", "A_lambda", "A_coupling")]
        self.assertEqual(full["surviving_pairs"], 1)
        self.assertTrue(full["conditional_non_strict_J_DHL_exported"])
        self.assertFalse(full["strict_source_exported"])
        self.assertEqual(full["classification"], "non_strict_axiom_augmented_selector")

    def test_obligations_and_exports(self):
        obligations = {row["obligation"]: row for row in self.payload["obligation_rows"]}
        self.assertEqual(obligations["absolute_origin"]["axiom_needed"], "A_origin")
        self.assertEqual(obligations["unpaired_lambda"]["axiom_needed"], "A_lambda")
        self.assertFalse(obligations["strict_source_provenance"]["non_strict_if_assumed"])
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["conditional_non_strict_J_DHL_exported_under_full_axioms"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_updated(self):
        self.assertIn("P3140/S2090", MD.read_text(encoding="utf-8"))
        self.assertIn("P3140/S2090", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3140/S2090", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3140", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
