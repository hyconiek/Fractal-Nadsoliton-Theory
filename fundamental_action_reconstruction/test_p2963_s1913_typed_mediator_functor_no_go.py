import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2963_s1913_typed_mediator_functor_no_go.py"
OUT = ROOT / "generated" / "p2963_s1913_typed_mediator_functor_no_go.json"
MD = ROOT / "generated" / "p2963_s1913_typed_mediator_functor_no_go.md"

class P2963TypedMediatorFunctorNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2963_TYPED_MEDIATOR_FUNCTOR_NO_GO_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2962"])

    def test_mediator_classification(self):
        cert = self.payload["mediator_certificate"]
        self.assertEqual(cert["mediator_count"], 8)
        self.assertEqual(cert["equalizing_mediators"], ["max_nonzero"])
        self.assertEqual(cert["accepted_typed_mediators"], [])
        self.assertFalse(cert["strict_typed_mediator_exported"])
        rows = {r["mediator"]: r for r in self.payload["constructed_theoretical_objects"]["mediator_rows"]}
        self.assertTrue(rows["max_nonzero"]["makes_K_and_C_equal"])
        self.assertFalse(rows["max_nonzero"]["does_not_erase_P2962_mismatch"])
        self.assertFalse(rows["identity_typed_signature"]["makes_K_and_C_equal"])
        self.assertTrue(rows["identity_typed_signature"]["does_not_erase_P2962_mismatch"])

    def test_obligations_and_matrix(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_mediator_menu_constructed"])
        self.assertTrue(obligations["some_mediator_makes_K_and_C_equal"])
        self.assertFalse(obligations["equalizing_mediator_preserves_support_and_provenance"])
        self.assertFalse(obligations["strict_nadsoliton_mediator_functor_exported"])
        matrix = self.payload["constructed_theoretical_objects"]["finite_acceptance_matrix"]
        self.assertEqual(len(matrix), 32)
        self.assertEqual(sum(1 for row in matrix if row["accepts_strict_typed_mediator"]), 1)

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2963/S1913", MD.read_text(encoding="utf-8"))
        self.assertIn("P2963/S1913", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2963/S1913", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2963", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
