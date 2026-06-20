import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2985_s1935_z12_leibniz_derivation_source_candidate_obstruction.py"
OUT = ROOT / "generated" / "p2985_s1935_z12_leibniz_derivation_source_candidate_obstruction.json"
MD = ROOT / "generated" / "p2985_s1935_z12_leibniz_derivation_source_candidate_obstruction.md"

class P2985Z12LeibnizDerivationSourceCandidateObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2985_Z12_LEIBNIZ_DERIVATION_SOURCE_CANDIDATE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2984"])

    def test_derivation_certificate(self):
        cert = self.payload["derivation_certificate"]
        self.assertEqual(cert["modulus"], 12)
        self.assertEqual(cert["candidate_count"], 12)
        self.assertEqual(cert["products_tested_per_candidate"], 144)
        self.assertEqual(cert["accepted_derivations"], ["D_0"])
        self.assertEqual(cert["nonzero_accepted_derivations"], [])
        self.assertEqual(cert["unit_product_defects"]["0"], 0)
        self.assertNotEqual(cert["unit_product_defects"]["1"], 0)
        self.assertEqual(cert["acceptance_matrix_rows"], 256)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_additive_endomorphism_scan"])
        self.assertTrue(obligations["full_Leibniz_product_scan"])
        self.assertTrue(obligations["derivation_algebra_computed"])
        self.assertFalse(obligations["nonzero_strict_flow_exported"])
        self.assertFalse(obligations["orientation_or_selector_source"])
        self.assertFalse(obligations["positive_scale_or_damping_source"])
        self.assertFalse(obligations["bridge_completion_map"])
        self.assertFalse(obligations["unit_bearing_action_density"])
        self.assertFalse(obligations["accepted_current_source_candidate"])
        rows = {r["candidate"]: r for r in self.payload["constructed_theoretical_objects"]["candidate_rows"]}
        self.assertTrue(rows["D_0"]["leibniz_holds_all_144_products"])
        self.assertFalse(rows["D_1"]["leibniz_holds_all_144_products"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2985/S1935", MD.read_text(encoding="utf-8"))
        self.assertIn("P2985/S1935", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2985/S1935", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2985", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
