import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2997_s1947_zero_divisor_graph_strict_provenance_obstruction.py"
OUT = ROOT / "generated" / "p2997_s1947_zero_divisor_graph_strict_provenance_obstruction.json"
MD = ROOT / "generated" / "p2997_s1947_zero_divisor_graph_strict_provenance_obstruction.md"

class P2997ZeroDivisorGraphStrictProvenanceObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2997_ZERO_DIVISOR_GRAPH_STRICT_PROVENANCE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2996"])

    def test_provenance_certificate(self):
        cert = self.payload["provenance_certificate"]
        self.assertEqual(cert["vertex_count"], 7)
        self.assertEqual(cert["edge_count"], 8)
        self.assertEqual(cert["edge_rows"], 8)
        self.assertEqual(cert["nonedge_rows"], 13)
        self.assertTrue(cert["all_edges_zero_products"])
        self.assertTrue(cert["all_nonedges_nonzero_products"])
        self.assertTrue(cert["all_unit_actions_preserve_graph"])
        self.assertEqual(cert["vertex_set_preserving_translations"], [0])
        self.assertFalse(cert["accepted_strict_provenance"])
        self.assertEqual(cert["acceptance_matrix_rows"], 64)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_zero_product_graph_exact"])
        self.assertTrue(obligations["unit_action_graph_symmetry_verified"])
        self.assertTrue(obligations["translation_noninvariance_witness"])
        self.assertFalse(obligations["strict_nadsoliton_source_map_exported"])
        self.assertFalse(obligations["nonpremise_internal_graph_provenance"])
        self.assertFalse(obligations["accepted_current_strict_provenance"])
        witness = obj["provenance_witness"]
        self.assertEqual(len(witness["unit_action_rows"]), 4)
        self.assertTrue(all(r["preserves_edges"] for r in witness["unit_action_rows"]))
        self.assertFalse(witness["accepted_strict_provenance"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2997/S1947", MD.read_text(encoding="utf-8"))
        self.assertIn("P2997/S1947", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2997/S1947", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2997", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
