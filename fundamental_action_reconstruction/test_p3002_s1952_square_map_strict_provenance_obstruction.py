import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3002_s1952_square_map_strict_provenance_obstruction.py"
OUT = ROOT / "generated" / "p3002_s1952_square_map_strict_provenance_obstruction.json"
MD = ROOT / "generated" / "p3002_s1952_square_map_strict_provenance_obstruction.md"

class P3002SquareMapStrictProvenanceObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3002_SQUARE_MAP_STRICT_PROVENANCE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3001"])

    def test_provenance_certificate(self):
        cert = self.payload["provenance_certificate"]
        self.assertEqual(cert["node_count"], 12)
        self.assertEqual(cert["edge_count"], 12)
        self.assertEqual(cert["multiplicative_pair_count"], 144)
        self.assertTrue(cert["all_multiplicative_pairs_compatible"])
        self.assertGreater(cert["additive_defect_count"], 0)
        self.assertEqual(cert["unit_action_row_count"], 48)
        self.assertTrue(cert["all_unit_actions_square_invariant"])
        self.assertEqual(cert["graph_preserving_translations"], [0])
        self.assertFalse(cert["accepted_strict_provenance"])
        self.assertEqual(cert["acceptance_matrix_rows"], 64)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_square_graph_recomputed"])
        self.assertTrue(obligations["multiplicative_endomap_verified"])
        self.assertTrue(obligations["unit_square_invariance_verified"])
        self.assertTrue(obligations["translation_nonprovenance_witness"])
        self.assertFalse(obligations["strict_nadsoliton_source_map_exported"])
        self.assertFalse(obligations["nonpremise_internal_square_law_exported"])
        self.assertFalse(obligations["accepted_current_strict_provenance"])
        witness = obj["provenance_witness"]
        self.assertEqual(len(witness["multiplicative_rows"]), 144)
        self.assertEqual(len(witness["unit_rows"]), 48)
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3002/S1952", MD.read_text(encoding="utf-8"))
        self.assertIn("P3002/S1952", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3002/S1952", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3002", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
