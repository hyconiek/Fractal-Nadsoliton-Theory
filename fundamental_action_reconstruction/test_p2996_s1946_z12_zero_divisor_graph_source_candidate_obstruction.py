import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2996_s1946_z12_zero_divisor_graph_source_candidate_obstruction.py"
OUT = ROOT / "generated" / "p2996_s1946_z12_zero_divisor_graph_source_candidate_obstruction.json"
MD = ROOT / "generated" / "p2996_s1946_z12_zero_divisor_graph_source_candidate_obstruction.md"

class P2996ZeroDivisorGraphSourceCandidateObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2996_Z12_ZERO_DIVISOR_GRAPH_SOURCE_CANDIDATE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2995"])

    def test_graph_certificate(self):
        cert = self.payload["graph_certificate"]
        self.assertEqual(cert["vertex_count"], 7)
        self.assertEqual(cert["edge_count"], 8)
        self.assertEqual(cert["degree_sequence"], [1, 1, 2, 2, 3, 3, 4])
        self.assertEqual(cert["component_count"], 1)
        self.assertEqual(cert["triangle_count"], 0)
        self.assertEqual(cert["automorphism_group_order"], 8)
        self.assertEqual(cert["accepted_strict_sources"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 128)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_zero_divisor_graph_enumerated"])
        self.assertTrue(obligations["exact_graph_invariants_computed"])
        self.assertFalse(obligations["strict_nadsoliton_provenance"])
        self.assertFalse(obligations["nonpremise_vertex_edge_localizer"])
        self.assertFalse(obligations["named_source_atom_coupling"])
        self.assertFalse(obligations["unit_bearing_action_installation"])
        self.assertFalse(obligations["accepted_current_strict_source"])
        witness = obj["zero_divisor_graph_witness"]
        self.assertEqual(witness["vertices"], [2, 3, 4, 6, 8, 9, 10])
        self.assertEqual(witness["edges"], [[2, 6], [3, 4], [3, 8], [4, 6], [4, 9], [6, 8], [6, 10], [8, 9]])
        self.assertFalse(any(r["accepted_strict_source"] for r in witness["adjacency_rows"]))

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2996/S1946", MD.read_text(encoding="utf-8"))
        self.assertIn("P2996/S1946", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2996/S1946", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2996", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
