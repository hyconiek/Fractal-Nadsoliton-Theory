import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3001_s1951_z12_square_map_functional_graph_source_candidate_obstruction.py"
OUT = ROOT / "generated" / "p3001_s1951_z12_square_map_functional_graph_source_candidate_obstruction.json"
MD = ROOT / "generated" / "p3001_s1951_z12_square_map_functional_graph_source_candidate_obstruction.md"

class P3001SquareMapFunctionalGraphSourceCandidateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3001_Z12_SQUARE_MAP_FUNCTIONAL_GRAPH_SOURCE_CANDIDATE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3000"])

    def test_functional_graph_certificate(self):
        cert = self.payload["functional_graph_certificate"]
        self.assertEqual(cert["node_count"], 12)
        self.assertEqual(cert["edge_count"], 12)
        self.assertEqual(cert["image"], [0, 1, 4, 9])
        self.assertEqual(cert["image_size"], 4)
        self.assertEqual(cert["fixed_points"], [0, 1, 4, 9])
        self.assertEqual(cert["fixed_point_count"], 4)
        self.assertEqual(cert["basin_sizes"], {"0": 2, "1": 4, "4": 4, "9": 2})
        self.assertEqual(cert["max_depth_to_attractor"], 1)
        self.assertTrue(cert["all_periods_one"])
        self.assertEqual(cert["accepted_strict_sources"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 256)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_square_map_functional_graph"])
        self.assertTrue(obligations["idempotent_attractor_decomposition"])
        self.assertTrue(obligations["fiber_basin_certificate"])
        self.assertFalse(obligations["strict_nadsoliton_provenance"])
        self.assertFalse(obligations["nonpremise_basin_attractor_localizer"])
        self.assertFalse(obligations["named_source_atom_coupling"])
        self.assertFalse(obligations["unit_bearing_action_installation"])
        self.assertFalse(obligations["accepted_current_strict_source"])
        rows = obj["square_map_witness"]["row_data"]
        self.assertEqual(len(rows), 12)
        self.assertFalse(any(r["accepted_strict_source"] for r in rows))

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P3001/S1951", MD.read_text(encoding="utf-8"))
        self.assertIn("P3001/S1951", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3001/S1951", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3001", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
