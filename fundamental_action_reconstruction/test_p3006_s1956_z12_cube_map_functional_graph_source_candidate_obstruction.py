import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3006_s1956_z12_cube_map_functional_graph_source_candidate_obstruction.py"
OUT = ROOT / "generated" / "p3006_s1956_z12_cube_map_functional_graph_source_candidate_obstruction.json"
MD = ROOT / "generated" / "p3006_s1956_z12_cube_map_functional_graph_source_candidate_obstruction.md"

class P3006Z12CubeMapFunctionalGraphSourceCandidateObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3006_Z12_CUBE_MAP_FUNCTIONAL_GRAPH_SOURCE_CANDIDATE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3005"])

    def test_cube_map_certificate(self):
        cert = self.payload["cube_map_certificate"]
        self.assertEqual(cert["node_count"], 12)
        self.assertEqual(cert["edge_count"], 12)
        self.assertEqual(cert["image"], [0, 1, 3, 4, 5, 7, 8, 9, 11])
        self.assertEqual(cert["image_size"], 9)
        self.assertEqual(cert["fixed_points"], [0, 1, 3, 4, 5, 7, 8, 9, 11])
        self.assertEqual(cert["fixed_point_count"], 9)
        self.assertEqual(cert["moved_points"], [2, 6, 10])
        self.assertEqual(cert["moved_point_count"], 3)
        self.assertEqual(cert["basin_sizes"], {"0": 2, "1": 1, "3": 1, "4": 2, "5": 1, "7": 1, "8": 2, "9": 1, "11": 1})
        self.assertEqual(cert["nontrivial_fibers"], {"0": [0, 6], "4": [4, 10], "8": [2, 8]})
        self.assertEqual(cert["max_depth_to_attractor"], 1)
        self.assertTrue(cert["all_periods_one"])
        self.assertEqual(cert["accepted_strict_sources"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 256)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_cube_map_functional_graph"])
        self.assertTrue(obligations["fixed_point_decomposition"])
        self.assertTrue(obligations["fiber_basin_certificate"])
        self.assertFalse(obligations["strict_nadsoliton_provenance"])
        self.assertFalse(obligations["nonpremise_basin_fixed_sector_localizer"])
        self.assertFalse(obligations["named_source_atom_coupling"])
        self.assertFalse(obligations["unit_bearing_action_installation"])
        self.assertFalse(obligations["accepted_current_strict_source"])
        witness = obj["cube_map_witness"]
        self.assertEqual(len(witness["row_data"]), 12)
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3006/S1956", MD.read_text(encoding="utf-8"))
        self.assertIn("P3006/S1956", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3006/S1956", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3006", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
