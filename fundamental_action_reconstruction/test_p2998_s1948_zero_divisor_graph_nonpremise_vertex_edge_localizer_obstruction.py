import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2998_s1948_zero_divisor_graph_nonpremise_vertex_edge_localizer_obstruction.py"
OUT = ROOT / "generated" / "p2998_s1948_zero_divisor_graph_nonpremise_vertex_edge_localizer_obstruction.json"
MD = ROOT / "generated" / "p2998_s1948_zero_divisor_graph_nonpremise_vertex_edge_localizer_obstruction.md"

class P2998ZeroDivisorGraphLocalizerObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2998_ZERO_DIVISOR_GRAPH_NONPREMISE_VERTEX_EDGE_LOCALIZER_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2996"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2997"])

    def test_localizer_certificate(self):
        cert = self.payload["localizer_certificate"]
        self.assertEqual(cert["vertex_count"], 7)
        self.assertEqual(cert["edge_count"], 8)
        self.assertEqual(cert["vertex_orbit_count"], 4)
        self.assertGreaterEqual(cert["edge_orbit_count"], 1)
        self.assertEqual(cert["singleton_vertex_orbits"], [6])
        self.assertEqual(cert["accepted_vertex_localizers"], [])
        self.assertEqual(cert["accepted_edge_localizers"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 64)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_vertex_edge_signatures_built"])
        self.assertTrue(obligations["automorphism_orbit_classification"])
        self.assertFalse(obligations["strict_provenance_available"])
        self.assertFalse(obligations["nonpremise_physical_sector"])
        self.assertFalse(obligations["accepted_vertex_localizer"])
        self.assertFalse(obligations["accepted_edge_localizer"])
        witness = obj["localizer_witness"]
        self.assertEqual(len(witness["vertex_rows"]), 7)
        self.assertEqual(len(witness["edge_rows"]), 8)
        self.assertFalse(any(r["accepted_vertex_localizer"] for r in witness["vertex_rows"]))
        self.assertFalse(any(r["accepted_edge_localizer"] for r in witness["edge_rows"]))

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2998/S1948", MD.read_text(encoding="utf-8"))
        self.assertIn("P2998/S1948", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2998/S1948", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2998", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
