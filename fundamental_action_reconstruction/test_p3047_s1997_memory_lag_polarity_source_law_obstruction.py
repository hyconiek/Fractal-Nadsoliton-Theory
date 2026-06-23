import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3047_s1997_memory_lag_polarity_source_law_obstruction.py"
OUT = ROOT / "generated" / "p3047_s1997_memory_lag_polarity_source_law_obstruction.json"
MD = ROOT / "generated" / "p3047_s1997_memory_lag_polarity_source_law_obstruction.md"

class P3047MemoryLagPolaritySourceLawObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3047_MEMORY_LAG_POLARITY_SOURCE_LAW_OBSTRUCTION_BOUNDED_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3046"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["source_domain_rows"], 5)
        self.assertEqual(cert["trivial_source_rows"], 4)
        self.assertEqual(cert["trivial_source_equivariant_maps"], 0)
        self.assertEqual(cert["odd_source_equivariant_maps"], 2)
        self.assertEqual(cert["candidate_law_rows"], 3)
        self.assertEqual(cert["accepted_polarity_source_law_rows"], 0)
        self.assertFalse(cert["strict_polarity_source_law_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "MemoryLagPolarity_SourceLawObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["p3046_coupling_polarity_gap_used"])
        self.assertTrue(obligations["trivial_source_no_section_theorem"])
        self.assertTrue(obligations["odd_source_representation_identified"])
        self.assertFalse(obligations["nonzero_strict_odd_source_value"])
        self.assertFalse(obligations["unique_polarity_selection_law"])
        self.assertFalse(obligations["selector_readout_installation"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3047/S1997", MD.read_text(encoding="utf-8"))
        self.assertIn("P3047/S1997", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3047/S1997", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3047", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
