import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3045_s1995_memory_lag_localizer_equivariance_obstruction.py"
OUT = ROOT / "generated" / "p3045_s1995_memory_lag_localizer_equivariance_obstruction.json"
MD = ROOT / "generated" / "p3045_s1995_memory_lag_localizer_equivariance_obstruction.md"

class P3045MemoryLagLocalizerEquivarianceObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3045_MEMORY_LAG_LOCALIZER_EQUIVARIANCE_OBSTRUCTION_BOUNDED_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3044"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["oriented_lag_rows"], 11)
        self.assertEqual(cert["translation_invariant_rows"], 11)
        self.assertEqual(cert["inversion_flip_rows"], 11)
        self.assertEqual(cert["localizer_candidate_rows"], 4)
        self.assertEqual(cert["accepted_chart_independent_lag_localizer_rows"], 0)
        self.assertEqual(cert["aut_compatible_candidate_rows"], 0)
        self.assertFalse(cert["chart_independent_lag_localizer_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "MemoryLagLocalizer_EquivarianceObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["p3044_commutator_used_without_replay"])
        self.assertTrue(obligations["oriented_lag_torsor_constructed"])
        self.assertTrue(obligations["translation_origin_invariance_checked"])
        self.assertTrue(obligations["inversion_action_checked"])
        self.assertFalse(obligations["aut_compatible_nonpremise_lag_localizer"])
        self.assertFalse(obligations["strict_source_law_for_lag"])
        self.assertFalse(obligations["selector_or_readout_coupling"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3045/S1995", MD.read_text(encoding="utf-8"))
        self.assertIn("P3045/S1995", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3045/S1995", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3045", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
