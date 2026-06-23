import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3044_s1994_memory_lag_commutator_source_candidate.py"
OUT = ROOT / "generated" / "p3044_s1994_memory_lag_commutator_source_candidate.json"
MD = ROOT / "generated" / "p3044_s1994_memory_lag_commutator_source_candidate.md"

class P3044MemoryLagCommutatorSourceCandidateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3044_MEMORY_LAG_COMMUTATOR_SOURCE_CANDIDATE_BOUNDED_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3043"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["lag_rows"], 6)
        self.assertEqual(cert["finite_nonzero_lag_rows"], 6)
        self.assertEqual(cert["signed_sum_nonzero_lag_rows"], 5)
        self.assertEqual(cert["exchange_antisymmetry_rows"], 6)
        self.assertEqual(cert["accepted_new_strict_source_law_rows"], 0)
        self.assertEqual(cert["p3043_predicates_satisfied"], 2)
        self.assertFalse(cert["new_strict_source_law_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "MemoryLagCommutator_SourceCandidateAcceptanceMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["p3043_new_source_intake_used"])
        self.assertTrue(obligations["outside_exhausted_receiver_classes"])
        self.assertTrue(obligations["finite_nonzero_signed_commutator"])
        self.assertTrue(obligations["exchange_antisymmetry_verified"])
        self.assertFalse(obligations["strict_nadsoliton_source_law"])
        self.assertFalse(obligations["chart_independent_lag_localizer"])
        self.assertFalse(obligations["selector_torsor_or_readout_coupling"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3044/S1994", MD.read_text(encoding="utf-8"))
        self.assertIn("P3044/S1994", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3044/S1994", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3044", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
