import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3046_s1996_memory_lag_torsor_coupling_polarity_audit.py"
OUT = ROOT / "generated" / "p3046_s1996_memory_lag_torsor_coupling_polarity_audit.json"
MD = ROOT / "generated" / "p3046_s1996_memory_lag_torsor_coupling_polarity_audit.md"

class P3046MemoryLagTorsorCouplingPolarityAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3046_MEMORY_LAG_TORSOR_COUPLING_POLARITY_AUDIT_BOUNDED_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3045"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["sign_torsor_size"], 2)
        self.assertEqual(cert["orientation_torsor_size"], 2)
        self.assertEqual(cert["candidate_coupling_rows"], 2)
        self.assertEqual(cert["aut_equivariant_coupling_rows"], 2)
        self.assertEqual(cert["polarity_selected_rows"], 0)
        self.assertEqual(cert["accepted_selector_readout_coupling_rows"], 0)
        self.assertEqual(cert["source_readout_rows"], 3)
        self.assertEqual(cert["accepted_source_readout_rows"], 0)
        self.assertFalse(cert["selector_readout_coupling_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "MemoryLagTorsor_CouplingPolarityAudit")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["p3045_signed_lag_torsor_used"])
        self.assertTrue(obligations["finite_equivariant_coupling_pair_enumerated"])
        self.assertTrue(obligations["coupling_polarity_pair_exposed"])
        self.assertFalse(obligations["strict_polarity_selection_law"])
        self.assertFalse(obligations["nonpremise_selector_readout_installation"])
        self.assertFalse(obligations["strict_memory_lag_source_law"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3046/S1996", MD.read_text(encoding="utf-8"))
        self.assertIn("P3046/S1996", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3046/S1996", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3046", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
