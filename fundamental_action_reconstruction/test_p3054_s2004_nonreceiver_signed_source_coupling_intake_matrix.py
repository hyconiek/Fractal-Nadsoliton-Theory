import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3054_s2004_nonreceiver_signed_source_coupling_intake_matrix.py"
OUT = ROOT / "generated" / "p3054_s2004_nonreceiver_signed_source_coupling_intake_matrix.json"
MD = ROOT / "generated" / "p3054_s2004_nonreceiver_signed_source_coupling_intake_matrix.md"

class P3054NonReceiverSignedSourceCouplingIntakeMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3054_NONRECEIVER_SIGNED_SOURCE_COUPLING_INTAKE_MATRIX_BOUNDED_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3053"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_patterns"], 2)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["candidate_classes"], 5)
        self.assertEqual(cert["accepted_candidate_classes"], 0)
        self.assertEqual(cert["p3053_odd_polarity_pairs"], 8)
        self.assertEqual(cert["coupling_pushout_rows"], 32)
        self.assertEqual(cert["artifact_selected_pushout_rows"], 0)
        self.assertEqual(cert["proof_obligations"], 6)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)
        self.assertFalse(cert["strict_nonreceiver_signed_source_coupling_law_exported"])

    def test_objects_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "NonReceiverSignedSource_CouplingIntakeMatrix")
        self.assertEqual(len(obj["coupling_pushout_rows"]), 32)
        self.assertFalse(any(row["accepted"] for row in obj["candidate_rows"]))
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["construct_coupling_pushout"])
        self.assertTrue(obligations["exhaust_p3053_odd_pair_couplings"])
        self.assertFalse(obligations["export_concrete_nonreceiver_signed_source"])
        self.assertFalse(obligations["install_selector_readout_or_ltotal"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3054/S2004", MD.read_text(encoding="utf-8"))
        self.assertIn("P3054/S2004", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3054/S2004", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3054", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
