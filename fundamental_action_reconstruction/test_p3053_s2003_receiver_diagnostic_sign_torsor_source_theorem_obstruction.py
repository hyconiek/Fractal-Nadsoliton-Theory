import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3053_s2003_receiver_diagnostic_sign_torsor_source_theorem_obstruction.py"
OUT = ROOT / "generated" / "p3053_s2003_receiver_diagnostic_sign_torsor_source_theorem_obstruction.json"
MD = ROOT / "generated" / "p3053_s2003_receiver_diagnostic_sign_torsor_source_theorem_obstruction.md"

class P3053ReceiverDiagnosticSignTorsorSourceTheoremObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3053_RECEIVER_DIAGNOSTIC_SIGN_TORSOR_SOURCE_THEOREM_OBSTRUCTION_BOUNDED_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3052"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["diagnostic_signs"], 3)
        self.assertEqual(cert["domain_states"], 8)
        self.assertEqual(cert["boolean_laws_enumerated"], 256)
        self.assertEqual(cert["invariant_laws"], 16)
        self.assertEqual(cert["odd_equivariant_laws"], 16)
        self.assertEqual(cert["odd_polarity_pairs"], 8)
        self.assertEqual(cert["invariant_laws_distinguishing_base_pair"], 0)
        self.assertEqual(cert["artifact_selected_odd_polarity_laws"], 0)
        self.assertEqual(cert["source_acceptance_criteria"], 8)
        self.assertEqual(cert["satisfied_source_acceptance_criteria"], 4)
        self.assertFalse(cert["strict_receiver_diagnostic_source_theorem_exported"])

    def test_obligations_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "ReceiverDiagnosticSignTorsor_SourceTheoremObstruction")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["construct_sign_torsor_source_theorem_object"])
        self.assertTrue(obligations["exhaust_boolean_selector_laws"])
        self.assertTrue(obligations["prove_invariant_laws_cannot_select_orientation"])
        self.assertFalse(obligations["export_unique_nonpremise_polarity"])
        self.assertFalse(obligations["selector_ltotal_bridge_toe_installation"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3053/S2003", MD.read_text(encoding="utf-8"))
        self.assertIn("P3053/S2003", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3053/S2003", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3053", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
