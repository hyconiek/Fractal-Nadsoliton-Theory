import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3027_s1977_dissipation_external_unit_source_acceptance_gate.py"
OUT = ROOT / "generated" / "p3027_s1977_dissipation_external_unit_source_acceptance_gate.json"
MD = ROOT / "generated" / "p3027_s1977_dissipation_external_unit_source_acceptance_gate.md"

class P3027DissipationExternalUnitSourceAcceptanceGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3027_DISSIPATION_EXTERNAL_UNIT_SOURCE_ACCEPTANCE_GATE_NO_ACCEPTED_SOURCE")
        self.assertIsNotNone(self.payload["input_hashes"]["P3026"])

    def test_acceptance_gate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["candidate_count"], 4)
        self.assertEqual(cert["obligation_count"], 6)
        self.assertEqual(cert["accepted_candidate_count"], 0)
        self.assertEqual(cert["best_pass_count"], 3)
        self.assertFalse(cert["imported_symbol_accepted"])

    def test_candidates_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "DissipationExternalUnitSource_AcceptanceGate")
        rows = {row["candidate"]: row for row in obj["candidate_rows"]}
        self.assertFalse(rows["P3023_monotone_dissipation_chain"]["accepted"])
        self.assertFalse(rows["P3024_chart_orbit_anchors"]["accepted"])
        self.assertFalse(rows["P3025_internal_tick_H_equals_S_over_tau_ratios"]["accepted"])
        self.assertFalse(rows["formal_imported_unit_symbol_U_ext"]["accepted"])
        self.assertIn("strict_nadsoliton_provenance", rows["formal_imported_unit_symbol_U_ext"]["failed_obligations"])
        self.assertIn("explicit_hamiltonian_coupling", rows["formal_imported_unit_symbol_U_ext"]["failed_obligations"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3027/S1977", MD.read_text(encoding="utf-8"))
        self.assertIn("P3027/S1977", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3027/S1977", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3027", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
