import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2890_s1840_fourier_phase_source_law_9_over_5_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2890_s1840_fourier_phase_source_law_9_over_5_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2890_s1840_fourier_phase_source_law_9_over_5_no_go_audit.md"


class P2890FourierPhaseSourceLawNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["fourier_phase_source_law_9_over_5_no_go_audit"]
        cls.power = cls.audit["phase_blind_power_audit"]
        cls.phase = cls.audit["phaseful_character_boundary"]
        cls.facts = cls.payload["acceptance_matrix"]["facts"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2890_FOURIER_PHASE_SOURCE_LAW_9_OVER_5_NO_GO_AUDIT_NO_CLOSURE")
        self.assertEqual(self.audit["input_status_rechecked"], "P2889_TRANSLATION_ORBIT_SOURCE_LAW_9_OVER_5_NO_GO_AUDIT_NO_CLOSURE")
        self.assertTrue(self.facts["p2889_rechecked"])

    def test_phase_blind_power_obstruction(self):
        self.assertEqual(self.power["translation_orbit_representative_count"], 50)
        self.assertEqual(self.power["phase_blind_autocorrelation_signature_count"], 29)
        self.assertEqual(self.power["phase_blind_signature_multiplicity_histogram"], {"1": 10, "2": 18, "4": 1})
        self.assertFalse(self.power["phase_blind_signatures_unique_for_all_orbits"])
        self.assertTrue(self.facts["phase_blind_power_not_unique_for_all_orbits"])

    def test_phaseful_character_boundary(self):
        self.assertEqual(self.phase["nontrivial_character_mode_count"], 11)
        self.assertEqual(len(self.phase["character_rows"]), 11)
        self.assertTrue(self.phase["all_nontrivial_characters_need_phase_pin_or_have_nontrivial_translation_phase_orbit"])
        self.assertTrue(self.facts["phaseful_characters_need_external_phase_pin"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_strict_translation_breaking_phase_source"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_unique_9_over_5_carrier_orbit_or_representative"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_nonimported_9_over_5_variational_chain_rule"])

    def test_documents_updated(self):
        self.assertIn("P2890/S1840", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2890/S1840", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2890/S1840", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2890", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
