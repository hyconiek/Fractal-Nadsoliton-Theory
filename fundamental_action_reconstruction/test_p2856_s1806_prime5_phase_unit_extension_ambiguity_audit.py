import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2855_SCRIPT = ROOT / "p2855_s1805_z12_rational_phase_lattice_source_candidate_audit.py"
SCRIPT = ROOT / "p2856_s1806_prime5_phase_unit_extension_ambiguity_audit.py"
JSON_PATH = ROOT / "generated" / "p2856_s1806_prime5_phase_unit_extension_ambiguity_audit.json"
MD_PATH = ROOT / "generated" / "p2856_s1806_prime5_phase_unit_extension_ambiguity_audit.md"


class P2856Prime5PhaseUnitExtensionAmbiguityAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2855_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["prime5_phase_unit_extension_ambiguity_audit"]

    def test_status_and_p2855_input(self):
        self.assertEqual(self.payload["status"], "P2856_PRIME5_PHASE_UNIT_EXTENSION_AMBIGUITY_AUDIT_NO_CLOSURE")
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2855_Z12_RATIONAL_PHASE_LATTICE_SOURCE_CANDIDATE_AUDIT_NO_CLOSURE",
        )

    def test_prime5_extension_represents_strict_tuple(self):
        self.assertEqual(self.audit["strict_omega"]["fraction"], "743/4000")
        self.assertEqual(self.audit["strict_phi"]["fraction"], "13/80")
        self.assertIn(4000, self.audit["omega_exact_supporting_denominators_first12"])
        self.assertIn(80, self.audit["phi_exact_supporting_denominators_first12"])
        facts = self.payload["acceptance_matrix"]["facts"]
        self.assertTrue(facts["prime5_extension_represents_omega"])
        self.assertTrue(facts["prime5_extension_represents_phi"])

    def test_local_phase_bit_ambiguity_blocks_selection(self):
        witnesses = self.audit["same_phase_bit_local_witnesses"]
        self.assertGreater(len(witnesses), 0)
        for row in witnesses:
            self.assertNotEqual((row["omega"]["fraction"], row["phi"]["fraction"]), ("743/4000", "13/80"))
        facts = self.payload["acceptance_matrix"]["facts"]
        self.assertTrue(facts["local_ambiguity_witness_exists"])

    def test_no_false_exports_and_documents_updated(self):
        self.assertEqual(self.audit["accepted_candidate_count"], 0)
        self.assertTrue(self.payload["acceptance_matrix"]["accepted_as_prime5_extension_ambiguity_audit"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_prime5_phase_unit_source_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_strict_phase_frequency_source_law"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["prime5_phase_unit_source_exported"])
        self.assertFalse(flags["strict_phase_frequency_source_law_exported"])
        self.assertFalse(flags["full_kernel_bridge_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2856/S1806", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2856/S1806", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2856/S1806", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2856", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
