import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2874_SCRIPT = ROOT / "p2874_s1824_dihedral_character_chiral_endpoint_source_no_go_audit.py"
SCRIPT = ROOT / "p2875_s1825_d12_two_dimensional_irrep_endpoint_localization_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2875_s1825_d12_two_dimensional_irrep_endpoint_localization_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2875_s1825_d12_two_dimensional_irrep_endpoint_localization_no_go_audit.md"


class P2875D12TwoDimensionalIrrepEndpointLocalizationNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2874_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["d12_two_dimensional_irrep_endpoint_localization_no_go_audit"]

    def test_status_and_p2874_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2875_D12_TWO_DIMENSIONAL_IRREP_ENDPOINT_LOCALIZATION_NO_GO_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2874_DIHEDRAL_CHARACTER_CHIRAL_ENDPOINT_SOURCE_NO_GO_AUDIT_NO_CLOSURE",
        )

    def test_all_five_two_dimensional_modes_are_enumerated(self):
        self.assertEqual(self.audit["mode_count"], 5)
        self.assertEqual([record["mode"] for record in self.audit["mode_records"]], [1, 2, 3, 4, 5])
        facts = self.payload["acceptance_matrix"]["facts"]
        self.assertTrue(facts["all_five_real_two_dimensional_d12_irrep_modes_enumerated"])
        self.assertTrue(facts["each_mode_has_reflection_fixed_axis_dimension_one"])

    def test_each_mode_is_global_not_singleton_11(self):
        for record in self.audit["mode_records"]:
            self.assertEqual(record["support"], list(range(12)))
            self.assertEqual(record["support_size"], 12)
            self.assertFalse(record["singleton_11_support"])
            self.assertAlmostEqual(record["target_11_norm_squared"], 1.0, places=9)
        facts = self.payload["acceptance_matrix"]["facts"]
        self.assertTrue(facts["each_mode_has_full_endpoint_support"])
        self.assertTrue(facts["no_single_mode_has_singleton_11_support"])

    def test_full_dft_reconstructs_delta_11_only_as_imported_representability(self):
        self.assertEqual(self.audit["delta_11_reconstruction"], [0.0] * 11 + [1.0])
        coefficients = self.audit["delta_11_dft_coefficients"]
        self.assertEqual(len(coefficients), 12)
        self.assertTrue(all(coefficient["abs"] == round(1 / 12, 12) for coefficient in coefficients))
        facts = self.payload["acceptance_matrix"]["facts"]
        self.assertTrue(facts["full_dft_reconstructs_delta_11_with_imported_phase_coefficients"])
        self.assertTrue(facts["dft_coefficients_depend_on_imported_target_11"])

    def test_no_false_exports_and_documents_updated(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["boundary_source_law_exported"])
        self.assertFalse(flags["orientation_source_law_exported"])
        self.assertFalse(flags["chiral_endpoint_11_source_law_exported"])
        self.assertFalse(flags["selector_or_localizer_source_exported"])
        self.assertFalse(flags["unit_bearing_coupling_localization_theorem_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2875/S1825", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2875/S1825", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2875/S1825", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2875", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
