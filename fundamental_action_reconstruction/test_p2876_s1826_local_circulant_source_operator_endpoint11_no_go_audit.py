import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2875_SCRIPT = ROOT / "p2875_s1825_d12_two_dimensional_irrep_endpoint_localization_no_go_audit.py"
SCRIPT = ROOT / "p2876_s1826_local_circulant_source_operator_endpoint11_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2876_s1826_local_circulant_source_operator_endpoint11_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2876_s1826_local_circulant_source_operator_endpoint11_no_go_audit.md"


class P2876LocalCirculantSourceOperatorEndpoint11NoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2875_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["local_circulant_source_operator_endpoint11_no_go_audit"]

    def test_status_and_p2875_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2876_LOCAL_CIRCULANT_SOURCE_OPERATOR_ENDPOINT11_NO_GO_AUDIT_NO_CLOSURE",
        )
        self.assertEqual(
            self.audit["input_status_rechecked"],
            "P2875_D12_TWO_DIMENSIONAL_IRREP_ENDPOINT_LOCALIZATION_NO_GO_AUDIT_NO_CLOSURE",
        )

    def test_stencil_families_are_enumerated(self):
        full = self.audit["full_family_summary"]
        reflected = self.audit["reflection_symmetric_summary"]
        self.assertEqual(full["stencil_count"], 243)
        self.assertEqual(reflected["stencil_count"], 27)
        self.assertEqual(full["operator_input_pair_count"], 243 * 7)
        self.assertEqual(reflected["operator_input_pair_count"], 27 * 7)
        facts = self.payload["acceptance_matrix"]["facts"]
        self.assertTrue(facts["all_ternary_radius2_circulant_stencils_enumerated"])
        self.assertTrue(facts["reflection_symmetric_subfamily_enumerated"])
        self.assertTrue(facts["input_modes_k0_through_k6_checked"])

    def test_no_singleton_11_outputs(self):
        full = self.audit["full_family_summary"]
        reflected = self.audit["reflection_symmetric_summary"]
        self.assertEqual(full["singleton_11_hit_count"], 0)
        self.assertEqual(reflected["singleton_11_hit_count"], 0)
        self.assertEqual(full["singleton_11_hits"], [])
        self.assertEqual(reflected["singleton_11_hits"], [])
        facts = self.payload["acceptance_matrix"]["facts"]
        self.assertTrue(facts["no_full_family_singleton_11_hits"])
        self.assertTrue(facts["no_reflection_symmetric_singleton_11_hits"])

    def test_outputs_are_zero_or_global_on_fourier_inputs(self):
        full = self.audit["full_family_summary"]
        self.assertEqual(full["partial_non_singleton_output_count"], 0)
        self.assertEqual(full["nonzero_global_output_count"] + full["zero_output_count"], full["operator_input_pair_count"])
        self.assertTrue(self.payload["acceptance_matrix"]["facts"]["outputs_are_only_zero_or_global_for_fourier_inputs"])
        self.assertTrue(self.audit["sample_records"])

    def test_no_false_exports_and_documents_updated(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["boundary_source_law_exported"])
        self.assertFalse(flags["orientation_source_law_exported"])
        self.assertFalse(flags["chiral_endpoint_11_source_law_exported"])
        self.assertFalse(flags["selector_or_localizer_source_exported"])
        self.assertFalse(flags["unit_bearing_coupling_localization_theorem_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2876/S1826", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2876/S1826", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2876/S1826", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2876", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
