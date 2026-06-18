import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2877_SCRIPT = ROOT / "p2877_s1827_pinned_non_circulant_defect_origin_import_no_go_audit.py"
SCRIPT = ROOT / "p2878_s1828_d12_equivariant_pinned_defect_origin_law_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2878_s1828_d12_equivariant_pinned_defect_origin_law_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2878_s1828_d12_equivariant_pinned_defect_origin_law_no_go_audit.md"


class P2878D12EquivariantPinnedDefectOriginLawNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2877_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["d12_equivariant_pinned_defect_origin_law_no_go_audit"]

    def test_status_and_p2877_input(self):
        self.assertEqual(self.payload["status"], "P2878_D12_EQUIVARIANT_PINNED_DEFECT_ORIGIN_LAW_NO_GO_AUDIT_NO_CLOSURE")
        self.assertEqual(self.audit["input_status_rechecked"], "P2877_PINNED_NON_CIRCULANT_DEFECT_ORIGIN_IMPORT_NO_GO_AUDIT_NO_CLOSURE")

    def test_group_and_density_family_are_finite_and_checked(self):
        self.assertEqual(self.audit["raw_pinned_defect_vector_count"], 324)
        self.assertGreater(self.audit["unique_density_vector_count"], 0)
        self.assertEqual(self.audit["d12_group_size"], 24)
        self.assertGreater(self.audit["orbit_count"], 0)
        facts = self.payload["acceptance_matrix"]["facts"]
        self.assertTrue(facts["all_pinned_defect_vectors_generated"])
        self.assertTrue(facts["d12_group_size_checked"])

    def test_no_orbit_forces_endpoint_11_or_unique_selector(self):
        self.assertEqual(self.audit["target_11_forcing_orbit_count"], 0)
        self.assertEqual(self.audit["unique_endpoint_selector_orbit_count"], 0)
        facts = self.payload["acceptance_matrix"]["facts"]
        self.assertTrue(facts["no_orbit_forces_endpoint_11"])
        self.assertTrue(facts["no_orbit_has_unique_endpoint_selector"])
        self.assertTrue(facts["singleton_orbits_are_not_endpoint_11_source_laws"])

    def test_orbit_records_carry_stabilizer_fixed_point_certificate(self):
        self.assertTrue(self.audit["sample_orbit_records"])
        for record in self.audit["sample_orbit_records"]:
            self.assertIn("stabilizer_size", record)
            self.assertIn("fixed_endpoint_options_for_equivariant_selector", record)
            self.assertFalse(record["target_11_forced"])
        self.assertIn("equivariant_map_criterion", self.audit["proof_certificate"])

    def test_no_false_exports_and_documents_updated(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["independent_origin_support_law_exported"])
        self.assertFalse(flags["boundary_source_law_exported"])
        self.assertFalse(flags["unit_bearing_9_over_5_coupling_theorem_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2878/S1828", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2878/S1828", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2878/S1828", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2878", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
