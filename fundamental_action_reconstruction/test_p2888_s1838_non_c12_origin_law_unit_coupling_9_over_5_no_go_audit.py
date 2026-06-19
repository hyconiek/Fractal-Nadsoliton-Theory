import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2888_s1838_non_c12_origin_law_unit_coupling_9_over_5_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2888_s1838_non_c12_origin_law_unit_coupling_9_over_5_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2888_s1838_non_c12_origin_law_unit_coupling_9_over_5_no_go_audit.md"


class P2888NonC12OriginLawUnitCouplingNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["non_c12_origin_law_unit_coupling_9_over_5_no_go_audit"]
        cls.support = cls.audit["support_origin_audit"]
        cls.coupling = cls.audit["unit_coupling_audit"]
        cls.facts = cls.payload["acceptance_matrix"]["facts"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2888_NON_C12_ORIGIN_LAW_UNIT_COUPLING_9_OVER_5_NO_GO_AUDIT_NO_CLOSURE")
        self.assertEqual(self.audit["input_status_rechecked"], "P2887_C12_UNIT_MEASURE_LOCALIZED_DENSITY_9_OVER_5_NO_GO_AUDIT_NO_CLOSURE")
        self.assertTrue(self.facts["p2887_rechecked"])

    def test_support_origin_enumeration(self):
        self.assertEqual(self.support["nonempty_binary_support_count"], 2**12 - 1)
        self.assertEqual(self.support["supports_with_unique_intrinsic_distance_profile_origin"], 360)
        self.assertTrue(self.facts["all_nonempty_binary_supports_checked"])
        self.assertTrue(self.facts["non_c12_unique_origins_exist"])

    def test_unit_coupling_representation_but_nonunique(self):
        self.assertEqual(self.coupling["bounded_density_value_set"], [-2, -1, 0, 1, 2])
        self.assertEqual(self.coupling["candidate_density_assignments_on_unique_origin_supports"], 689443560)
        self.assertEqual(self.coupling["target_9_over_5_record_count"], 600)
        self.assertEqual(self.coupling["target_9_over_5_support_count"], 120)
        self.assertEqual(self.coupling["target_9_over_5_record_count_by_support_size"], {"5": 600})
        self.assertTrue(self.facts["target_9_over_5_representable_in_bounded_non_c12_family"])
        self.assertTrue(self.facts["target_remains_nonunique"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_non_c12_strict_origin_source_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_unit_bearing_9_over_5_coupling_source"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_variational_chain_rule_to_ltotal"])

    def test_documents_updated(self):
        self.assertIn("P2888/S1838", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2888/S1838", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2888/S1838", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2888", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
