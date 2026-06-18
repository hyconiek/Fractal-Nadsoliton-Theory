import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2886_SCRIPT = ROOT / "p2886_s1836_external_unit_measure_action_density_inventory_no_go_audit.py"
SCRIPT = ROOT / "p2887_s1837_c12_unit_measure_localized_density_9_over_5_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2887_s1837_c12_unit_measure_localized_density_9_over_5_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2887_s1837_c12_unit_measure_localized_density_9_over_5_no_go_audit.md"


class P2887C12UnitMeasureLocalizedDensityNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2886_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["c12_unit_measure_localized_density_9_over_5_no_go_audit"]
        cls.support = cls.audit["support_audit"]
        cls.density = cls.audit["density_audit"]
        cls.facts = cls.payload["acceptance_matrix"]["facts"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2887_C12_UNIT_MEASURE_LOCALIZED_DENSITY_9_OVER_5_NO_GO_AUDIT_NO_CLOSURE")
        self.assertEqual(self.audit["input_status_rechecked"], "P2886_EXTERNAL_UNIT_MEASURE_ACTION_DENSITY_INVENTORY_NO_GO_AUDIT_NO_CLOSURE")
        self.assertTrue(self.facts["p2886_rechecked"])

    def test_support_enumeration(self):
        self.assertEqual(self.support["binary_support_count"], 2**12)
        self.assertEqual(self.support["c12_invariant_support_count"], 2)
        self.assertEqual(self.support["c12_invariant_support_sizes"], [0, 12])
        self.assertEqual(self.support["singleton_support_count"], 12)
        self.assertEqual(self.support["c12_invariant_singleton_support_count"], 0)
        self.assertTrue(self.facts["c12_neutral_supports_are_only_empty_or_full"])

    def test_density_enumeration_and_integrals(self):
        self.assertEqual(self.density["ternary_density_count"], 3**12)
        self.assertEqual(self.density["c12_invariant_density_count"], 3)
        self.assertEqual(self.density["singleton_density_count"], 24)
        self.assertEqual(self.density["c12_invariant_singleton_density_count"], 0)
        self.assertEqual(self.density["uniform_unit_measure_integrals_of_invariant_densities"], ["-1", "0", "1"])
        self.assertEqual(self.density["target_9_over_5_integral_count"], 0)
        self.assertTrue(self.facts["uniform_c12_invariant_density_integral_never_9_over_5"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_strict_unit_measure_localized_density_source"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_unit_bearing_9_over_5_coupling"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_variational_chain_rule_to_ltotal"])

    def test_documents_updated(self):
        self.assertIn("P2887/S1837", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2887/S1837", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2887/S1837", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2887", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
