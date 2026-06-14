import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2731_s1681_local_time_arrow_variational_law_no_go.py"
OUT = ROOT / "generated" / "p2731_s1681_local_time_arrow_variational_law_no_go.json"
MD = ROOT / "generated" / "p2731_s1681_local_time_arrow_variational_law_no_go.md"


class P2731LocalTimeArrowVariationalLawNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.audit = cls.payload["local_variational_time_arrow_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_exhausts_local_pair_energy_laws_and_fields(self):
        self.assertEqual(self.payload["status"], "P2731_LOCAL_TIME_ARROW_VARIATIONAL_LAW_NO_GO")
        self.assertEqual(self.audit["local_pair_energy_table_count"], 81)
        self.assertEqual(self.audit["field_count"], 4096)
        self.assertEqual(self.audit["time_reversal_even_law_count"], 9)

    def test_even_laws_have_no_unpaired_tau_selector(self):
        self.assertEqual(self.audit["time_reversal_even_unpaired_selector_count"], 0)
        self.assertGreater(self.audit["non_time_reversal_even_unpaired_selector_count"], 0)
        self.assertIn("Every time-reversal-even", self.audit["obstruction"])

    def test_acceptance_matrix_blocks_promotion(self):
        self.assertFalse(self.acceptance["accepted_as_strict_time_arrow_source_law"])
        self.assertIn("time_reversal_even_unpaired_selector_exists", self.acceptance["missing_criteria"])
        self.assertIn("strict_nonpremise_time_arrow_variational_source_exported", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_requires_odd_source_term_not_replay(self):
        recommendation = self.payload["decision"]["next_honest_step"]
        self.assertIn("time-reversal-odd source term", recommendation)
        self.assertIn("P2697-P2731", recommendation)

    def test_documentation_updated(self):
        self.assertIn("P2731/S1681", MD.read_text(encoding="utf-8"))
        self.assertIn("P2731/S1681", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2731/S1681", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2731", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
