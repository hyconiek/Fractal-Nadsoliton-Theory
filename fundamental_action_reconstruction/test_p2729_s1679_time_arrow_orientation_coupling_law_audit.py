from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2729_s1679_time_arrow_orientation_coupling_law_audit.py"
OUT = ROOT / "generated" / "p2729_s1679_time_arrow_orientation_coupling_law_audit.json"
MD = ROOT / "generated" / "p2729_s1679_time_arrow_orientation_coupling_law_audit.md"


class P2729TimeArrowOrientationCouplingLawAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_exhausts_time_arrow_laws(self):
        self.assertEqual(self.payload["status"], "P2729_TIME_ARROW_CONDITIONAL_POLARITY_BUT_NO_STRICT_SOURCE")
        audit = self.payload["time_arrow_audit"]
        self.assertEqual(audit["law_count"], 16)
        self.assertEqual(audit["checked_transition_rows"], 9216)
        self.assertGreater(len(audit["tau_dependent_equivariant_law_names"]), 0)
        self.assertIn("follow_time_arrow_tau", audit["tau_dependent_equivariant_law_names"])

    def test_fixed_tau_is_conditional_not_unconditional(self):
        audit = self.payload["time_arrow_audit"]
        self.assertGreater(len(audit["fixed_tau_polarity_selector_names"]), 0)
        self.assertEqual(audit["unconditional_time_reversal_equivariant_selector_count"], 0)
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["facts"]["fixed_tau_can_conditionally_select_polarity"])
        self.assertFalse(acceptance["facts"]["strict_time_arrow_source_value_exported"])
        self.assertFalse(acceptance["accepted_as_strict_time_arrow_source"])

    def test_no_closure_flags_are_exported(self):
        decision = self.payload["decision"]
        self.assertFalse(decision["strict_time_arrow_source_value_exported"])
        self.assertFalse(decision["time_arrow_selected_as_nonpremise_tau"])
        self.assertFalse(decision["time_arrow_p2721_polarity_coupling_exported"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_written(self):
        self.assertIn("P2729/S1679", MD.read_text(encoding="utf-8"))
        self.assertIn("P2729/S1679", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2729/S1679", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2729", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
