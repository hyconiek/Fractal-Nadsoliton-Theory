from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2722_s1672_clock_phase_polarity_selection_source_law_audit.py"
OUT = ROOT / "generated" / "p2722_s1672_clock_phase_polarity_selection_source_law_audit.json"
MD = ROOT / "generated" / "p2722_s1672_clock_phase_polarity_selection_source_law_audit.md"


class P2722ClockPhasePolaritySelectionSourceLawAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_clock_phase_candidates_do_not_select_strict_polarity(self):
        self.assertEqual(self.payload["status"], "P2722_CLOCK_PHASE_POLARITY_SOURCE_LAW_AUDIT_NO_UNLOCK")
        audit = self.payload["clock_phase_audit"]
        self.assertEqual(audit["accepted_candidate_count"], 0)
        self.assertTrue(audit["p2721_conditional_pair_reused"])
        self.assertFalse(audit["strict_polarity_selection_source_law_exported"])
        self.assertEqual(len(audit["candidate_rows"]), 4)
        self.assertTrue(all(not row["accepted_as_strict_polarity_source"] for row in audit["candidate_rows"]))

    def test_moving_clock_hand_is_not_a_nonpremise_origin(self):
        rows = {row["name"]: row for row in self.payload["clock_phase_audit"]["candidate_rows"]}
        self.assertIn("origin-relative", rows["moving_clock_hand_position_theta_t"]["blocker"])
        self.assertIn("external face-label convention", rows["fixed_clock_face_zero"]["blocker"])
        self.assertIn("chiral/time-orientation datum", rows["angular_velocity_sign_dtheta_dt"]["blocker"])
        self.assertIn("conditional pair", rows["p2721_equivariant_coupling_pair"]["blocker"])

    def test_no_closure_flags_are_exported(self):
        decision = self.payload["decision"]
        self.assertFalse(decision["strict_polarity_selection_source_law_exported"])
        self.assertFalse(decision["clock_phase_motion_source_exported"])
        self.assertFalse(decision["continuous_phase_origin_exported"])
        self.assertFalse(decision["canonical_coupling_polarity_selected"])
        self.assertFalse(decision["strict_mechanism_fixing_lambda_exported"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_written(self):
        self.assertIn("P2722/S1672", MD.read_text(encoding="utf-8"))
        self.assertIn("P2722/S1672", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2722/S1672", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2722", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
