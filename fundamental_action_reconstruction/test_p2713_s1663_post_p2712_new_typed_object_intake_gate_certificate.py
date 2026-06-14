from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2713_s1663_post_p2712_new_typed_object_intake_gate_certificate.py"
OUT = ROOT / "generated" / "p2713_s1663_post_p2712_new_typed_object_intake_gate_certificate.json"
MD = ROOT / "generated" / "p2713_s1663_post_p2712_new_typed_object_intake_gate_certificate.md"


class P2713PostP2712NewTypedObjectIntakeGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_no_candidate_certificate_preserves_prior_no_frontier(self):
        self.assertEqual(self.payload["status"], "P2713_NEW_TYPED_OBJECT_INTAKE_GATE_NO_CANDIDATE_CERTIFICATE")
        self.assertTrue(self.payload["prior_certificates_hold"])
        self.assertEqual(self.payload["admitted_candidate_count"], 0)
        self.assertTrue(self.payload["decision"]["no_new_live_frontier_certificate_preserved"])
        self.assertTrue(self.payload["decision"]["replay_lanes_frozen"])

    def test_intake_gate_requires_new_object_and_keeps_closure_flags_false(self):
        intake = self.payload["candidate_intake"]
        self.assertEqual(len(intake), 1)
        self.assertEqual(intake[0]["candidate_id"], "NO_NEW_TYPED_OBJECT_SUPPLIED")
        self.assertFalse(intake[0]["admitted_for_next_test"])
        criteria = {row["criterion"] for row in self.payload["acceptance_criteria"]}
        self.assertIn("strict_source_or_typed_object_is_new", criteria)
        self.assertIn("non_premise_non_convention_export", criteria)
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_closed_replay_lanes_and_docs_written(self):
        lanes = set(self.payload["closed_replay_lanes"])
        self.assertIn("selector_sign_replay_or_qw2191_discharge_without_new_source", lanes)
        self.assertIn("ltotal_or_toe_promotion", lanes)
        self.assertIn("P2713/S1663", MD.read_text(encoding="utf-8"))
        self.assertIn("P2713/S1663", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2713/S1663", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2713", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
