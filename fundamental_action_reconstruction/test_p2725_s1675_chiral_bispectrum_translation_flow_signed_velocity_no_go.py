from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2725_s1675_chiral_bispectrum_translation_flow_signed_velocity_no_go.py"
OUT = ROOT / "generated" / "p2725_s1675_chiral_bispectrum_translation_flow_signed_velocity_no_go.json"
MD = ROOT / "generated" / "p2725_s1675_chiral_bispectrum_translation_flow_signed_velocity_no_go.md"


class P2725ChiralBispectrumTranslationFlowSignedVelocityNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_all_translation_flow_velocities_are_zero(self):
        self.assertEqual(self.payload["status"], "P2725_TRANSLATION_FLOW_SIGNED_VELOCITY_NO_GO_NO_STRICT_DYNAMIC_CHIRAL_SOURCE")
        audit = self.payload["translation_flow_audit"]
        self.assertEqual(audit["checked_flow_rows"], 264)
        self.assertEqual(audit["velocity_count"], 11)
        self.assertEqual(audit["nonzero_signed_velocity_count"], 0)
        self.assertTrue(audit["all_translation_flow_deltas_zero"])
        self.assertTrue(all(row["all_zero"] for row in audit["velocity_summary_rows"]))

    def test_acceptance_matrix_blocks_dynamic_source(self):
        matrix = self.payload["acceptance_matrix"]
        self.assertTrue(matrix["facts"]["explicit_dynamic_lift_defined"])
        self.assertTrue(matrix["facts"]["computable_signed_value_exported"])
        self.assertFalse(matrix["facts"]["nonzero_signed_value_exported"])
        self.assertFalse(matrix["facts"]["coupled_to_p2721_polarity_pair"])
        self.assertFalse(matrix["accepted_as_strict_dynamic_chiral_source"])

    def test_no_closure_flags_are_exported(self):
        decision = self.payload["decision"]
        self.assertFalse(decision["strict_dynamic_chiral_source_artifact_exported"])
        self.assertFalse(decision["p2721_polarity_selected"])
        self.assertFalse(decision["strict_mechanism_fixing_lambda_exported"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_written(self):
        self.assertIn("P2725/S1675", MD.read_text(encoding="utf-8"))
        self.assertIn("P2725/S1675", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2725/S1675", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2725", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
