from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2711_s1661_inversion_odd_character_source_law_sign_coupling_audit.py"
OUT = ROOT / "generated" / "p2711_s1661_inversion_odd_character_source_law_sign_coupling_audit.json"
MD = ROOT / "generated" / "p2711_s1661_inversion_odd_character_source_law_sign_coupling_audit.md"


class P2711InversionOddCharacterSourceLawSignCouplingAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_finite_source_law_candidates_are_enumerated(self):
        counts = self.payload["finite_counts"]
        self.assertEqual(counts["anti_inversion_characters"], 2)
        self.assertEqual(counts["source_law_candidates"], 4)
        self.assertEqual(counts["sign_degenerate_pairs"], 2)

    def test_lambda_sign_degeneracy_remains(self):
        for row in self.payload["sign_pair_degeneracy"]:
            self.assertEqual(row["lambda_values"], [-1, 1])
            self.assertTrue(row["strict_sign_source_required"])
            self.assertFalse(row["strict_sign_source_exported"])
            self.assertIn("lambda -> -lambda", row["degeneracy"])

    def test_no_strict_source_law_or_closure_exported(self):
        self.assertEqual(self.payload["status"], "P2711_SOURCE_LAW_SIGN_COUPLING_DEGENERACY_NO_STRICT_SIGN_SOURCE")
        self.assertTrue(all(not row["exports_nonpremise_coupling_sign"] for row in self.payload["artifact_sign_source_scan"]))
        decision = self.payload["decision"]
        self.assertTrue(decision["finite_source_law_candidates_enumerated"])
        self.assertFalse(decision["source_law_coupling_sign_exported"])
        self.assertFalse(decision["strict_source_law_selects_plus_omega"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_outputs_and_docs_written(self):
        self.assertIn("P2711/S1661", MD.read_text(encoding="utf-8"))
        self.assertIn("P2711/S1661", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2711/S1661", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2711", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
