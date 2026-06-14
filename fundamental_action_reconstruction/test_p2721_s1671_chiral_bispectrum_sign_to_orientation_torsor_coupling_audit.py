from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.py"
OUT = ROOT / "generated" / "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.json"
MD = ROOT / "generated" / "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.md"


class P2721ChiralBispectrumSignToOrientationTorsorCouplingAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_two_equivariant_couplings_exist_but_are_not_unconditional(self):
        self.assertEqual(self.payload["status"], "P2721_CONDITIONAL_SIGN_TORSOR_COUPLINGS_EXIST_BUT_NO_CANONICAL_POLARITY")
        audit = self.payload["coupling_audit"]
        self.assertEqual(audit["candidate_count"], 4)
        self.assertEqual(audit["aut_equivariant_coupling_count"], 2)
        self.assertFalse(audit["accepted_unconditional_coupling"])
        facts = audit["facts"]
        self.assertTrue(facts["marker_sign_torsor_defined"])
        self.assertTrue(facts["aut_equivariant_couplings_exist"])
        self.assertTrue(facts["coupling_polarity_pair_remains"])
        self.assertFalse(facts["phase_origin_source_localizer_exported"])
        self.assertFalse(facts["canonical_coupling_polarity_selected"])

    def test_no_closure_flags_are_exported(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["conditional_equivariant_coupling_exported"])
        self.assertFalse(decision["canonical_coupling_polarity_selected"])
        self.assertFalse(decision["unconditional_torsor_coupling_theorem_exported"])
        self.assertFalse(decision["strict_chiral_bispectrum_source_exported"])
        self.assertFalse(decision["strict_mechanism_fixing_lambda_exported"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_written(self):
        self.assertIn("P2721/S1671", MD.read_text(encoding="utf-8"))
        self.assertIn("P2721/S1671", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2721/S1671", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2721", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
