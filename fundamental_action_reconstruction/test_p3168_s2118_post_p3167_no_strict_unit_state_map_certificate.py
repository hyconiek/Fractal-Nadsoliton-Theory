import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3168_s2118_post_p3167_no_strict_unit_state_map_certificate.py"
OUT = ROOT / "generated" / "p3168_s2118_post_p3167_no_strict_unit_state_map_certificate.json"
MD = ROOT / "generated" / "p3168_s2118_post_p3167_no_strict_unit_state_map_certificate.md"


class P3168PostP3167NoStrictUnitStateMapCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_certificate_counts(self):
        self.assertEqual(self.payload["status"], "P3168_POST_P3167_NO_STRICT_UNIT_NO_NEW_LIVE_FRONTIER_CERTIFICATE")
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["replay_classes_checked"], 8)
        self.assertEqual(cert["strict_source_obligations"], 8)
        self.assertEqual(cert["obligation_rows"], 64)
        self.assertEqual(cert["accepted_strict_unit_or_origin_sources"], 0)
        self.assertFalse(cert["strict_S_plus_exported"])
        self.assertFalse(cert["strict_Lambda_origin_exported"])
        self.assertTrue(cert["conditional_CA_SA_schema_exported"])

    def test_state_map_preserves_negative_flags(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertEqual(objs["accepted_rows"], [])
        self.assertIn("W0", objs["minimal_CA_SA_architecture_schema"])
        self.assertIn("CA", objs["minimal_CA_SA_architecture_schema"])
        self.assertIn("SA", objs["minimal_CA_SA_architecture_schema"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))
        self.assertIn("nonzero scale-charged strict S_+", self.payload["decision"]["next_honest_step"])
        self.assertIn("translation-breaking strict Lambda_origin", self.payload["decision"]["next_honest_step"])

    def test_inputs_and_docs_updated(self):
        statuses = self.payload["input_statuses"]
        self.assertIn("P3167_S_PLUS", statuses["P3167"])
        self.assertIn("P3140_AXIOM", statuses["P3140"])
        self.assertIn("P3146_AXIOM", statuses["P3146"])
        self.assertIn("P3168/S2118", MD.read_text(encoding="utf-8"))
        self.assertIn("P3168/S2118", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3168/S2118", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3168", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
