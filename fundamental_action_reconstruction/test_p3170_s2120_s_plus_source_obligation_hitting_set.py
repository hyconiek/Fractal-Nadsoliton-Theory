import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3170_s2120_s_plus_source_obligation_hitting_set.py"
OUT = ROOT / "generated" / "p3170_s2120_s_plus_source_obligation_hitting_set.json"
MD = ROOT / "generated" / "p3170_s2120_s_plus_source_obligation_hitting_set.md"


class P3170SPlusSourceObligationHittingSetTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_hitting_set_certificate(self):
        self.assertEqual(self.payload["status"], "P3170_S_PLUS_SOURCE_OBLIGATION_HITTING_SET_THEOREM")
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["candidate_routes"], 12)
        self.assertEqual(cert["distinct_failed_obligations"], 8)
        self.assertEqual(cert["minimal_hitting_sets"], 1)
        self.assertEqual(cert["minimal_hitting_set_size"], 1)
        self.assertEqual(cert["common_failed_atoms_count"], 1)
        self.assertEqual(cert["admissible_schema_missing_atoms"], 2)
        self.assertEqual(cert["accepted_S_plus_sources"], 0)

    def test_minimal_blocker_and_repair_frontier(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertEqual(objs["minimal_universal_blocker_cut"], [["strict_nadsoliton_source_exported"]])
        self.assertEqual(objs["common_failed_atoms"], ["strict_nadsoliton_source_exported"])
        schema_rows = [row for row in objs["least_repair_frontier"] if row["candidate"] == "new_strict_scale_charged_datum_schema"]
        self.assertEqual(len(schema_rows), 1)
        self.assertTrue(schema_rows[0]["recommended_for_repair"])
        self.assertEqual(schema_rows[0]["missing_atoms"], ["nonzero_value_exported", "strict_nadsoliton_source_exported"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_docs_updated_and_next_step(self):
        self.assertIn("P3170/S2120", MD.read_text(encoding="utf-8"))
        self.assertIn("P3170/S2120", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3170/S2120", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3170", (REPO / "AGENTS.md").read_text(encoding="utf-8"))
        self.assertIn("Source_S_plus", self.payload["decision"]["next_honest_step"])
        self.assertIn("s > 0", self.payload["decision"]["next_honest_step"])


if __name__ == "__main__":
    unittest.main()
