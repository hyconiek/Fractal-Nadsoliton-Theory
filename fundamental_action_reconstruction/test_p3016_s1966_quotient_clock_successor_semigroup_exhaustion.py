import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3016_s1966_quotient_clock_successor_semigroup_exhaustion.py"
OUT = ROOT / "generated" / "p3016_s1966_quotient_clock_successor_semigroup_exhaustion.json"
MD = ROOT / "generated" / "p3016_s1966_quotient_clock_successor_semigroup_exhaustion.md"

class P3016QuotientClockSuccessorSemigroupTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3016_QUOTIENT_CLOCK_SUCCESSOR_SEMIGROUP_EXHAUSTION_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3015"])

    def test_exhaustion_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["orbit_count"], 6)
        self.assertEqual(cert["label_constraint_count"], 12)
        self.assertEqual(cert["conflicting_source_orbit_count"], 4)
        self.assertEqual(cert["total_candidate_maps_exhausted"], 46656)
        self.assertEqual(cert["max_satisfied_label_constraints"], 6)
        self.assertEqual(cert["required_label_constraints"], 12)
        self.assertFalse(cert["accepted_as_directed_successor_semigroup"])

    def test_obligations_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "QuotientClockSuccessorSemigroup_ExhaustionNoGoMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["unit_quotient_domain"])
        self.assertTrue(obligations["total_successor_map_on_quotient"])
        self.assertFalse(obligations["represents_all_label_successors"])
        self.assertFalse(obligations["nonconflicting_source_orbit_targets"])
        self.assertFalse(obligations["directed_semigroup_time_evolution_export"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3016/S1966", MD.read_text(encoding="utf-8"))
        self.assertIn("P3016/S1966", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3016/S1966", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3016", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
