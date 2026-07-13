import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3148_s2098_sm_representation_registry_completion_audit.py"
OUT = ROOT / "generated" / "p3148_s2098_sm_representation_registry_completion_audit.json"
MD = ROOT / "generated" / "p3148_s2098_sm_representation_registry_completion_audit.md"


class P3148SmRepresentationRegistryCompletionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_counts(self):
        self.assertEqual(self.payload["status"], "P3148_SM_REPRESENTATION_REGISTRY_COMPLETION_ALGEBRAIC_PASS_CONDITIONAL_NO_STRICT_SOURCE")
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["registry_rows"], 6)
        self.assertEqual(counts["fermion_rows"], 5)
        self.assertEqual(counts["higgs_rows"], 1)
        self.assertEqual(counts["representation_rows_passing"], 6)
        self.assertEqual(counts["total_representation_failure_slots"], 0)
        self.assertEqual(counts["anomaly_rows_vanishing"], 4)
        self.assertEqual(counts["yukawa_rows_u1_invariant"], 3)
        self.assertEqual(counts["strict_nadsoliton_source_rows"], 0)
        self.assertEqual(counts["unit_bearing_ltotal_rows"], 0)

    def test_registry_and_anomalies(self):
        fields = {row["field"] for row in self.payload["registry_rows"]}
        self.assertEqual(fields, {"Q_L", "u_c", "d_c", "L_L", "e_c", "H"})
        self.assertTrue(all(row["all_representation_checks_pass"] for row in self.payload["representation_commutator_rows"]))
        self.assertTrue(all(row["vanishes"] for row in self.payload["anomaly_sums"].values()))
        self.assertTrue(all(row["u1_invariant"] for row in self.payload["yukawa_hypercharge_checks"].values()))

    def test_decision_preserves_strict_boundaries(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["accepted_as_algebraic_registry_completion"])
        self.assertIn("algebraic registry", decision["bounded_result"])
        self.assertIn("SM ansatz", decision["why_not_strict"])
        self.assertIn("P3149", decision["next_honest_step"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_updated(self):
        self.assertIn("P3148/S2098", MD.read_text(encoding="utf-8"))
        self.assertIn("P3148/S2098", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3148/S2098", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3148", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
