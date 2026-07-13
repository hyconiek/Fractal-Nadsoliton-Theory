import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3153_s2103_axiom_gr_eh_nonproxy_coupling_audit.py"
OUT = ROOT / "generated" / "p3153_s2103_axiom_gr_eh_nonproxy_coupling_audit.json"
MD = ROOT / "generated" / "p3153_s2103_axiom_gr_eh_nonproxy_coupling_audit.md"


class P3153AxiomGrEhNonproxyCouplingAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_finite_counts_show_source_gap(self):
        self.assertEqual(self.payload["status"], "P3153_AXIOM_GR_EH_NONPROXY_COUPLING_AUDIT_RESIDUAL_SOURCE_GAP")
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["frw_metric_candidates"], 5)
        self.assertEqual(counts["vacuum_zero_rows"], 1)
        self.assertEqual(counts["nonflat_rows_needing_source"], 4)
        self.assertEqual(counts["accepted_full_eh_coupling_rows"], 0)

    def test_residual_rows_are_exact_and_nonflat_need_source(self):
        rows = {row["candidate_metric"]: row for row in self.payload["frw_einstein_residual_rows"]}
        self.assertTrue(rows["Minkowski_static"]["vacuum_EH_residual_zero"])
        self.assertEqual(rows["radiation_like_power_p_1_2"]["G_00"], "3/(4*t**2)")
        self.assertEqual(rows["matter_like_power_p_2_3"]["G_00"], "4/(3*t**2)")
        self.assertTrue(all(row["source_needed_for_nonflat_solution"] for name, row in rows.items() if name != "Minkowski_static"))

    def test_no_closure_and_docs_updated(self):
        self.assertFalse(any(row["metric_bundle_exported"] and row["stress_energy_tensor_exported"] and row["newton_or_action_unit_exported"] and row["nonproxy_variation_exported"] and row["noncircular_strict_source"] for row in self.payload["source_interface_rows"]))
        self.assertIn("P3154", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))
        self.assertIn("P3153/S2103", MD.read_text(encoding="utf-8"))
        self.assertIn("P3153/S2103", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3153/S2103", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3153", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
