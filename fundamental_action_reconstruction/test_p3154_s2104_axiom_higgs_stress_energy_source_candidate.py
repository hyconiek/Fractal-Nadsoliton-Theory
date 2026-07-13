import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3154_s2104_axiom_higgs_stress_energy_source_candidate.py"
OUT = ROOT / "generated" / "p3154_s2104_axiom_higgs_stress_energy_source_candidate.json"
MD = ROOT / "generated" / "p3154_s2104_axiom_higgs_stress_energy_source_candidate.md"


class P3154AxiomHiggsStressEnergySourceCandidateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_symbolic_tensor_and_identity_counts(self):
        self.assertEqual(self.payload["status"], "P3154_AXIOM_HIGGS_STRESS_ENERGY_SYMBOLIC_CANDIDATE_NO_STRICT_SOURCE")
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["symbolic_tensor_components"], 2)
        self.assertEqual(counts["conservation_identity_residual_zero"], 1)
        self.assertEqual(counts["state_candidate_rows"], 3)
        self.assertEqual(counts["source_gate_rows"], 5)
        self.assertEqual(counts["source_gates_satisfied"], 1)
        self.assertEqual(counts["accepted_strict_stress_energy_sources"], 0)

    def test_conservation_identity_and_candidate_rows(self):
        tensor = self.payload["symbolic_stress_energy"]
        self.assertEqual(tensor["conservation_minus_hdot_times_kg"], "0")
        rows = {row["state_candidate"]: row for row in self.payload["state_candidate_rows"]}
        self.assertTrue(rows["zero_higgs_state"]["conserved_without_extra_condition"])
        self.assertFalse(rows["zero_higgs_state"]["nonzero_stress_energy"])
        self.assertTrue(rows["constant_imported_vev_v0"]["nonzero_stress_energy"])
        self.assertFalse(rows["constant_imported_vev_v0"]["strict_state_or_vev_source_exported"])
        self.assertEqual(rows["rolling_log_profile_free_potential"]["conservation_residual"], "k**2*(3*p - 1)/t**3")

    def test_no_strict_source_and_docs_updated(self):
        self.assertFalse(any(row["strict_state_or_vev_source_exported"] and row["unit_dimension_source_exported"] for row in self.payload["state_candidate_rows"]))
        self.assertIn("P3155", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))
        self.assertIn("P3154/S2104", MD.read_text(encoding="utf-8"))
        self.assertIn("P3154/S2104", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3154/S2104", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3154", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
