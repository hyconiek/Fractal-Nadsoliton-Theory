import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3122_s2072_iota_irrev_irreversible_defect_source_audit.py"
OUT = ROOT / "generated" / "p3122_s2072_iota_irrev_irreversible_defect_source_audit.json"
MD = ROOT / "generated" / "p3122_s2072_iota_irrev_irreversible_defect_source_audit.md"


class P3122IotaIrrevIrreversibleDefectSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3122_IOTA_IRREV_IRREVERSIBLE_DEFECT_SOURCE_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3121"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3121_accepted_Kappa_cycle_sources"], 0)
        self.assertEqual(cert["candidate_Iota_irrev_sources"], 15)
        self.assertEqual(cert["defect_law_rows"], 165)
        self.assertEqual(cert["signed_witness_rows"], 120)
        self.assertEqual(cert["Kappa_Tau_Xi_R_coupling_rows"], 135)
        self.assertEqual(cert["candidate_gate_rows"], 255)
        self.assertEqual(cert["accepted_Iota_irrev_sources"], 0)

    def test_candidate_matrix_and_rejections(self):
        objs = self.payload["constructed_theoretical_objects"]
        candidates = objs["candidate_Iota_irrev_sources"]
        self.assertTrue(any(row["candidate"] == "entropy_increment_defect" and row["signed_defect_value_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "z12_commutator_defect" and not row["gauge_invariant_not_label_choice"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "thermodynamic_entropy_production_defect" and not row["strict_nadsoliton_data_only"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "selector_signed_defect" and not row["selector_bridge_ltotal_toe_free"] for row in candidates))
        self.assertTrue(all(not row["accepted_Iota_irrev_source"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["accepted_coupling_chain"] for row in objs["Kappa_Tau_Xi_R_coupling_rows"]))

    def test_negative_flags_docs_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["candidate_Iota_irrev_sources_constructed"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("strict asymmetric-transition source object Delta_asym", decision["next_honest_step"])
        self.assertIn("forward/reverse asymmetry witness", decision["next_honest_step"])
        self.assertIn("P3122/S2072", MD.read_text(encoding="utf-8"))
        self.assertIn("P3122/S2072", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3122/S2072", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3122", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
