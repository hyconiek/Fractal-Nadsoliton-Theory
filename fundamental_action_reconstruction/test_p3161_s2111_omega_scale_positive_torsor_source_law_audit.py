import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3161_s2111_omega_scale_positive_torsor_source_law_audit.py"
OUT = ROOT / "generated" / "p3161_s2111_omega_scale_positive_torsor_source_law_audit.json"
MD = ROOT / "generated" / "p3161_s2111_omega_scale_positive_torsor_source_law_audit.md"


class P3161OmegaScalePositiveTorsorSourceLawAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_equivariant_obstruction_and_counts(self):
        self.assertEqual(self.payload["status"], "P3161_OMEGA_SCALE_POSITIVE_TORSOR_SOURCE_LAW_EQUIVARIANT_OBSTRUCTION")
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["candidate_Omega_scale_sources"], 12)
        self.assertEqual(cert["scale_action_factors"], 7)
        self.assertEqual(cert["equivariant_obstruction_rows"], 84)
        self.assertEqual(cert["candidate_gate_rows"], 132)
        self.assertEqual(cert["accepted_Omega_scale_sources"], 0)
        self.assertGreater(cert["obstructed_invariant_rows"], 0)

    def test_candidate_matrix_and_no_exports(self):
        objs = self.payload["constructed_theoretical_objects"]
        candidates = objs["candidate_Omega_scale_sources"]
        self.assertTrue(any(row["candidate"] == "hypothetical_scale_charged_source" and row["equivariant_section_possible"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "dimensionless_alpha_phi_section" and not row["scale_charged_input_exported"] for row in candidates))
        self.assertTrue(all(not row["accepted_Omega_scale_source"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["next_source_shape_S_plus_identified"])

    def test_docs_updated(self):
        self.assertIn("P3161/S2111", MD.read_text(encoding="utf-8"))
        self.assertIn("P3161/S2111", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3161/S2111", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3161", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
