import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3114_s2064_dimensional_triad_source_law_audit.py"
OUT = ROOT / "generated" / "p3114_s2064_dimensional_triad_source_law_audit.json"
MD = ROOT / "generated" / "p3114_s2064_dimensional_triad_source_law_audit.md"


class P3114DimensionalTriadSourceLawAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3114_DIMENSIONAL_TRIAD_SOURCE_LAW_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3113"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3113_accepted_U_action_source_laws"], 0)
        self.assertEqual(cert["candidate_D_phi_source_laws"], 7)
        self.assertEqual(cert["scale_orbit_section_rows"], 35)
        self.assertEqual(cert["carrier_axis_rows"], 21)
        self.assertEqual(cert["relation_coupling_rows"], 28)
        self.assertEqual(cert["candidate_gate_rows"], 77)
        self.assertEqual(cert["accepted_D_phi_source_laws"], 0)

    def test_candidate_matrix_and_rejections(self):
        objs = self.payload["constructed_theoretical_objects"]
        candidates = objs["candidate_D_phi_source_laws"]
        self.assertTrue(any(row["candidate"] == "declared_dimension_symbol_triad" and row["nonzero_action_carrier_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "imported_planck_light_triad" and not row["standard_physics_import_free"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "internal_dispersion_velocity_triad" and row["action_length_time_relation_exported"] for row in candidates))
        self.assertTrue(all(not row["accepted_D_phi_source_law"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["relation_accepted"] for row in objs["relation_coupling_rows"]))

    def test_negative_flags_docs_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["candidate_D_phi_source_laws_constructed"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("scale-section theorem", decision["next_honest_step"])
        self.assertIn("Sigma_dim", decision["next_honest_step"])
        self.assertIn("P3114/S2064", MD.read_text(encoding="utf-8"))
        self.assertIn("P3114/S2064", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3114/S2064", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3114", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
