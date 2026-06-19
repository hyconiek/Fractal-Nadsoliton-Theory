import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2951_s1901_ratio_package_strict_source_normal_form_lattice.py"
OUT = ROOT / "generated" / "p2951_s1901_ratio_package_strict_source_normal_form_lattice.json"
MD = ROOT / "generated" / "p2951_s1901_ratio_package_strict_source_normal_form_lattice.md"


class P2951RatioPackageStrictSourceNormalFormLatticeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2951_RATIO_PACKAGE_STRICT_SOURCE_NORMAL_FORM_LATTICE_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2948"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2949"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2950"])

    def test_normal_form_certificate(self):
        cert = self.payload["normal_form_certificate"]
        self.assertEqual(cert["obligation_count"], 4)
        self.assertEqual(cert["truth_table_row_count"], 16)
        self.assertEqual(cert["accepted_row_count"], 1)
        self.assertTrue(cert["unique_accepting_row_requires_all_atoms"])
        self.assertTrue(cert["all_atoms_essential"])
        self.assertEqual(cert["current_present_count"], 0)
        self.assertEqual(len(cert["current_missing_atoms"]), 4)
        self.assertFalse(cert["current_artifact_accepts_strict_ratio_damping_source"])

    def test_lattice_rows_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        rows = obj["truth_table_rows"]
        accepted = [row for row in rows if row["accepted_strict_ratio_damping_source"]]
        self.assertEqual(len(accepted), 1)
        self.assertEqual(accepted[0]["present_count"], 4)
        self.assertTrue(all(row["essential"] for row in obj["essentiality_rows"]))
        self.assertEqual(obj["current_artifact_row"]["present_count"], 0)
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_ratio_package_source_theorem_exported", "strict_torsion_character_source_theorem_exported", "strict_delta_numerator_semantics_exported", "strict_positive_beta_scale_source_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2951/S1901", MD.read_text(encoding="utf-8"))
        self.assertIn("P2951/S1901", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2951/S1901", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2951", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
