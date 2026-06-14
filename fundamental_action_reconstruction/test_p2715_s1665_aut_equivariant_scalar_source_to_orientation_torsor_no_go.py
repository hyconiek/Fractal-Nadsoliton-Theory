from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2715_s1665_aut_equivariant_scalar_source_to_orientation_torsor_no_go.py"
OUT = ROOT / "generated" / "p2715_s1665_aut_equivariant_scalar_source_to_orientation_torsor_no_go.json"
MD = ROOT / "generated" / "p2715_s1665_aut_equivariant_scalar_source_to_orientation_torsor_no_go.md"


class P2715AutEquivariantScalarSourceToOrientationTorsorNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_no_equivariant_maps_from_trivial_scalar_domains(self):
        self.assertEqual(self.payload["status"], "P2715_AUT_TRIVIAL_SCALAR_SOURCE_TO_ORIENTATION_TORSOR_NO_GO")
        summary = self.payload["finite_equivariance_summary"]
        self.assertEqual(summary["orientation_reversing_units"], [7, 11])
        self.assertEqual(summary["one_point_candidate_maps"], 2)
        self.assertEqual(summary["one_point_equivariant_maps"], 0)
        self.assertEqual(summary["two_point_trivial_branch_candidate_maps"], 4)
        self.assertEqual(summary["two_point_trivial_branch_equivariant_maps"], 0)

    def test_scalar_source_matrix_is_orientation_blind(self):
        classes = {row["source_class"] for row in self.payload["scalar_source_matrix"]}
        self.assertEqual(classes, {
            "entropy_or_uv_scale_scalar",
            "alpha_geo_amplitude_normalization_scalar",
            "positive_beta_or_z_beta_damping_scalar",
            "scalar_lagrangian_density_or_metric_residual",
        })
        self.assertTrue(all(row["aut_action"] == "trivial_scalar" for row in self.payload["scalar_source_matrix"]))
        self.assertTrue(all(not row["breaks_orientation_torsor"] for row in self.payload["scalar_source_matrix"]))
        decision = self.payload["decision"]
        self.assertFalse(decision["scalar_source_breaks_orientation_torsor"])
        self.assertFalse(decision["strict_pseudoscalar_or_chiral_source_exported"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_written(self):
        self.assertIn("P2715/S1665", MD.read_text(encoding="utf-8"))
        self.assertIn("P2715/S1665", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2715/S1665", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2715", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
