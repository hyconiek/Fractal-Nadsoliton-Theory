from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2680_s1630_legacy_strict_bridge_source_inventory_no_selector_replay_audit.py"
OUT = ROOT / "generated" / "p2680_s1630_legacy_strict_bridge_source_inventory_no_selector_replay_audit.json"
MD = ROOT / "generated" / "p2680_s1630_legacy_strict_bridge_source_inventory_no_selector_replay_audit.md"


class P2680LegacyStrictBridgeSourceInventoryNoSelectorReplayAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_covers_nonselector_bridge_sources(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "legacy_strict_kernel_content",
            "amplitude_normalization_content",
            "damping_compression_content",
            "positive_beta_source_content",
            "canonical_unit_content",
            "forbidden_replay_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_upstream_p2679_pivot_respected(self) -> None:
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2679_blocks_old_lanes"])
        self.assertTrue(upstream["p2679_no_new_reopening_evidence"])
        self.assertTrue(upstream["p2679_allowed_bridge_source_pivot"])

    def test_source_atom_inventory_separates_formal_material_from_sources(self) -> None:
        atoms = {row["atom"]: row for row in self.payload["source_atom_inventory"]}
        self.assertTrue(atoms["alpha_geo_scalar_shape_normalization"]["formal_material_present"])
        self.assertFalse(atoms["alpha_geo_scalar_shape_normalization"]["role_safe_for_completion"])
        self.assertTrue(atoms["fractal_pushforward_linear_to_power_damping_shape"]["formal_material_present"])
        self.assertFalse(atoms["target_independent_positive_beta_or_z_beta_source"]["source_theorem_exported"])
        self.assertFalse(atoms["canonical_length_or_uv_unit_source"]["source_theorem_exported"])
        self.assertFalse(atoms["selector_phase_orientation_source"]["source_theorem_exported"])

    def test_bridge_completion_lattice(self) -> None:
        lattice = self.payload["bridge_completion_lattice"]
        self.assertEqual(lattice["total_states"], 256)
        self.assertEqual(lattice["passing_states"], 1)
        self.assertEqual(lattice["hamming_distance_to_pass"], 3)
        self.assertFalse(lattice["full_bridge_completed_now"])
        self.assertIn("amplitude_role_safe_source_exported", lattice["missing_current_obligations"])
        self.assertIn("positive_beta_z_beta_source_exported", lattice["missing_current_obligations"])
        self.assertIn("canonical_length_uv_unit_exported", lattice["missing_current_obligations"])

    def test_no_selector_replay_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertFalse(decision["full_bridge_completed_now"])
        self.assertFalse(decision["selector_replay_used"])
        self.assertFalse(decision["beta_tors_chi11_reopened_now"])
        self.assertFalse(decision["role_transfer_allowed_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2680/S1630", MD.read_text(encoding="utf-8"))
        self.assertIn("P2680/S1630", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2680/S1630", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
