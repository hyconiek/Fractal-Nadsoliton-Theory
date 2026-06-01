from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2364_s1314_bridge_completed_frw_scalar_gauge_gravity_residual_table_probe.py"
OUT = ROOT / "generated" / "p2364_s1314_bridge_completed_frw_scalar_gauge_gravity_residual_table_probe.json"
MD = ROOT / "generated" / "p2364_s1314_bridge_completed_frw_scalar_gauge_gravity_residual_table_probe.md"

PREREQ_SCRIPTS = [
    ROOT / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.py",
    ROOT / "p1867_s817_strict_covariant_eom_residual_witness_scaffold.py",
    ROOT / "p1870_s820_strict_frw_background_metric_residual_probe.py",
    ROOT / "p1872_s822_strict_frw_full_lagrangian_stress_energy_export_probe.py",
    ROOT / "p1873_s823_strict_frw_multisector_stress_energy_scaffold.py",
    ROOT / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.py",
    ROOT / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2364BridgeCompletedFrwResidualTableProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["bridge_completed_frw_scalar_gauge_gravity_residual_table_probe"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2364_s1314_v1")
        self.assertEqual(self.payload["packet_id"], "P2364")
        self.assertEqual(self.payload["stage_id"], "S1314")
        self.assertEqual(
            self.payload["result_kind"],
            "BRIDGE_COMPLETED_FRW_SCALAR_GAUGE_GRAVITY_RESIDUAL_TABLE",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_bridge_completed_coefficient_source(self) -> None:
        coeff = self.probe["coefficient_source"]
        self.assertEqual(coeff["source"], "P2363 bridge-completed APD moment transport")
        self.assertFalse(coeff["raw_legacy_moments_used"])
        self.assertEqual(set(coeff["moments"]), {"c0", "c1", "c2"})
        self.assertEqual(
            set(coeff["effective_couplings"]),
            {"m2_eff", "lam_eff", "g_eff", "xi_eff"},
        )
        ratios = coeff["canonical_numeric_coupling_ratios"]
        self.assertGreater(ratios["m2_eff_over_m2"], 1.0)
        self.assertGreater(ratios["lam_eff_over_lam"], 1.0)
        self.assertGreater(ratios["g_eff_over_g"], 1.0)

    def test_residual_table_rows_and_zero_residuals(self) -> None:
        rows = self.probe["residual_table"]
        self.assertEqual(
            {row["row_id"] for row in rows},
            {"E_phi_FRW", "E_A_FRW", "R_g_00_FRW", "R_g_ii_FRW"},
        )
        self.assertTrue(all(not row["selector_required"] for row in rows))
        self.assertEqual(
            self.probe["normal_form_residuals"],
            {
                "E_phi_after_phiddot_substitution": "0",
                "E_A_after_Addot_substitution": "0",
            },
        )
        self.assertEqual(
            self.probe["metric_solution"]["solution_residuals"],
            {"R_g_00": "0", "R_g_ii": "0"},
        )
        self.assertTrue(
            all(value == "0" for value in self.probe["residual_vector_after_closure_substitution"].values())
        )

    def test_gap_movement_preserves_open_limits(self) -> None:
        gap = self.probe["gap_movement"]
        self.assertEqual(
            gap["advanced_gap"],
            "nonproxy_covariant_sector_export_with_background_residual_tables",
        )
        self.assertIn("still_open_reason", gap)
        theorem = self.probe["theorem_export"]
        self.assertIn(
            "full nonproxy tensor-resolved metric variation with nonminimal stress-energy terms",
            theorem["not_licensed"],
        )
        self.assertIn("background-independence atlas lift beyond this FRW family", theorem["not_licensed"])

    def test_selector_separation_gatekeepers_and_fingerprint(self) -> None:
        separation = self.probe["separation_from_selector"]
        self.assertTrue(separation["eom_lagrangian_track_is_selector_independent"])
        self.assertTrue(separation["selector_qw2191_remains_parallel"])
        self.assertFalse(separation["selector_required_for_residual_table"])
        self.assertFalse(separation["legacy_role_transfer_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("P2363 bridge-completed moments", theorem["claim"])
        self.assertIn("selector premise or QW-2191 discharge", theorem["not_licensed"])
        self.assertIn("full tensor-resolved variation", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
