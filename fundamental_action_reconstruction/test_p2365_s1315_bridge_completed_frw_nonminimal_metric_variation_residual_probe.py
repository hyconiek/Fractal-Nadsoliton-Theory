from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2365_s1315_bridge_completed_frw_nonminimal_metric_variation_residual_probe.py"
OUT = ROOT / "generated" / "p2365_s1315_bridge_completed_frw_nonminimal_metric_variation_residual_probe.json"
MD = ROOT / "generated" / "p2365_s1315_bridge_completed_frw_nonminimal_metric_variation_residual_probe.md"

PREREQ_SCRIPTS = [
    ROOT / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.py",
    ROOT / "p1872_s822_strict_frw_full_lagrangian_stress_energy_export_probe.py",
    ROOT / "p1873_s823_strict_frw_full_component_stress_energy_table_probe.py",
    ROOT / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.py",
    ROOT / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.py",
    ROOT / "p2364_s1314_bridge_completed_frw_scalar_gauge_gravity_residual_table_probe.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2365BridgeCompletedFrwNonminimalMetricVariationProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["bridge_completed_frw_nonminimal_metric_variation_residual_probe"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2365_s1315_v1")
        self.assertEqual(self.payload["packet_id"], "P2365")
        self.assertEqual(self.payload["stage_id"], "S1315")
        self.assertEqual(
            self.payload["result_kind"],
            "BRIDGE_COMPLETED_FRW_NONMINIMAL_METRIC_VARIATION_RESIDUAL_LIFT",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_coefficient_source_and_nonminimal_operator(self) -> None:
        coeff = self.probe["coefficient_source"]
        self.assertEqual(coeff["source"], "P2363 bridge-completed APD moment transport")
        self.assertFalse(coeff["raw_legacy_moments_used"])
        self.assertEqual(
            set(coeff["effective_couplings"]),
            {"m2_eff", "lam_eff", "g_eff", "xi_eff"},
        )
        convention = self.probe["frw_nonminimal_variation_convention"]
        self.assertEqual(convention["term"], "xi_eff*Phi^2*R")
        self.assertIn("F*G_mn", convention["variation_operator"])
        self.assertIn("phidot", convention["Fdot"])
        self.assertIn("phiddot", convention["Fddot"])
        self.assertIn("2*kappa2", convention["normalization"])

    def test_metric_residual_solution_and_minimal_limit(self) -> None:
        solution = self.probe["nonminimal_metric_solution"]
        self.assertEqual(
            solution["solution_residuals"],
            {"R_g_00_nonminimal": "0", "R_g_ii_nonminimal": "0"},
        )
        limit = self.probe["minimal_limit_against_p2364"]["minimal_limit_residuals"]
        self.assertTrue(all(value == "0" for value in limit.values()))
        rows = self.probe["nonminimal_metric_residual_rows"]
        self.assertIn("Delta_00_nonminimal", rows)
        self.assertIn("Delta_ii_nonminimal", rows)
        self.assertIn("xi", rows["Delta_00_nonminimal"])
        self.assertIn("xi", rows["Delta_ii_nonminimal"])

    def test_gap_movement_and_limits(self) -> None:
        gap = self.probe["gap_movement"]
        self.assertEqual(
            gap["advanced_gap"],
            "full tensor-resolved nonminimal metric variation on named FRW family",
        )
        self.assertIn("Bianchi-I/anisotropic replay", gap["still_open"])
        theorem = self.probe["theorem_export"]
        self.assertIn("anisotropic Bianchi-I replay", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 discharge", theorem["not_licensed"])
        self.assertIn("legacy physical-role transfer", theorem["not_licensed"])

    def test_selector_separation_gatekeepers_and_fingerprint(self) -> None:
        separation = self.probe["separation_from_selector"]
        self.assertTrue(separation["eom_lagrangian_track_is_selector_independent"])
        self.assertTrue(separation["selector_qw2191_remains_parallel"])
        self.assertFalse(separation["selector_required_for_metric_variation_lift"])
        self.assertFalse(separation["legacy_role_transfer_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("xi_eff*Phi^2*R", theorem["claim"])
        self.assertIn("diagonal Bianchi-I family", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
