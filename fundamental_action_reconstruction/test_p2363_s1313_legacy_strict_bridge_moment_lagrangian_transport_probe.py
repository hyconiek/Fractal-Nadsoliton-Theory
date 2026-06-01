from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.py"
OUT = ROOT / "generated" / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.json"
MD = ROOT / "generated" / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.md"

PREREQ_SCRIPTS = [
    ROOT / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.py",
    ROOT / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2363LegacyStrictBridgeMomentLagrangianTransportProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["legacy_strict_bridge_moment_lagrangian_transport_probe"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2363_s1313_v1")
        self.assertEqual(self.payload["packet_id"], "P2363")
        self.assertEqual(self.payload["stage_id"], "S1313")
        self.assertEqual(
            self.payload["result_kind"],
            "LEGACY_STRICT_BRIDGE_COMPLETED_MOMENT_LAGRANGIAN_TRANSPORT_CERTIFICATE",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_symbolic_bridge_and_moment_transport(self) -> None:
        bridge = self.probe["symbolic_bridge"]
        self.assertEqual(bridge["completed_minus_strict_residual"], "0")
        moments = self.probe["moment_transport"]
        self.assertEqual(moments["completed_minus_strict_residuals"], {"c0": "0", "c1": "0", "c2": "0"})
        self.assertEqual(moments["p1866_strict_export_residuals"], {"c0": "0", "c1": "0", "c2": "0"})
        self.assertIn("alpha_geo*cos", bridge["legacy_kernel"])
        self.assertIn("Q_APD", bridge["completion_factor"])

    def test_effective_coupling_transport_matches_p1866(self) -> None:
        couplings = self.probe["effective_coupling_transport"]
        self.assertEqual(
            set(couplings["map"]),
            {"m2_eff", "lam_eff", "y_eff", "g_eff", "xi_eff"},
        )
        self.assertTrue(all(value == "0" for value in couplings["completed_minus_strict_residuals"].values()))
        self.assertTrue(all(value == "0" for value in couplings["p1866_strict_export_residuals"].values()))
        self.assertEqual(couplings["map"]["m2_eff"], "m2*(1 + c0)")
        self.assertEqual(couplings["map"]["lam_eff"], "lam*(1 + c1**2)")

    def test_negative_control_blocks_silent_legacy_substitution(self) -> None:
        numeric = self.probe["canonical_numeric_witness"]
        self.assertLess(numeric["completed_numeric_max_abs_residual"], 1e-12)
        self.assertGreater(numeric["legacy_raw_numeric_max_abs_mismatch"], 1e-2)
        self.assertGreater(numeric["legacy_amplitude_normalized_numeric_max_abs_mismatch"], 1e-2)
        self.assertNotEqual(numeric["legacy_raw_moments"]["c0"], numeric["strict_moments"]["c0"])
        self.assertNotEqual(numeric["legacy_amplitude_normalized_moments"]["c1"], numeric["strict_moments"]["c1"])

    def test_selector_separation_and_gatekeepers(self) -> None:
        separation = self.probe["separation_from_selector"]
        self.assertTrue(separation["eom_lagrangian_track_is_selector_independent"])
        self.assertTrue(separation["selector_closure_is_parallel_problem"])
        self.assertFalse(separation["selector_required_for_moment_transport"])
        self.assertEqual(
            separation["selector_qw2191_status"],
            "OPEN_PARALLEL_TRACK_NOT_USED_AS_EOM_INPUT",
        )
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_theorem_fingerprint_and_limits(self) -> None:
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("Raw legacy moments", theorem["claim"])
        self.assertIn("selector premise or QW-2191 discharge", theorem["not_licensed"])
        self.assertIn("legacy physical-role transfer", theorem["not_licensed"])
        self.assertIn("full 4D EOM theorem-grade closure", theorem["not_licensed"])
        self.assertIn("bridge-completed moments", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
