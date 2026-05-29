from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2343_s1293_strict_track_b_flat_round_d4_boundary_functional_derivation.py"
OUT = ROOT / "generated" / "p2343_s1293_strict_track_b_flat_round_d4_boundary_functional_derivation.json"
MD = ROOT / "generated" / "p2343_s1293_strict_track_b_flat_round_d4_boundary_functional_derivation.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2343TrackBFlatRoundD4BoundaryFunctionalTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_flat_round_d4_boundary_functional_derivation_probe"]
        cls.derivation = cls.probe["track_B_flat_round_D4_boundary_functional_derivation"]
        cls.flat_round = cls.derivation["flat_round_boundary_class"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2343_s1293_v1")
        self.assertEqual(self.payload["packet_id"], "P2343")
        self.assertEqual(self.payload["stage_id"], "S1293")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_flat_round_boundary_functional_derivation_values(self) -> None:
        self.assertEqual(self.derivation["ledger_coefficient_b_GB"], "13152087*log(2)/(320000000*pi**2)")
        self.assertEqual(self.derivation["per_euler_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.flat_round["bulk_gb_fixture_value"], "0")
        self.assertEqual(self.flat_round["boundary_area"], "2*pi**2*R**3")
        self.assertEqual(self.flat_round["shape_operator_principal_curvatures"], ["1/R", "1/R", "1/R"])
        self.assertEqual(self.flat_round["sigma3_shape_operator"], "R**(-3)")
        self.assertEqual(self.flat_round["boundary_density_normalization_factor"], "16")
        self.assertEqual(self.flat_round["boundary_functional_B_flat_round_R"], "32*pi**2")
        self.assertEqual(self.flat_round["target_topological_number_32pi2_chi"], "32*pi**2")
        self.assertEqual(self.flat_round["boundary_functional_residual"], "0")
        self.assertEqual(self.flat_round["radius_derivative_residual"], "0")
        self.assertEqual(self.flat_round["unit_radius_value"], "32*pi**2")
        self.assertEqual(self.flat_round["normalized_track_B_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.flat_round["normalized_pairing_residual"], "0")
        self.assertEqual(self.flat_round["two_copy_boundary_functional"], "64*pi**2")
        self.assertEqual(self.flat_round["two_copy_residual_against_p2342_S4"], "0")

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        deps = self.probe["current_export_dependencies"]
        self.assertTrue(deps["p2341_boundary_value_rederived_at_unit_radius"])
        self.assertTrue(deps["p2341_pairing_rederived_at_unit_radius"])
        self.assertTrue(deps["p2342_two_copy_s4_number_rederived"])
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("boundary functional for non-round or curved-boundary four-manifolds", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
