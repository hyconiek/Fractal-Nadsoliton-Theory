from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2347_s1297_strict_track_b_spherical_cap_chern_boundary_form_derivation.py"
OUT = ROOT / "generated" / "p2347_s1297_strict_track_b_spherical_cap_chern_boundary_form_derivation.json"
MD = ROOT / "generated" / "p2347_s1297_strict_track_b_spherical_cap_chern_boundary_form_derivation.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2347SphericalCapChernBoundaryFormTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_track_b_spherical_cap_chern_boundary_form_derivation_probe"]
        cls.witness = cls.probe["track_B_spherical_cap_chern_boundary_form_derivation"]
        cls.derivation = cls.witness["chern_form_derivation"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2347_s1297_v1")
        self.assertEqual(self.payload["packet_id"], "P2347")
        self.assertEqual(self.payload["stage_id"], "S1297")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_chern_boundary_form_derivation_values(self) -> None:
        self.assertEqual(self.witness["ledger_coefficient_b_GB"], "13152087*log(2)/(320000000*pi**2)")
        self.assertEqual(self.witness["per_euler_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.derivation["boundary_area"], "2*pi**2*(1 - c**2)**(3/2)")
        self.assertEqual(self.derivation["geodesic_principal_curvature"], "c/sqrt(1-c**2)")
        self.assertEqual(self.derivation["integrated_boundary_transgression"], "-16*pi**2*c*(c**2 - 3)")
        self.assertEqual(self.derivation["target_boundary_transgression"], "-16*pi**2*c*(c**2 - 3)")
        self.assertEqual(self.derivation["chern_boundary_residual"], "0")
        self.assertEqual(self.derivation["bulk_from_p2346_formula"], "16*pi**2*(c - 1)**2*(c + 2)")
        self.assertEqual(self.derivation["total_with_derived_boundary"], "32*pi**2")
        self.assertEqual(self.derivation["total_residual"], "0")
        self.assertEqual(self.derivation["total_derivative_residual"], "0")
        self.assertEqual(self.derivation["flat_limit_boundary_c1"], "32*pi**2")
        self.assertEqual(self.derivation["hemisphere_boundary_c0"], "0")
        self.assertEqual(self.derivation["sample_c_half_total"], "32*pi**2")
        self.assertEqual(self.derivation["normalized_track_B_pairing"], "13152087*log(2)/10000000")
        self.assertEqual(self.derivation["normalized_pairing_residual"], "0")

    def test_dependencies_gatekeepers_and_fingerprint(self) -> None:
        deps = self.probe["current_export_dependencies"]
        self.assertTrue(deps["p2346_boundary_transgression_derived"])
        self.assertTrue(deps["p2346_total_replayed"])
        for key, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, key)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("general Chern boundary form for arbitrary hypersurfaces", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])


if __name__ == "__main__":
    unittest.main()
