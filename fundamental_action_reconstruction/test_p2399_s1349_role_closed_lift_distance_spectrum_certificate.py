from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2399_s1349_role_closed_lift_distance_spectrum_certificate.py"
OUT = ROOT / "generated" / "p2399_s1349_role_closed_lift_distance_spectrum_certificate.json"
MD = ROOT / "generated" / "p2399_s1349_role_closed_lift_distance_spectrum_certificate.md"
PREREQ = ROOT / "p2398_s1348_role_closed_quotient_anf_certificate.py"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2399RoleClosedLiftDistanceSpectrumCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(PREREQ)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["role_closed_lift_distance_spectrum_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2399_s1349_v1")
        self.assertEqual(self.payload["packet_id"], "P2399")
        self.assertEqual(self.payload["stage_id"], "S1349")
        self.assertEqual(self.payload["result_kind"], "ROLE_CLOSED_LIFT_DISTANCE_SPECTRUM_CERTIFICATE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_lift_distance_spectra(self) -> None:
        self.assertEqual(self.theorem["minimum_role_transfer_lift_distance"], 3)
        self.assertEqual(self.theorem["minimum_toe_lift_distance"], 3)
        self.assertEqual(self.theorem["role_transfer_distance_spectrum"], {"3": 8, "4": 8})
        self.assertEqual(self.theorem["toe_distance_spectrum"], {"3": 1, "4": 4, "5": 6, "6": 4, "7": 1})
        self.assertEqual(self.theorem["nearest_toe_rows"], [15])

    def test_missing_matrix_ranks_and_nearest_row(self) -> None:
        spectrum = self.cert["lift_distance_spectrum"]
        self.assertEqual(spectrum["missing_matrix_ranks_mod2"]["role_transfer_theorem_level_closure"], 2)
        self.assertEqual(spectrum["missing_matrix_ranks_mod2"]["toe_closure"], 5)
        nearest = spectrum["lift_rows"][15]
        self.assertTrue(all(nearest["nonrole_assignment"].values()))
        self.assertEqual(nearest["toe_distance_decomposition"], {
            "forced_role_atoms_missing": 3,
            "nonrole_atoms_missing": 0,
            "total_distance": 3,
        })

    def test_fingerprint_gatekeepers_docs_and_limits(self) -> None:
        self.assertEqual(self.cert["theorem_fingerprint_sha256"], sha256_json(self.theorem))
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertTrue(self.theorem["p2398_role_transfer_zero_polynomial"])
        self.assertTrue(self.theorem["p2398_toe_zero_polynomial"])
        self.assertIn("No role-transfer theorem is created by a finite lift distance.", self.theorem["not_licensed"])
        self.assertIn("role-closed lift-distance spectrum", MD.read_text(encoding="utf-8"))
        self.assertIn("P2399/S1349 role-closed", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2399/S1349 lift-distance", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
