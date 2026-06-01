from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2400_s1350_nearest_lift_role_successor_lattice_certificate.py"
OUT = ROOT / "generated" / "p2400_s1350_nearest_lift_role_successor_lattice_certificate.json"
MD = ROOT / "generated" / "p2400_s1350_nearest_lift_role_successor_lattice_certificate.md"
PREREQ = ROOT / "p2399_s1349_role_closed_lift_distance_spectrum_certificate.py"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2400NearestLiftRoleSuccessorLatticeCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(PREREQ)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["nearest_lift_role_successor_lattice_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2400_s1350_v1")
        self.assertEqual(self.payload["packet_id"], "P2400")
        self.assertEqual(self.payload["stage_id"], "S1350")
        self.assertEqual(self.payload["result_kind"], "NEAREST_LIFT_ROLE_SUCCESSOR_LATTICE_CERTIFICATE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_only_full_role_subset_closes(self) -> None:
        self.assertEqual(self.theorem["role_transfer_true_masks"], [7])
        self.assertEqual(self.theorem["toe_true_masks"], [7])
        self.assertEqual(self.theorem["proper_subset_fail_count"], 7)
        self.assertEqual(self.theorem["nearest_one_role_missing_masks"], [3, 5, 6])

    def test_role_anf_and_derivatives(self) -> None:
        self.assertEqual(self.theorem["role_transfer_role_anf"], self.theorem["toe_role_anf"])
        self.assertEqual(len(self.theorem["toe_role_anf"]), 1)
        self.assertEqual(self.theorem["toe_role_anf"][0]["degree"], 3)
        self.assertTrue(all(len(support) == 1 for support in self.theorem["toe_boolean_derivative_support_by_role_atom"].values()))

    def test_lattice_context_and_limits(self) -> None:
        lattice = self.cert["role_successor_lattice"]
        self.assertEqual(len(lattice["rows"]), 8)
        self.assertEqual(lattice["target_truth_vectors_by_mask_0_to_7"]["toe_closure"], [0, 0, 0, 0, 0, 0, 0, 1])
        self.assertEqual(self.theorem["p2399_minimum_toe_lift_distance"], 3)
        self.assertIn("No one-role or two-role subset licenses role transfer.", self.theorem["not_licensed"])

    def test_fingerprint_gatekeepers_and_docs(self) -> None:
        self.assertEqual(self.cert["theorem_fingerprint_sha256"], sha256_json(self.theorem))
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("nearest-lift role-successor lattice", MD.read_text(encoding="utf-8"))
        self.assertIn("P2400/S1350 nearest-lift", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2400/S1350 role-successor", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
