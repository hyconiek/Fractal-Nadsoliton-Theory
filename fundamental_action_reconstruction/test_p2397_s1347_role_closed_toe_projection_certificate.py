from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2397_s1347_role_closed_toe_projection_certificate.py"
OUT = ROOT / "generated" / "p2397_s1347_role_closed_toe_projection_certificate.json"
MD = ROOT / "generated" / "p2397_s1347_role_closed_toe_projection_certificate.md"
PREREQ = ROOT / "p2396_s1346_role_package_negative_closure_rebase_certificate.py"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2397RoleClosedToeProjectionCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(PREREQ)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["role_closed_toe_projection_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2397_s1347_v1")
        self.assertEqual(self.payload["packet_id"], "P2397")
        self.assertEqual(self.payload["stage_id"], "S1347")
        self.assertEqual(self.payload["result_kind"], "ROLE_CLOSED_TOE_PROJECTION_CERTIFICATE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_role_closed_slice_counts(self) -> None:
        slice_cert = self.cert["role_closed_slice"]
        self.assertEqual(slice_cert["slice_row_count"], 16)
        self.assertEqual(slice_cert["role_transfer_true_count"], 0)
        self.assertEqual(slice_cert["toe_true_count"], 0)
        self.assertEqual(slice_cert["bridge_true_count"], 2)
        self.assertEqual(slice_cert["selector_true_count"], 8)
        self.assertEqual(slice_cert["signature_rank_mod2"], 2)

    def test_role_atoms_forced_false_and_nonrole_progress_separated(self) -> None:
        slice_cert = self.cert["role_closed_slice"]
        self.assertEqual(
            slice_cert["role_atoms_forced_false"],
            [
                "alpha_geo_electroweak_role_theorem",
                "beta_tors_strict_role_theorem",
                "beta_power_hierarchy_successor_theorem",
            ],
        )
        signatures = slice_cert["signature_counts"]
        self.assertEqual(sum(signatures.values()), 16)
        self.assertIn("0000", signatures)
        self.assertIn("1010", signatures)
        self.assertTrue(all(row["target_values"]["toe_closure"] is False for row in slice_cert["rows"]))

    def test_fingerprint_gatekeepers_docs_and_limits(self) -> None:
        self.assertEqual(self.cert["theorem_fingerprint_sha256"], sha256_json(self.theorem))
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("No role-transfer theorem is recovered by non-role atoms alone.", self.theorem["not_licensed"])
        self.assertIn("No forever impossibility of future explicit role-successor evidence is claimed.", self.theorem["not_licensed"])
        self.assertIn("role-closed ToE projection", MD.read_text(encoding="utf-8"))
        self.assertIn("P2397/S1347 role-closed", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2397/S1347 role-closed", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
