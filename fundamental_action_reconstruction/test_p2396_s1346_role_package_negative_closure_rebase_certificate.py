from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2396_s1346_role_package_negative_closure_rebase_certificate.py"
OUT = ROOT / "generated" / "p2396_s1346_role_package_negative_closure_rebase_certificate.json"
MD = ROOT / "generated" / "p2396_s1346_role_package_negative_closure_rebase_certificate.md"
PREREQ = ROOT / "p2395_s1345_post_bridge_retained_negative_successor_frontier_certificate.py"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2396RolePackageNegativeClosureRebaseCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(PREREQ)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["role_package_negative_closure_rebase_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2396_s1346_v1")
        self.assertEqual(self.payload["packet_id"], "P2396")
        self.assertEqual(self.payload["stage_id"], "S1346")
        self.assertEqual(self.payload["result_kind"], "ROLE_PACKAGE_NEGATIVE_CLOSURE_REBASE_CERTIFICATE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_all_claim_specific_closures_rebased(self) -> None:
        rows = self.cert["role_rows"]
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(row["retained_branch_closed_negative"] for row in rows))
        self.assertTrue(all(row["replaced_branch_closed_negative"] for row in rows))
        self.assertTrue(all(row["full_claim_specific_closed_negative"] for row in rows))
        self.assertTrue(all(row["current_state_role_transfer_verdict"] == "CLOSED_NEGATIVE_CURRENT_REPO_STATE" for row in rows))

    def test_p2395_flags_are_future_only_not_current_transfer(self) -> None:
        rows = self.cert["role_rows"]
        self.assertTrue(all(row["p2395_modified_successor_flag"] for row in rows))
        self.assertTrue(
            all(row["p2395_modified_successor_rebased_verdict"] == "FUTURE_ONLY_CONDITIONAL_NOT_CURRENT_EXPORT" for row in rows)
        )
        self.assertTrue(all(row["future_successor_not_forbidden_forever"] for row in rows))
        self.assertTrue(all(row["current_transfer_licensed"] is False for row in rows))

    def test_package_matrix_and_gatekeepers(self) -> None:
        package = self.cert["package_certificate"]
        self.assertTrue(package["n116_package_closed_negative"])
        self.assertTrue(package["all_current_role_transfer_closed_negative"])
        self.assertEqual(package["current_licensed_transfer_count"], 0)
        self.assertEqual(package["p2395_future_only_flag_count"], 3)
        self.assertEqual(package["future_successor_not_forbidden_count"], 3)
        self.assertEqual(package["matrix_rank_mod2"], 1)
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_fingerprint_docs_and_limits(self) -> None:
        self.assertEqual(self.cert["theorem_fingerprint_sha256"], sha256_json(self.theorem))
        self.assertIn("No P2395 modified-successor flag is a current exported successor theorem.", self.theorem["not_licensed"])
        self.assertIn("No forever impossibility of future strict successor evidence is claimed.", self.theorem["not_licensed"])
        self.assertIn("role-package negative closure", MD.read_text(encoding="utf-8"))
        self.assertIn("P2396/S1346 role-package", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2396/S1346 role-package", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
