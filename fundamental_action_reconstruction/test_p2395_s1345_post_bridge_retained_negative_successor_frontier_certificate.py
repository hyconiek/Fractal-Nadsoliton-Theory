from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2395_s1345_post_bridge_retained_negative_successor_frontier_certificate.py"
OUT = ROOT / "generated" / "p2395_s1345_post_bridge_retained_negative_successor_frontier_certificate.json"
MD = ROOT / "generated" / "p2395_s1345_post_bridge_retained_negative_successor_frontier_certificate.md"
PREREQ = ROOT / "p2394_s1344_apd_bridge_chi11_rebased_role_frontier_certificate.py"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2395PostBridgeRetainedNegativeSuccessorFrontierCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(PREREQ)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["post_bridge_retained_negative_successor_frontier_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2395_s1345_v1")
        self.assertEqual(self.payload["packet_id"], "P2395")
        self.assertEqual(self.payload["stage_id"], "S1345")
        self.assertEqual(self.payload["result_kind"], "POST_BRIDGE_RETAINED_NEGATIVE_SUCCESSOR_FRONTIER_CERTIFICATE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_context_inherits_apd_and_chi11_rebase(self) -> None:
        context = self.cert["closed_context"]
        self.assertTrue(context["apd_bridge_found"])
        self.assertTrue(context["chi11_selector_found"])
        self.assertEqual(context["p2394_current_assignment"]["closed_role_count"], 0)

    def test_retained_negative_closed_but_successor_open(self) -> None:
        rows = self.cert["role_rows"]
        self.assertEqual(len(rows), 3)
        self.assertTrue(
            all(row["retained_unchanged_transfer_verdict"] == "CLOSED_NEGATIVE_CURRENT_REPO_STATE" for row in rows)
        )
        self.assertTrue(all(row["modified_successor_verdict"] == "OPEN_ACTIVE_SUCCESSOR_FRONTIER" for row in rows))
        self.assertTrue(all(row["rejection_forever_verdict"] == "NOT_PROVEN_FOREVER_REJECTION" for row in rows))
        self.assertTrue(all(row["current_transfer_licensed"] is False for row in rows))

    def test_successor_matrix_and_minimal_atoms(self) -> None:
        frontier = self.cert["successor_frontier"]
        self.assertEqual(self.cert["matrix_rank_mod2"], 3)
        self.assertEqual(frontier["retained_negative_closed_count"], 3)
        self.assertEqual(frontier["current_licensed_role_count"], 0)
        self.assertEqual(frontier["open_modified_successor_count"], 3)
        self.assertEqual(
            frontier["successor_atom_universe"],
            [
                "alpha_geo_electroweak_role_theorem",
                "beta_power_hierarchy_successor_theorem",
                "beta_tors_strict_role_theorem",
            ],
        )
        self.assertEqual(frontier["successor_atom_degrees"]["beta_tors_strict_role_theorem"], 2)

    def test_fingerprint_gatekeepers_and_docs(self) -> None:
        self.assertEqual(self.cert["theorem_fingerprint_sha256"], sha256_json(self.theorem))
        self.assertIn("No unchanged legacy physical-role transfer is reopened.", self.theorem["not_licensed"])
        self.assertIn("No forever rejection of all strict successor roles is claimed.", self.theorem["not_licensed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("post-bridge retained-negative", MD.read_text(encoding="utf-8"))
        self.assertIn("P2395/S1345 post-bridge", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2395/S1345 retained-negative", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
