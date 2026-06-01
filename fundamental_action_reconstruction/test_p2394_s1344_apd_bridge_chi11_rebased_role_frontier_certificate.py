from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2394_s1344_apd_bridge_chi11_rebased_role_frontier_certificate.py"
OUT = ROOT / "generated" / "p2394_s1344_apd_bridge_chi11_rebased_role_frontier_certificate.json"
MD = ROOT / "generated" / "p2394_s1344_apd_bridge_chi11_rebased_role_frontier_certificate.md"
PREREQ = ROOT / "p2393_s1343_kernel_completion_boundary_residual_certificate.py"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2394APDBridgeChi11RebasedRoleFrontierCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(PREREQ)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["apd_bridge_chi11_rebased_role_frontier_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2394_s1344_v1")
        self.assertEqual(self.payload["packet_id"], "P2394")
        self.assertEqual(self.payload["stage_id"], "S1344")
        self.assertEqual(self.payload["result_kind"], "APD_BRIDGE_CHI11_REBASED_ROLE_FRONTIER_CERTIFICATE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_apd_bridge_and_chi11_are_rebased_as_found(self) -> None:
        apd = self.cert["apd_bridge_status"]
        selector = self.cert["selector_status"]
        self.assertTrue(apd["apd_bridge_found_in_existing_repo"])
        self.assertTrue(selector["strict_selector_found_in_declared_scope"])
        self.assertTrue(selector["selector_realized_without_beta_tors_chi11"])
        self.assertEqual(self.theorem["active_apd_reproof_obligation_count"], 0)
        self.assertEqual(self.theorem["active_missing_selector_route_count"], 0)

    def test_p2393_is_reinterpreted_as_boundary_negative_control(self) -> None:
        correction = self.cert["p2393_boundary_correction"]
        self.assertTrue(correction["p2393_boundary_identity_kept_as_negative_control"])
        self.assertTrue(correction["p2393_current_residual_nonzero"])
        self.assertTrue(correction["p2393_residual_reinterpreted_as_eta1_slice_insufficiency_not_missing_apd_bridge"])
        self.assertIn("not be read as evidence that the APD bridge was unfound", correction["correction"])

    def test_role_frontier_truth_table_after_rebase(self) -> None:
        frontier = self.cert["rebased_role_frontier"]
        self.assertEqual(frontier["truth_table_size"], 8)
        self.assertTrue(frontier["chi11_not_in_active_role_atom_set"])
        self.assertEqual(frontier["current_assignment"]["closed_role_count"], 0)
        self.assertEqual(
            frontier["minimal_supports"]["all_physical_roles_minimal_support"],
            [
                "alpha_geo_electroweak_role_theorem",
                "beta_power_hierarchy_successor_theorem",
                "beta_tors_strict_role_theorem",
            ],
        )
        singleton = {row["atom"]: row["closed_role_count"] for row in frontier["singleton_unlocks"]}
        self.assertEqual(singleton["alpha_geo_electroweak_role_theorem"], 1)
        self.assertEqual(singleton["beta_tors_strict_role_theorem"], 0)
        self.assertEqual(singleton["beta_power_hierarchy_successor_theorem"], 0)

    def test_fingerprint_gatekeepers_and_docs(self) -> None:
        self.assertEqual(self.cert["theorem_fingerprint_sha256"], sha256_json(self.theorem))
        self.assertIn("No beta_tors -> chi11 selector-search route is reopened.", self.theorem["not_licensed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("APD bridge found", MD.read_text(encoding="utf-8"))
        self.assertIn("P2394/S1344 APD bridge found", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2394/S1344 APD bridge rebase", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
