from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2403_s1353_strict_kernel_primary_physics_generation_rebase_certificate.py"
OUT = ROOT / "generated" / "p2403_s1353_strict_kernel_primary_physics_generation_rebase_certificate.json"
MD = ROOT / "generated" / "p2403_s1353_strict_kernel_primary_physics_generation_rebase_certificate.md"
PREREQ = ROOT / "p2402_s1352_role_successor_marginal_credit_certificate.py"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2403StrictKernelPrimaryPhysicsGenerationRebaseCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(PREREQ)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_kernel_primary_physics_generation_rebase_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2403_s1353_v1")
        self.assertEqual(self.payload["packet_id"], "P2403")
        self.assertEqual(self.payload["stage_id"], "S1353")
        self.assertEqual(self.payload["result_kind"], "STRICT_KERNEL_PRIMARY_PHYSICS_GENERATION_REBASE_CERTIFICATE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_characteristic_matrix_rebases_strict_as_primary(self) -> None:
        summary = self.theorem["characteristic_matrix_summary"]
        self.assertEqual(summary["legacy_characteristic_count"], 1)
        self.assertEqual(summary["strict_characteristic_count"], 5)
        self.assertEqual(summary["strict_addition_count"], 4)
        self.assertEqual(summary["strict_structural_dominance_delta"], 4)
        self.assertIn("nonlinear_damping_compression", summary["strict_additions"])

    def test_lane_rebase_and_inherited_bridge_state(self) -> None:
        lane_summary = self.theorem["lane_rebase_summary"]
        self.assertTrue(lane_summary["all_lanes_rebased_to_strict_primary"])
        self.assertTrue(lane_summary["all_lanes_keep_legacy_as_bridge_not_final_kernel"])
        self.assertTrue(lane_summary["all_lanes_block_silent_role_transfer"])
        self.assertTrue(self.theorem["p2394_apd_bridge_found"])
        self.assertTrue(self.theorem["p2394_chi11_selector_found"])

    def test_fingerprint_gatekeepers_docs_and_limits(self) -> None:
        self.assertEqual(self.cert["theorem_fingerprint_sha256"], sha256_json(self.theorem))
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("No theorem that strict already generates SM/GR physical roles is exported here.", self.theorem["not_licensed"])
        self.assertIn("strict-kernel primary physics-generation", MD.read_text(encoding="utf-8"))
        self.assertIn("P2403/S1353 strict-primary", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2403/S1353 strict-primary", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
