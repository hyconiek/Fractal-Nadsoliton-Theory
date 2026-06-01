from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2391_s1341_selector_epoch_rebased_bridge_gap_matrix.py"
OUT = ROOT / "generated" / "p2391_s1341_selector_epoch_rebased_bridge_gap_matrix.json"
MD = ROOT / "generated" / "p2391_s1341_selector_epoch_rebased_bridge_gap_matrix.md"
PREREQ = ROOT / "p2390_s1340_selector_qualified_beta_tors_chi11_role_audit.py"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2391SelectorEpochRebasedBridgeGapMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(PREREQ)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["selector_epoch_rebased_bridge_gap_matrix"]
        cls.rows = {row["component"]: row for row in cls.probe["rows"]}

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2391_s1341_v1")
        self.assertEqual(self.payload["packet_id"], "P2391")
        self.assertEqual(self.payload["stage_id"], "S1341")
        self.assertEqual(self.payload["result_kind"], "SELECTOR_EPOCH_REBASED_BRIDGE_GAP_MATRIX")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_rg_nonduplication_and_epoch_delta(self) -> None:
        audit = self.probe["rg_nonduplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertGreater(audit["patterns"]["P2390|S1340|selector-qualified"]["count"], 0)
        self.assertGreater(audit["patterns"]["selector_source_gap_remains|bit_source_not_exported"]["count"], 0)
        self.assertIn("rebases the matrix by epoch", audit["finding"])
        self.assertEqual(self.probe["selector_hamming_delta_old_to_rebased"], 1)
        self.assertEqual(self.probe["old_selector_vector"], [1, 0, 1, 0, 0])
        self.assertEqual(self.probe["new_selector_vector"], [1, 0, 1, 1, 0])

    def test_chi11_selector_rebased_but_transport_blocked(self) -> None:
        chi11 = self.rows["topological_phase_bit_chi11"]
        self.assertTrue(chi11["epoch_delta"]["selector_bit_flipped_by_P1343_P1348"])
        self.assertTrue(chi11["rebased_bits"]["generic_strict_selector_exported"])
        self.assertFalse(chi11["rebased_bits"]["explicit_beta_tors_to_chi11_transport"])
        self.assertFalse(chi11["rebased_bits"]["role_transfer_allowed_now"])
        self.assertEqual(
            chi11["status"],
            "generic_selector_rebased_present__beta_tors_transport_and_role_transfer_absent",
        )

    def test_matrix_rank_fingerprint_and_gatekeepers(self) -> None:
        self.assertEqual(self.probe["matrix_columns"], [
            "finite_certificate_exported",
            "generic_strict_selector_exported",
            "explicit_beta_tors_to_chi11_transport",
            "strict_dynamic_source_exported",
            "role_transfer_allowed_now",
        ])
        self.assertEqual(self.probe["matrix_rank_mod2"], 2)
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("No explicit beta_tors -> chi11 transport theorem is exported.", theorem["not_licensed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))


if __name__ == "__main__":
    unittest.main()
