from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2392_s1342_auxiliary_beta_tors_chi11_retirement_certificate.py"
OUT = ROOT / "generated" / "p2392_s1342_auxiliary_beta_tors_chi11_retirement_certificate.json"
MD = ROOT / "generated" / "p2392_s1342_auxiliary_beta_tors_chi11_retirement_certificate.md"
PREREQS = [
    ROOT / "p2390_s1340_selector_qualified_beta_tors_chi11_role_audit.py",
    ROOT / "p2391_s1341_selector_epoch_rebased_bridge_gap_matrix.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2392AuxiliaryBetaTorsChi11RetirementCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for prereq in PREREQS:
            subprocess.run([sys.executable, str(prereq)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["auxiliary_beta_tors_chi11_retirement_certificate"]
        cls.theorem = cls.probe["theorem_export"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2392_s1342_v1")
        self.assertEqual(self.payload["packet_id"], "P2392")
        self.assertEqual(self.payload["stage_id"], "S1342")
        self.assertEqual(
            self.payload["result_kind"],
            "AUXILIARY_BETA_TORS_CHI11_SELECTOR_ASSUMPTION_RETIREMENT_CERTIFICATE",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_minimal_support_retires_beta_for_selector(self) -> None:
        available = self.theorem["available_atoms"]
        self.assertTrue(available["strict_internal_selector_P1343_P1348"])
        self.assertFalse(available["auxiliary_beta_tors_to_chi11"])
        selector_support = self.theorem["support_evaluation"]["selector_mechanism"]
        self.assertTrue(selector_support["target_satisfied"])
        self.assertFalse(selector_support["uses_auxiliary_beta_tors_in_realized_support"])
        self.assertEqual(selector_support["realized_minimal_supports"], [["strict_internal_selector_P1343_P1348"]])
        self.assertEqual(self.theorem["active_beta_tors_chi11_obligation_count"], 0)
        self.assertEqual(self.probe["active_beta_tors_chi11_obligations"], [])

    def test_role_transfer_still_open_but_not_selector_blocker(self) -> None:
        rows = {row["obligation"]: row for row in self.probe["obligation_rows"]}
        self.assertEqual(
            rows["auxiliary beta_tors -> chi11 selector-search hypothesis"]["status"],
            "RETIRED_AS_SELECTOR_ROUTE_ASSUMPTION",
        )
        self.assertFalse(rows["auxiliary beta_tors -> chi11 selector-search hypothesis"]["active_next_target"])
        self.assertEqual(rows["legacy -> strict bridge completion"]["status"], "OPEN_ACTIVE_BRIDGE_TARGET")
        role_support = self.theorem["support_evaluation"]["legacy_role_transfer_permission"]
        self.assertFalse(role_support["target_satisfied"])

    def test_fingerprint_gatekeepers_and_docs(self) -> None:
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(self.theorem))
        self.assertIn("No beta_tors -> chi11 theorem is claimed or needed here.", self.theorem["not_licensed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        md = MD.read_text(encoding="utf-8")
        self.assertIn("retired from the active selector blocker list", md)


if __name__ == "__main__":
    unittest.main()
