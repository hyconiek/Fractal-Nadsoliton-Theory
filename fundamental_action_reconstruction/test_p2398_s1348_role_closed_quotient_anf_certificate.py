from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2398_s1348_role_closed_quotient_anf_certificate.py"
OUT = ROOT / "generated" / "p2398_s1348_role_closed_quotient_anf_certificate.json"
MD = ROOT / "generated" / "p2398_s1348_role_closed_quotient_anf_certificate.md"
PREREQ = ROOT / "p2397_s1347_role_closed_toe_projection_certificate.py"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2398RoleClosedQuotientANFCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(PREREQ)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["role_closed_quotient_anf_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2398_s1348_v1")
        self.assertEqual(self.payload["packet_id"], "P2398")
        self.assertEqual(self.payload["stage_id"], "S1348")
        self.assertEqual(self.payload["result_kind"], "ROLE_CLOSED_QUOTIENT_ANF_CERTIFICATE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_quotient_anf_shapes(self) -> None:
        rows = self.cert["quotient_anf"]["target_rows"]
        bridge = rows["bridge_theorem_level_closure"]
        selector = rows["selector_qw2191_closure"]
        role = rows["role_transfer_theorem_level_closure"]
        toe = rows["toe_closure"]
        self.assertEqual(bridge["anf_monomial_count"], 1)
        self.assertEqual(bridge["anf_degree"], 3)
        self.assertIn("strict_damping_beta_eta_source", bridge["anf_monomials"][0]["monomial"])
        self.assertEqual(selector["anf_monomial_count"], 1)
        self.assertEqual(selector["anf_degree"], 1)
        self.assertEqual(selector["anf_monomials"][0]["monomial"], "chi11_selector_source")
        self.assertTrue(role["is_zero_polynomial"])
        self.assertTrue(toe["is_zero_polynomial"])

    def test_minimal_supports_and_truth_vectors(self) -> None:
        rows = self.cert["quotient_anf"]["target_rows"]
        self.assertEqual(rows["bridge_theorem_level_closure"]["minimal_true_supports"], [[
            "strict_dynamical_source_for_A_P_D",
            "strict_phase_frequency_source",
            "strict_damping_beta_eta_source",
        ]])
        self.assertEqual(rows["selector_qw2191_closure"]["minimal_true_supports"], [["chi11_selector_source"]])
        self.assertEqual(rows["role_transfer_theorem_level_closure"]["truth_vector_by_mask_0_to_15"], [0] * 16)
        self.assertEqual(rows["toe_closure"]["truth_vector_by_mask_0_to_15"], [0] * 16)

    def test_fingerprint_gatekeepers_docs_and_limits(self) -> None:
        self.assertEqual(self.cert["theorem_fingerprint_sha256"], sha256_json(self.theorem))
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("No role-transfer theorem is available on the role-closed quotient.", self.theorem["not_licensed"])
        self.assertIn("No forever impossibility of future role-successor evidence is claimed.", self.theorem["not_licensed"])
        self.assertIn("role-closed quotient ANF", MD.read_text(encoding="utf-8"))
        self.assertIn("P2398/S1348 role-closed", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2398/S1348 quotient", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
