from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2393_s1343_kernel_completion_boundary_residual_certificate.py"
OUT = ROOT / "generated" / "p2393_s1343_kernel_completion_boundary_residual_certificate.json"
MD = ROOT / "generated" / "p2393_s1343_kernel_completion_boundary_residual_certificate.md"
PREREQ = ROOT / "p2392_s1342_auxiliary_beta_tors_chi11_retirement_certificate.py"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2393KernelCompletionBoundaryResidualCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(PREREQ)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["kernel_completion_boundary_residual_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2393_s1343_v1")
        self.assertEqual(self.payload["packet_id"], "P2393")
        self.assertEqual(self.payload["stage_id"], "S1343")
        self.assertEqual(self.payload["result_kind"], "KERNEL_COMPLETION_BOUNDARY_RESIDUAL_CERTIFICATE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_boundary_identity_is_exact_finite_replay(self) -> None:
        boundary = self.cert["boundary_identity_audit"]
        self.assertTrue(boundary["allclose_at_1e_minus_15"])
        self.assertLess(boundary["max_abs_error_float_replay"], 1.0e-15)
        self.assertEqual(boundary["boundary_parameter_substitution"]["eta"], 1.0)
        self.assertEqual(len(boundary["rows"]), 13)

    def test_current_target_has_nonzero_completion_residual(self) -> None:
        residual = self.cert["current_strict_residual_audit"]
        self.assertTrue(residual["nonzero_residual_on_current_target"])
        self.assertGreater(residual["l2_residual"], 0.0)
        self.assertGreater(residual["linf_residual_row"]["abs_residual"], 0.0)
        self.assertIn("phase/frequency", residual["interpretation"])
        self.assertEqual(residual["strict_target_parameters"]["eta"], 9.0 / 5.0)

    def test_obligation_matrix_keeps_bridge_and_role_transfer_open(self) -> None:
        matrix = self.cert["bridge_obligation_matrix"]
        self.assertEqual(matrix["rank_mod2"], 3)
        self.assertTrue(matrix["all_current_role_transfer_bits_zero"])
        self.assertIn("phase/frequency passage legacy canonical to strict target", matrix["open_completion_rows"])
        self.assertIn("linear torsion damping to nonlinear strict compression", matrix["open_completion_rows"])
        statuses = {row["obligation"]: row["status"] for row in matrix["rows"]}
        self.assertEqual(
            statuses["selector route after P2392"],
            "SELECTOR_MECHANISM_AVAILABLE_BETA_TORS_CHI11_ROUTE_RETIRED",
        )

    def test_fingerprint_gatekeepers_and_docs(self) -> None:
        self.assertEqual(self.cert["theorem_fingerprint_sha256"], sha256_json(self.theorem))
        self.assertIn("No beta_tors -> chi11 selector target is reopened.", self.theorem["not_licensed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        md = MD.read_text(encoding="utf-8")
        self.assertIn("normalized eta=1 boundary embedding", md)
        eq = (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8")
        lag = (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8")
        self.assertIn("P2393/S1343 normalized kernel boundary", eq)
        self.assertIn("P2393/S1343 normalized kernel boundary residual", lag)


if __name__ == "__main__":
    unittest.main()
