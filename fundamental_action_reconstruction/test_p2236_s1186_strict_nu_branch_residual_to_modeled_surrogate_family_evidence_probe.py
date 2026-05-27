from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2236S1186StrictNuBranchResidualToModeledSurrogateFamilyEvidenceProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2236_s1186_strict_nu_branch_residual_to_modeled_surrogate_family_evidence_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2236_s1186_strict_nu_branch_residual_to_modeled_surrogate_family_evidence_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2236_s1186_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["surrogate_family_evidence_exported"])
        self.assertTrue(g["affine_derivative_positive"])
        self.assertTrue(g["quadratic_derivative_positive_on_unit_interval"])


if __name__ == "__main__":
    unittest.main()
