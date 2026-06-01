from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2384_s1334_symbolic_bathtub_inequality_proof_packet.py"
OUT = ROOT / "generated" / "p2384_s1334_symbolic_bathtub_inequality_proof_packet.json"
MD = ROOT / "generated" / "p2384_s1334_symbolic_bathtub_inequality_proof_packet.md"

PREREQ_SCRIPTS = [
    ROOT / "p2383_s1333_closed_form_bathtub_corner_reduction_theorem.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2384SymbolicBathtubInequalityProofPacketTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["symbolic_bathtub_inequality_proof_packet"]
        cls.ratio = cls.probe["symbolic_ratio_proof"]
        cls.eta = cls.probe["eta_derivative_sign_proof"]
        cls.beta = cls.probe["beta_tors_derivative_sign_proof"]
        cls.cap = cls.probe["cap_derivative_sign_proof"]
        cls.target = cls.probe["source_target_translation"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2384_s1334_v1")
        self.assertEqual(self.payload["packet_id"], "P2384")
        self.assertEqual(self.payload["stage_id"], "S1334")
        self.assertEqual(self.payload["result_kind"], "SYMBOLIC_BATHTUB_INEQUALITY_PROOF_PACKET")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_symbolic_ratio_proof_margin(self) -> None:
        self.assertTrue(self.ratio["passes"])
        self.assertGreater(float(self.ratio["sqrt125_gt_threshold"]["margin"]), 4.7)
        self.assertEqual(
            self.ratio["sqrt125_gt_threshold"]["integer_certificate"]["integer_margin_after_squaring"],
            10384,
        )
        self.assertGreater(float(self.ratio["rmin_lower_square_minus_3"]), 0.0)
        self.assertIn("without a grid-first proof premise", self.ratio["logical_chain"][-1])

    def test_derivative_sign_proofs_have_positive_coarse_gaps(self) -> None:
        self.assertTrue(self.eta["passes"])
        self.assertIn("dW1/deta=0", self.eta["d1_term"])
        self.assertTrue(self.beta["passes"])
        self.assertGreater(float(self.beta["positive_gap_N5_minus_3N1"]), 0.7)
        self.assertIn("N5>3*N1", self.beta["conclusion"])
        self.assertTrue(self.cap["passes"])
        self.assertGreater(float(self.cap["positive_gap_h5_minus_3h1"]), 0.54)
        self.assertIn("h'(x)=x/(1+x)^2>0", self.cap["monotonic_h_note"])

    def test_source_target_and_gatekeepers(self) -> None:
        self.assertEqual(self.target["early_interval_length"], "0.625")
        self.assertEqual(self.target["early_half_mass"], "0.8")
        self.assertEqual(self.target["barycenter"], "0.3125")
        self.assertEqual(self.target["barycenter_shift_from_uniform"], "0.1875")
        self.assertIn("non-strict premise", self.target["source_theorem_obligation"])
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("symbolic inequality core", theorem["claim"])
        self.assertIn("strict source theorem deriving the cap M or bang-bang density", theorem["not_licensed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))


if __name__ == "__main__":
    unittest.main()
