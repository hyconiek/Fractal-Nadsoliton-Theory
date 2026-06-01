from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2375_s1325_damping_compression_polarity_interval_robustness_theorem.py"
OUT = ROOT / "generated" / "p2375_s1325_damping_compression_polarity_interval_robustness_theorem.json"
MD = ROOT / "generated" / "p2375_s1325_damping_compression_polarity_interval_robustness_theorem.md"
SCRATCH_DAMPING = ROOT / "scratch" / "bridge_strict_completion_legacy_to_strict_damping_compression_separation_certificate_probe.py"

PREREQ_SCRIPTS = [
    ROOT / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.py",
    ROOT / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.py",
    ROOT / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.py",
    ROOT / "p2366_s1316_selector_candidate_phase_reference_chi11_audit_probe.py",
    ROOT / "p2367_s1317_selector_phase_origin_admissibility_no_go_probe.py",
    ROOT / "p2368_s1318_self_recorded_endpoint_anchor_selector_candidate_probe.py",
    ROOT / "p2369_s1319_self_recorded_ledger_closed_form_uniqueness_theorem.py",
    ROOT / "p2370_s1320_d5_bandpass_support_closed_form_theorem.py",
    ROOT / "p2371_s1321_aut_invariant_unit_bandpass_obstruction_theorem.py",
    ROOT / "p2372_s1322_bridge_kernel_direct_band_polarity_audit.py",
    ROOT / "p2373_s1323_bridge_kernel_polarity_correction_cone_theorem.py",
    ROOT / "p2374_s1324_damping_compression_band_polarity_candidate_theorem.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2375DampingCompressionPolarityIntervalRobustnessTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRATCH_DAMPING)], check=True, cwd=ROOT.parent)
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["damping_compression_polarity_interval_robustness_theorem"]
        cls.certificate = cls.probe["robustness_certificate"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2375_s1325_v1")
        self.assertEqual(self.payload["packet_id"], "P2375")
        self.assertEqual(self.payload["stage_id"], "S1325")
        self.assertEqual(self.payload["result_kind"], "DAMPING_COMPRESSION_POLARITY_INTERVAL_ROBUSTNESS_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_interval_proof_certificate(self) -> None:
        symbolic = self.certificate["symbolic_inequality"]
        proof = self.certificate["proof_summary"]
        self.assertGreater(symbolic["margin_min_at_zero"], 0)
        self.assertGreater(symbolic["derivative_min_at_zero"], 0)
        self.assertTrue(proof["derivative_positive_on_interval"])
        self.assertTrue(proof["margin_positive_on_interval"])
        self.assertTrue(proof["therefore_C1_over_C5_lt_one_third_for_all_interval"])

    def test_sample_support_scans(self) -> None:
        rows = self.certificate["sample_score_audits"]
        self.assertEqual([row["beta_tors"] for row in rows], [0.0, 0.01, 0.05, 0.1])
        for row in rows:
            self.assertLess(row["C1_over_C5"], 1 / 3)
            self.assertGreater(row["margin_C5_minus_3C1"], 0)
            self.assertEqual(row["score_audit"]["maximizer_h1_h5_pair_distribution"], {"0,4": 12})

    def test_gatekeepers_fingerprint_limits_and_next_step(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("not a fine-tuned artifact", theorem["claim"])
        self.assertIn("QW-2191 discharge", theorem["not_licensed"])
        self.assertIn("variational/transport source theorem", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
