from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2374_s1324_damping_compression_band_polarity_candidate_theorem.py"
OUT = ROOT / "generated" / "p2374_s1324_damping_compression_band_polarity_candidate_theorem.json"
MD = ROOT / "generated" / "p2374_s1324_damping_compression_band_polarity_candidate_theorem.md"

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
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2374DampingCompressionBandPolarityCandidateTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["damping_compression_band_polarity_candidate_theorem"]
        cls.certificate = cls.probe["damping_compression_certificate"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2374_s1324_v1")
        self.assertEqual(self.payload["packet_id"], "P2374")
        self.assertEqual(self.payload["stage_id"], "S1324")
        self.assertEqual(self.payload["result_kind"], "DAMPING_COMPRESSION_BAND_POLARITY_CANDIDATE_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_compression_chamber_and_blend_threshold(self) -> None:
        compression = self.certificate["compression_log_weights"]
        chamber = self.certificate["closed_form_chamber_tests"]
        self.assertLess(compression["C1_over_C5"], 1 / 3)
        self.assertTrue(chamber["compression_alone_condition_holds"])
        self.assertGreater(chamber["blend_denominator_C5_minus_3C1"], 0)
        self.assertGreater(chamber["blend_tau_gt"], 0)

    def test_score_audits(self) -> None:
        audits = self.certificate["score_audits"]
        self.assertEqual(audits["strict_direct_kernel"]["maximizer_h1_h5_pair_distribution"], {"4,0": 12})
        self.assertEqual(audits["compression_log_weight_alone"]["maximizer_h1_h5_pair_distribution"], {"0,4": 12})
        self.assertEqual(
            audits["strict_plus_minimal_compression_blend_plus_epsilon"]["maximizer_h1_h5_pair_distribution"],
            {"0,4": 12},
        )

    def test_gatekeepers_fingerprint_limits_and_recommendation(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("candidate source direction", theorem["claim"])
        self.assertIn("QW-2191 discharge", theorem["not_licensed"])
        self.assertIn("variational or transport-level theorem", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
