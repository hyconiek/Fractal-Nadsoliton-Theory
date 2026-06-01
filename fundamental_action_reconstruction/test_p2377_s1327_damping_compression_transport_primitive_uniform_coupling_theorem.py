from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.py"
OUT = ROOT / "generated" / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.json"
MD = ROOT / "generated" / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.md"

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
    ROOT / "p2375_s1325_damping_compression_polarity_interval_robustness_theorem.py",
    ROOT / "p2376_s1326_damping_compression_eta_beta_rectangle_robustness_theorem.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2377DampingCompressionTransportPrimitiveUniformCouplingTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["damping_compression_transport_primitive_uniform_coupling_theorem"]
        cls.certificate = cls.probe["transport_primitive_uniform_coupling_certificate"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2377_s1327_v1")
        self.assertEqual(self.payload["packet_id"], "P2377")
        self.assertEqual(self.payload["stage_id"], "S1327")
        self.assertEqual(self.payload["result_kind"], "DAMPING_COMPRESSION_TRANSPORT_PRIMITIVE_UNIFORM_COUPLING_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_transport_primitive_identity_and_uniform_threshold(self) -> None:
        identity = self.certificate["transport_identity"]
        threshold = self.certificate["uniform_threshold_certificate"]
        self.assertIn("integral_0^1", identity["primitive"])
        self.assertLess(self.certificate["max_midpoint_integral_abs_error"], 1e-8)
        self.assertEqual(threshold["denominator"]["minimum_corner"], {"eta": 1.8, "beta_tors": 0.1})
        self.assertGreater(threshold["denominator"]["minimum_corner_value"], 0)
        self.assertGreater(threshold["denominator"]["minimum_corner_dD_deta"], 0)
        self.assertLessEqual(threshold["denominator"]["right_endpoint_dD_dx"], 0)
        self.assertGreater(threshold["uniform_coupling"]["tau_gt_uniform"], 0)

    def test_grid_blend_support_scans(self) -> None:
        rows = self.certificate["sample_transport_and_support_audits"]
        self.assertEqual(len(rows), 9)
        threshold = self.certificate["uniform_threshold_certificate"]["uniform_coupling"]["tau_gt_uniform"]
        canonical = [row for row in rows if row["eta"] == 1.8 and row["beta_tors"] == 0.01][0]
        self.assertGreater(threshold, 1.6259309081595656)
        self.assertLess(canonical["blend_denominator_C5_minus_3C1"], rows[-1]["blend_denominator_C5_minus_3C1"])
        for row in rows:
            self.assertGreater(row["blend_denominator_C5_minus_3C1"], 0)
            self.assertTrue(row["blended_score_audit"]["d5_chamber"])
            self.assertEqual(row["blended_score_audit"]["maximizer_h1_h5_pair_distribution"], {"0,4": 12})

    def test_gatekeepers_fingerprint_limits_and_next_step(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("endpoint primitive", theorem["claim"])
        self.assertIn("strict variational source theorem fixing the scalar coupling tau", theorem["not_licensed"])
        self.assertIn("fixes a coupling tau above the threshold", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
