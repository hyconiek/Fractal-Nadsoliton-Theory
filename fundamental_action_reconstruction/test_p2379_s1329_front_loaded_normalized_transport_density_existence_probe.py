from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2379_s1329_front_loaded_normalized_transport_density_existence_probe.py"
OUT = ROOT / "generated" / "p2379_s1329_front_loaded_normalized_transport_density_existence_probe.json"
MD = ROOT / "generated" / "p2379_s1329_front_loaded_normalized_transport_density_existence_probe.md"

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
    ROOT / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.py",
    ROOT / "p2378_s1328_unit_normalized_transport_coupling_insufficiency_theorem.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2379FrontLoadedNormalizedTransportDensityExistenceProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["front_loaded_normalized_transport_density_existence_probe"]
        cls.certificate = cls.probe["front_loaded_normalized_transport_density_certificate"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2379_s1329_v1")
        self.assertEqual(self.payload["packet_id"], "P2379")
        self.assertEqual(self.payload["stage_id"], "S1329")
        self.assertEqual(self.payload["result_kind"], "FRONT_LOADED_NORMALIZED_TRANSPORT_DENSITY_EXISTENCE_PROBE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_closed_form_density_and_lattice_threshold(self) -> None:
        identity = self.certificate["weighted_transport_identity"]
        self.assertIn("int_0^1 rho_lambda(s) ds = 1", identity["density_mass"])
        self.assertIn("rho_lambda(s)>=0", identity["density_positivity"])
        self.assertIn("C(d)+lambda*B(d)", identity["closed_form"])
        lattice = self.certificate["lattice_threshold_audit"]
        self.assertEqual(lattice["lattice_points_per_axis"], 81)
        self.assertTrue(lattice["all_lattice_thresholds_between_0_and_1"])
        self.assertTrue(lattice["lambda_test_exceeds_lattice_max"])
        self.assertAlmostEqual(lattice["max_threshold_row"]["lambda_threshold_gt"], 0.7916644842269442)
        self.assertEqual(lattice["max_threshold_row"]["eta"], 1.8)
        self.assertEqual(lattice["max_threshold_row"]["beta_tors"], 0.1)

    def test_uniform_failure_front_loaded_success_and_integral_replay(self) -> None:
        rows = self.certificate["sample_front_loaded_support_audits"]
        self.assertEqual(len(rows), 9)
        self.assertLess(self.certificate["max_midpoint_integral_abs_error"], 2e-7)
        for row in rows:
            self.assertEqual(row["rho_total_mass"], 1.0)
            self.assertGreaterEqual(row["rho_min_at_s1"], 0.0)
            self.assertLess(row["lambda_threshold_gt"], row["lambda_test"])
            self.assertFalse(row["uniform_score_audit"]["d5_chamber"])
            self.assertEqual(row["uniform_score_audit"]["maximizer_h1_h5_pair_distribution"], {"3,3": 24})
            self.assertTrue(row["front_loaded_score_audit"]["d5_chamber"])
            self.assertEqual(row["front_loaded_score_audit"]["maximizer_h1_h5_pair_distribution"], {"0,4": 12})

    def test_gatekeepers_fingerprint_limits_and_next_step(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("P2378 rules out the unit-uniform endpoint primitive", theorem["claim"])
        self.assertIn("strict variational source theorem deriving the front-loaded density rho_lambda", theorem["not_licensed"])
        self.assertIn("derive either super-unit mass or a sufficiently front-loaded normalized density", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
