from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_kernel_completion_ladder_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_kernel_completion_ladder_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_kernel_completion_ladder_report.md"


class StrictKernelCompletionLadderProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_constants_and_stage_definitions(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_KERNEL_COMPLETION_LADDER__STAGEWISE_CERTIFICATE_NOT_DERIVATION",
        )
        self.assertEqual(payload["status"], "legacy-carrier-completed-to-current-strict-kernel-by-explicit-stage-ladder")
        constants = payload["constants"]
        self.assertAlmostEqual(constants["alpha_geo"], 2.772588722239781)
        self.assertEqual(constants["domain_d_values"], list(range(12)))
        self.assertEqual(constants["target_q_power"], "256/243")
        self.assertEqual(constants["target_eta"], "9/5")
        stages = {row["stage"]: row for row in payload["stage_definitions"]}
        self.assertEqual(
            set(stages),
            {
                "stage0_legacy_full",
                "stage1_alpha_removed_legacy_carrier",
                "stage2_phase_frequency_transported",
                "stage3_damping_compressed_strict_kernel",
            },
        )
        self.assertIn("alpha_geo", stages["stage0_legacy_full"]["formula"])
        self.assertIn("alpha_geo^{-1}", stages["stage1_alpha_removed_legacy_carrier"]["formula"])
        self.assertIn("cos(omega_S", stages["stage2_phase_frequency_transported"]["formula"])
        self.assertIn("beta*d^eta", stages["stage3_damping_compressed_strict_kernel"]["formula"])

    def test_completion_summary_and_rows(self):
        summary = self.payload["completion_summary"]
        self.assertLess(summary["max_abs_final_residual"], 1e-15)
        self.assertTrue(summary["final_residual_tolerance_pass"])
        self.assertAlmostEqual(summary["amplitude_factor"], 0.36067376022224085)
        self.assertTrue(summary["damping_factor_positive"])
        self.assertTrue(summary["damping_factor_nonincreasing"])
        self.assertAlmostEqual(summary["damping_factor_d0"], 1.0)
        self.assertAlmostEqual(summary["damping_factor_d11"], 0.014623674671644927)
        self.assertEqual(summary["phase_factor_sign_pattern"], [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1])
        self.assertGreater(summary["stage_l1_delta_alpha"], 0.0)
        self.assertGreater(summary["stage_l1_delta_phase"], 0.0)
        self.assertGreater(summary["stage_l1_delta_damping"], 0.0)
        self.assertGreater(summary["stage_norms"]["stage0_legacy_full_l2"], summary["stage_norms"]["stage3_strict_l2"])

        rows = self.payload["completion_ladder_rows"]
        self.assertEqual(len(rows), 12)
        self.assertEqual(rows[0]["d"], 0)
        self.assertEqual(rows[-1]["d"], 11)
        for row in rows:
            self.assertAlmostEqual(row["stage1_alpha_removed_legacy_carrier"], row["amplitude_factor_alpha_removal"] * row["stage0_legacy_full"], delta=1e-15)
            self.assertAlmostEqual(row["stage2_phase_frequency_transported"], row["phase_frequency_transport_factor"] * row["stage1_alpha_removed_legacy_carrier"], delta=1e-15)
            self.assertAlmostEqual(row["stage3_damping_compressed_strict_kernel"], row["damping_compression_factor"] * row["stage2_phase_frequency_transported"], delta=1e-15)
            self.assertAlmostEqual(row["stage3_damping_compressed_strict_kernel"], row["strict_kernel_target"], delta=1e-15)
            self.assertGreater(row["damping_compression_factor"], 0.0)

    def test_blocker_context_proof_guardrails_and_markdown(self):
        blockers = self.payload["blocker_context"]
        self.assertIn("strict-kernel-current-full-form", blockers["factorization_status"])
        self.assertIn("minimal-premise-lattice", blockers["blocker_lattice_status"])
        self.assertIn("orientation_chi11_source", blockers["still_open_for_theorem_level_bridge"])
        self.assertIn("orientation_chi11_source", blockers["main_one_bit_blocker"])
        self.assertIn("strict_transport_derivation", blockers["transport_blocker"])

        proof = self.payload["exact_proof_certificate"]
        self.assertIn("stage3_damping_compressed_strict_kernel equals K_strict_gate", proof["stage_identity"])
        self.assertIn("alpha removal", proof["how_legacy_is_completed"])
        self.assertIn("not another orbit", proof["nonduplication"])
        self.assertIn("not why nadsoliton dynamics", proof["theoretical_limit"])

        ontology = self.payload["ontology_guardrail"]
        self.assertIn("nadsoliton itself", ontology["allowed_reading"])
        self.assertIn("No separate informational layer", ontology["forbidden_reading"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("K_strict_gate is the current live/full kernel", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No proof derives the completion factors", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No QW-2191 discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
