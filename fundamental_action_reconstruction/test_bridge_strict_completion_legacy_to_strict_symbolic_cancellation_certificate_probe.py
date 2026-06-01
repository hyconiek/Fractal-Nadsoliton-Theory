import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_symbolic_cancellation_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_symbolic_cancellation_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_legacy_to_strict_symbolic_cancellation_certificate_report.md"


class LegacyToStrictSymbolicCancellationCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_symbolic_identity(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_LEGACY_TO_STRICT_SYMBOLIC_CANCELLATION_CERTIFICATE__ANSATZ_LEVEL_ONLY",
        )
        self.assertIn("symbolic-apd-cancellation", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "necessity",
            "finite_bridge_assembly",
            "finite_diagonal_completion_map",
            "amplitude_scalar_normalization",
            "phase_frequency_affine_transport",
            "damping_parameter_identifiability",
            "legacy_bridge_guardrail",
        })
        self.assertIn("symbolic bridge cancellation", payload["grep_disambiguation"]["searched_terms"])
        self.assertIn("formula-level cancellation ledger", payload["grep_disambiguation"]["finding"])

        symbolic = payload["symbolic_identity"]
        self.assertEqual(symbolic["legacy_kernel"], "K_L=alpha*cL/Lden")
        self.assertEqual(symbolic["strict_kernel"], "K_S=cS/Sden")
        self.assertIn("(1/alpha)*(cS/cL)*(Lden/Sden)", symbolic["completion_factor"])
        self.assertEqual(symbolic["reduced_product"], "cS/Sden=K_S")
        self.assertEqual(symbolic["admissibility_conditions"], ["alpha != 0", "cL != 0", "Lden != 0", "Sden != 0"])

    def test_admissibility_cancellation_summary_and_cross_checks(self):
        payload = self.payload
        rows = payload["admissibility_rows"]
        self.assertEqual(len(rows), 12)
        self.assertEqual([row["d"] for row in rows], list(range(12)))
        self.assertTrue(all(row["alpha_geo_nonzero"] for row in rows))
        self.assertTrue(all(row["legacy_cos_nonzero"] for row in rows))
        self.assertTrue(all(row["legacy_denominator_positive"] for row in rows))
        self.assertTrue(all(row["strict_denominator_positive"] for row in rows))

        cancellation_steps = [row["step"] for row in payload["cancellation_rows"]]
        self.assertEqual(cancellation_steps, [
            "expand",
            "cancel_alpha",
            "cancel_legacy_cos",
            "cancel_legacy_denominator",
            "identify_strict_kernel",
        ])
        self.assertIn("cS/Sden", payload["cancellation_rows"][-1]["expression"])

        summary = payload["symbolic_cancellation_summary"]
        self.assertTrue(summary["symbolic_cancellation_formula_exported"])
        self.assertTrue(summary["alpha_geo_nonzero_on_domain"])
        self.assertTrue(summary["legacy_cos_nonzero_on_domain"])
        self.assertTrue(summary["legacy_denominator_positive_on_domain"])
        self.assertTrue(summary["strict_denominator_positive_on_domain"])
        self.assertTrue(summary["finite_assembly_reconstruction_inherited"])
        self.assertTrue(summary["finite_assembly_matches_diagonal_inherited"])
        self.assertTrue(summary["diagonal_uniqueness_inherited"])
        self.assertTrue(summary["phase_transport_source_still_open"])
        self.assertTrue(summary["damping_source_still_open"])
        self.assertTrue(summary["role_transfer_required_after_full_bridge"])
        self.assertFalse(summary["strict_dynamic_source_exported"])
        self.assertFalse(summary["selector_source_exported"])
        self.assertFalse(summary["legacy_role_transfer_exported"])
        self.assertFalse(summary["raw_kernel_identity_claimed"])
        self.assertFalse(summary["full_bridge_theorem_exported"])
        self.assertFalse(summary["toe_closure_claimed"])
        self.assertTrue(all(payload["cross_checks"].values()))

    def test_proof_and_hard_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("rg was used", proof["grep_step"])
        self.assertIn("cancels to cS/Sden=K_S", proof["formal_step"])
        self.assertIn("alpha_geo != 0", proof["admissibility_step"])
        self.assertIn("unique diagonal Q=diag(K_strict/K_legacy)", proof["finite_consistency_step"])
        self.assertIn("ansatz-level formula identity", proof["scope_step"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No unqualified raw identity", hard_limits)
        self.assertIn("No strict dynamical derivation", hard_limits)
        self.assertIn("No beta_tors -> beta/eta", hard_limits)
        self.assertIn("No legacy physical-role transfer", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
