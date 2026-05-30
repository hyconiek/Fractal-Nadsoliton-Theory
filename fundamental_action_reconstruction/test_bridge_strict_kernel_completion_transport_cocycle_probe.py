from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_kernel_completion_transport_cocycle_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_kernel_completion_transport_cocycle_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_kernel_completion_transport_cocycle_report.md"


class StrictKernelCompletionTransportCocycleProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_constants_and_definitions(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_KERNEL_COMPLETION_TRANSPORT_COCYCLE__FINITE_Z12_WITNESS_NOT_DYNAMICAL_DERIVATION",
        )
        self.assertIn("unique-finite-node-transport", payload["status"])
        constants = payload["constants"]
        self.assertAlmostEqual(constants["alpha_geo"], 2.772588722239781)
        self.assertEqual(constants["domain_d_values"], list(range(12)))
        self.assertEqual(constants["edge_d_values"], list(range(11)))
        definitions = payload["transport_definitions"]
        self.assertIn("K_strict_gate", definitions["node_transport"])
        self.assertIn("R(d->d+1)", definitions["edge_cocycle"])
        self.assertIn("Delta log|P|", definitions["additive_log_split"])
        self.assertIn("phase-sign flips", definitions["sign_cocycle"])

    def test_cocycle_summary(self):
        summary = self.payload["cocycle_summary"]
        self.assertLess(summary["max_abs_node_factorization_residual"], 1e-15)
        self.assertLess(summary["max_abs_reconstruction_residual"], 1e-15)
        self.assertLess(summary["max_abs_edge_log_split_residual"], 1e-15)
        self.assertEqual(summary["interval_count"], 66)
        self.assertLess(summary["max_abs_interval_multiplicative_residual"], 1e-15)
        self.assertLess(summary["max_abs_interval_additive_log_residual"], 1e-14)
        self.assertEqual(summary["transport_sign_pattern"], [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1])
        self.assertEqual(summary["phase_sign_flip_edges"], ["1->2", "5->6", "7->8", "9->10"])
        self.assertTrue(summary["alpha_drops_out_of_edge_cocycle"])
        self.assertAlmostEqual(summary["T0"], 0.4109835665709375)
        self.assertAlmostEqual(summary["T11"], 0.0032388039359289075)
        self.assertLess(summary["T11_over_T0"], 0.008)

    def test_nodes_edges_reconstruction_and_intervals(self):
        nodes = self.payload["node_transport_rows"]
        edges = self.payload["edge_cocycle_rows"]
        recon = self.payload["node_reconstruction_rows"]
        intervals = self.payload["all_interval_cocycle_rows"]
        self.assertEqual(len(nodes), 12)
        self.assertEqual(len(edges), 11)
        self.assertEqual(len(recon), 12)
        self.assertEqual(len(intervals), 66)
        self.assertEqual(nodes[0]["d"], 0)
        self.assertEqual(nodes[-1]["d"], 11)
        for node in nodes:
            self.assertAlmostEqual(node["transport_T"], node["factor_product_A_P_D"], delta=1e-15)
        for edge in edges:
            self.assertAlmostEqual(
                edge["log_abs_R"],
                edge["delta_log_abs_phase"] + edge["delta_log_damping"],
                delta=1e-15,
            )
        for row in recon:
            self.assertAlmostEqual(row["reconstructed_T"], row["actual_T"], delta=1e-15)
        for row in intervals:
            self.assertAlmostEqual(row["edge_product"], row["endpoint_ratio"], delta=1e-15)
            self.assertAlmostEqual(row["edge_log_sum"], row["endpoint_log_delta"], delta=1e-14)
            self.assertEqual(row["edge_sign_product"], row["endpoint_sign_ratio"])

    def test_blocker_context_proof_guardrails_and_markdown(self):
        blockers = self.payload["blocker_context"]
        self.assertIn("all-three-explicit-factors", blockers["necessity_status"])
        self.assertIn("legacy-carrier-completed", blockers["ladder_status"])
        self.assertIn("minimal-premise-lattice", blockers["blocker_lattice_status"])
        self.assertIn("still_open", blockers["strict_transport_derivation_status"])
        self.assertIn("orientation_chi11_source", blockers["still_open"])
        self.assertIn("role_transfer_theorem", blockers["still_open"])

        proof = self.payload["proof_certificate"]
        self.assertIn("unique node value", proof["uniqueness_step"])
        self.assertIn("reconstructs all T(d)", proof["edge_step"])
        self.assertIn("66 nontrivial intervals", proof["interval_step"])
        self.assertIn("alpha is constant", proof["split_step"])
        self.assertIn("four phase sign flips", proof["sign_step"])
        self.assertIn("not another subset-necessity", proof["nonduplication"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("K_strict_gate remains the current live/full", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No proof derives the transport cocycle", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
