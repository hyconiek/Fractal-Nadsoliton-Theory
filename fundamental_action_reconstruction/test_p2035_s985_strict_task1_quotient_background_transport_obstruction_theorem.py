from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent

P1676 = ROOT / "p1676_renormalization_background_independence_integration_matrix.py"
P1678 = ROOT / "p1678_ur_plus_background_independence_globalization.py"
P1879 = ROOT / "p1879_s829_strict_frw_to_bianchiI_transport_contract_probe.py"
P1882 = ROOT / "p1882_s832_strict_kernel_to_qg_closure_obligation_ledger_probe.py"
P1935 = ROOT / "p1935_s885_strict_po3_epsilon_bound_and_scheme_transfer_attempt_probe.py"
P1848 = ROOT / "p1848_s798_strict_gravity_componentwise_variation_and_counterterm_witness_checkpoint.py"
P2027 = ROOT / "p2027_s977_strict_b1_rank3_gb_null_adaptive_quadrature_witness.py"
P2028 = ROOT / "p2028_s978_strict_b1_gb_quotient_counterterm_identifiability_theorem.py"
P2029 = ROOT / "p2029_s979_strict_task1_renormalization_quotient_ledger_update.py"
P2030 = ROOT / "p2030_s980_strict_tensor_resolved_projection_source_audit.py"
P2031 = ROOT / "p2031_s981_strict_b1_tensor_component_profile_table_scaffold.py"
P2032 = ROOT / "p2032_s982_strict_b1_metric_gauge_component_projection_rule_audit.py"
P2033 = ROOT / "p2033_s983_strict_curved_b1_metric_ansatz_nonavailability_theorem.py"
P2034 = ROOT / "p2034_s984_strict_task1_quotient_only_renormalization_theorem.py"
P2035 = ROOT / "p2035_s985_strict_task1_quotient_background_transport_obstruction_theorem.py"

OUT2035 = ROOT / "generated" / "p2035_s985_strict_task1_quotient_background_transport_obstruction_theorem.json"


class TestStrictTask1QuotientBackgroundTransportObstructionTheorem(unittest.TestCase):
    def test_p2035_blocks_background_globalization_of_local_b1_quotient(self) -> None:
        for script in (
            P1676,
            P1678,
            P1879,
            P1882,
            P1935,
            P1848,
            P2027,
            P2028,
            P2029,
            P2030,
            P2031,
            P2032,
            P2033,
            P2034,
            P2035,
        ):
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT2035.read_text(encoding="utf-8"))

        self.assertEqual(data["schema_version"], "p2035_s985_v1")
        self.assertEqual(data["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "OPEN_LOCAL_B1_QUOTIENT_BACKGROUND_TRANSPORT_OBSTRUCTION_WITH_TRACE",
        )
        self.assertEqual(
            data["obstruction_verdict"],
            "PASS_CURRENT_EXPORT_NONTRANSPORTABILITY_WITH_TRACE",
        )

        self.assertFalse(data["grep_audit_decision"]["p2035_already_done_before_this_packet"])
        self.assertEqual(
            data["professor_decision"]["decision"],
            "P2035_BACKGROUND_TRANSPORT_OBSTRUCTION_DO_NOT_GLOBALIZE_LOCAL_B1_QUOTIENT",
        )
        self.assertEqual(
            data["professor_decision"]["rejected_branch_for_now"],
            "silently transport [a]_B1 to FRW/BianchiI or global atlas data",
        )

        theorem = data["background_transport_obstruction_theorem"]
        self.assertEqual(theorem["status"], "PASS_WITH_TRACE")
        self.assertTrue(theorem["not_a_no_go_theorem"])
        self.assertEqual(theorem["source_null_direction_R2_Ric2_Riem2_GB"], [1.0, -4.0, 1.0, -1.0])
        self.assertIn("cannot be promoted", theorem["statement"])

        contract = data["required_transport_contract"]
        self.assertEqual(contract["contract_id"], "Task1_B1_quotient_background_transport_contract_v1")
        self.assertIn("B1_scalar", contract["background_family_set"])
        self.assertIn("FRW", contract["background_family_set"])
        self.assertIn("BianchiI", contract["background_family_set"])
        self.assertIn(
            "basis_preserving_quotient_map_between_backgrounds",
            contract["required_missing_exports"],
        )
        self.assertIn(
            "finite_part_scheme_transport_map_on_same_operator_basis",
            contract["required_missing_exports"],
        )
        self.assertIn(
            "global_atlas_cocycle_for_renormalized_quotient_data",
            contract["required_missing_exports"],
        )

        sources = {row["source"] for row in data["obstruction_trace"]}
        for source in {"P2034", "P2033", "P1879", "P1882", "P1935", "P1963", "P1676", "P1678"}:
            self.assertIn(source, sources)

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["local_b1_quotient_pass_from_p2034"])
        self.assertTrue(checks["p2034_no_background_globalization_claimed"])
        self.assertTrue(checks["p2034_no_tensor_projection_claimed"])
        self.assertTrue(checks["curved_b1_component_rule_unavailable"])
        self.assertTrue(checks["p1879_transport_contract_open"])
        self.assertTrue(checks["p1882_b1_background_independence_open"])
        self.assertTrue(checks["p1935_scheme_transfer_partial"])
        self.assertTrue(checks["p1963_global_background_independence_open"])
        self.assertTrue(checks["p1676_background_gate_theorems_missing"])
        self.assertTrue(checks["p1678_globalization_theorems_missing"])
        self.assertFalse(checks["basis_preserving_quotient_transport_map_exported"])
        self.assertFalse(checks["finite_part_scheme_transport_map_exported"])
        self.assertFalse(checks["gb_null_direction_transport_exported"])
        self.assertFalse(checks["background_globalization_exported"])
        self.assertTrue(checks["obstruction_theorem_pass"])
        self.assertTrue(checks["no_global_renormalization_claimed"])
        self.assertTrue(checks["no_tensor_projection_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])

        not_licensed = data["theorem_scope"]["not_licensed"]
        self.assertIn("FRW/BianchiI/global-atlas renormalization of the P2034 quotient class", not_licensed)
        self.assertIn("tensor-component B1 renormalization", not_licensed)
        self.assertIn("ToE closure", not_licensed)


if __name__ == "__main__":
    unittest.main()
