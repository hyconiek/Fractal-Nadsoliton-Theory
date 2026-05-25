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
P2036 = ROOT / "p2036_s986_strict_task1_quotient_background_transport_candidate_contract.py"

OUT2036 = ROOT / "generated" / "p2036_s986_strict_task1_quotient_background_transport_candidate_contract.json"


class TestStrictTask1QuotientBackgroundTransportCandidateContract(unittest.TestCase):
    def test_p2036_exports_only_nonproof_transport_contract_candidate(self) -> None:
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
            P2036,
        ):
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT2036.read_text(encoding="utf-8"))

        self.assertEqual(data["schema_version"], "p2036_s986_v1")
        self.assertEqual(data["status"], "OPEN_CANDIDATE_CONTRACT_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "OPEN_NEW_PREMISE_BACKGROUND_TRANSPORT_CONTRACT_CANDIDATE_EXPORTED__NO_GLOBALIZATION",
        )
        self.assertEqual(
            data["professor_decision"]["decision"],
            "EXPORT_NONSTRICT_BACKGROUND_TRANSPORT_CONTRACT_CANDIDATE_KEEP_TASK1_LOCAL",
        )
        self.assertFalse(data["grep_audit_decision"]["p2036_already_done_before_this_packet"])

        contract = data["candidate_contract"]
        self.assertEqual(contract["contract_id"], "Task1_B1_quotient_background_transport_candidate_v1")
        self.assertEqual(contract["contract_status"], "CANDIDATE_EXPORTED_NOT_A_THEOREM")
        self.assertEqual(contract["new_premise_status"], "NEW_PREMISE_CANDIDATE_NOT_ADMITTED_STRICT_PROOF")
        self.assertEqual(contract["quotient_basis"], ["R2_bar", "Ric2_bar", "Riem2_bar"])
        self.assertEqual(contract["source_null_direction_R2_Ric2_Riem2_GB"], [1.0, -4.0, 1.0, -1.0])
        self.assertEqual(
            contract["quotient_projection_matrix_Q_B1_from_R4_to_R3"],
            [[1.0, 0.0, 0.0, 1.0], [0.0, 1.0, 0.0, -4.0], [0.0, 0.0, 1.0, 1.0]],
        )
        self.assertEqual(contract["zeroth_order_transport_seed"]["identity_seed_rank"], 3)
        self.assertEqual(contract["zeroth_order_transport_seed"]["identity_seed_determinant"], 1.0)
        self.assertIn("not a transport theorem", contract["zeroth_order_transport_seed"]["warning"])

        gate_status = {row["gate_id"]: row["current_status"] for row in data["acceptance_gates"]}
        self.assertEqual(gate_status["C1_QUOTIENT_SOURCE"], "SATISFIED_BY_P2034_LOCAL_SOURCE")
        self.assertEqual(gate_status["C2_BASIS_PRESERVING_MAP"], "OPEN_SYMBOLIC_PLACEHOLDER")
        self.assertEqual(gate_status["C3_FINITE_PART_SCHEME_TRANSPORT"], "OPEN_SYMBOLIC_PLACEHOLDER")
        self.assertEqual(gate_status["C4_GB_NULL_TRANSPORT"], "OPEN_TOPOLOGICAL_CLASSIFICATION_MISSING")
        self.assertEqual(gate_status["C5_ANISOTROPIC_WITNESS_DATA"], "OPEN_LOOP_DATA_MISSING")
        self.assertEqual(gate_status["C6_GLOBAL_COCYCLE"], "OPEN_COCYCLE_THEOREM_MISSING")

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["p2034_local_source_ready"])
        self.assertTrue(checks["p2035_obstruction_ready"])
        self.assertTrue(checks["p1879_anisotropy_contract_available"])
        self.assertTrue(checks["p1935_scheme_transfer_still_partial"])
        self.assertTrue(checks["candidate_contract_syntactically_complete"])
        self.assertTrue(checks["identity_seed_rank3"])
        self.assertFalse(checks["all_acceptance_gates_passed"])
        self.assertFalse(checks["strict_transport_theorem_exported"])
        self.assertFalse(checks["background_globalization_exported"])
        self.assertFalse(checks["finite_part_scheme_transport_proven"])
        self.assertFalse(checks["gb_null_transport_classified"])
        self.assertFalse(checks["global_cocycle_proven"])
        self.assertTrue(checks["no_global_renormalization_claimed"])
        self.assertTrue(checks["no_tensor_projection_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])

        not_licensed = data["theorem_scope"]["not_licensed"]
        self.assertIn("background-global Task-1 renormalization", not_licensed)
        self.assertIn("finite-part scheme transport theorem", not_licensed)
        self.assertIn("ToE closure", not_licensed)


if __name__ == "__main__":
    unittest.main()
