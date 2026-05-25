from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
P2037 = ROOT / "p2037_s987_strict_task1_same_scheme_finite_part_map_seed.py"
P2038 = ROOT / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.py"
P2039 = ROOT / "p2039_s989_strict_same_scheme_finite_part_candidate_uncertainty_bound_probe.py"
P2040 = ROOT / "p2040_s990_strict_same_scheme_subtraction_compatibility_residual_audit.py"
P2041 = ROOT / "p2041_s991_strict_same_scheme_operator_level_subtraction_identity_witness.py"
P2042 = ROOT / "p2042_s992_strict_same_scheme_operator_level_residual_envelope_sweep.py"
P2043 = ROOT / "p2043_s993_strict_same_scheme_stratified_residual_envelope_audit.py"
OUT2043 = ROOT / "generated" / "p2043_s993_strict_same_scheme_stratified_residual_envelope_audit.json"


class TestP2043StratifiedResidualEnvelopeAudit(unittest.TestCase):
    def test_p2043_exports_stratified_envelope_and_keeps_c3_open(self) -> None:
        subprocess.run([sys.executable, str(P2037)], check=True)
        subprocess.run([sys.executable, str(P2038)], check=True)
        subprocess.run([sys.executable, str(P2039)], check=True)
        subprocess.run([sys.executable, str(P2040)], check=True)
        subprocess.run([sys.executable, str(P2041)], check=True)
        subprocess.run([sys.executable, str(P2042)], check=True)
        subprocess.run([sys.executable, str(P2043)], check=True)

        data = json.loads(OUT2043.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2043_s993_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_STRATIFIED_RESIDUAL_ENVELOPE_AUDIT_WITH_BOOTSTRAP_CI__C3_STILL_OPEN",
        )

        strat = data["stratified_envelope"]
        self.assertEqual(len(strat["per_bucket"]), 4)
        self.assertEqual(len(strat["norm_bucket_order"]), 4)
        self.assertEqual(len(strat["sup_sequence"]), 4)
        self.assertGreater(strat["global_sup_residual_with_buffer_linf"], 0.0)

        bci = data["bootstrap_ci"]
        self.assertEqual(bci["bootstrap_seed"], 2043)
        self.assertEqual(bci["replicates"], 600)
        self.assertLessEqual(bci["ci_low"], bci["ci_high"])

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["per_bucket_nonempty"])
        self.assertTrue(checks["bootstrap_ci_finite"])
        self.assertTrue(checks["global_sup_finite"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])

        c3 = data["c3_gate_update"]
        self.assertEqual(c3["C3_transport_theorem"], "OPEN")
        self.assertEqual(c3["C3_discharge_status"], "NOT_DISCHARGED")


if __name__ == "__main__":
    unittest.main()
