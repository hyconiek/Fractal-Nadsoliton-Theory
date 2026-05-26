from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p1950_s900_strict_renormalization_exact_integration.py",
    ROOT / "p2094_s1044_strict_b1_quotient_renormalization_rank_repair.py",
    ROOT / "p2095_s1045_strict_b1_gb_derived_channel_certificate.py",
    ROOT / "p2096_s1046_strict_b1_renormalization_closure_contract.py",
    ROOT / "p2097_s1047_strict_b1_quotient_closure_stability_minigrid.py",
    ROOT / "p2098_s1048_strict_precutkosky_readiness_contract.py",
    ROOT / "p2099_s1049_strict_u1_same_scheme_lock_witness.py",
    ROOT / "p2100_s1050_strict_u2_phase_space_quadrature_witness.py",
    ROOT / "p2101_s1051_strict_u3_residue_positivity_uncertainty_witness.py",
    ROOT / "p2102_s1052_strict_task2_entry_gate_summary.py",
    ROOT / "p2103_s1053_strict_dressed_discontinuity_backend_import_contract.py",
    ROOT / "p2104_s1054_strict_d1_dressed_pole_residue_source_object.py",
    ROOT / "p2105_s1055_strict_d2_disc_dressed_vs_cutsum_comparator_contract.py",
    ROOT / "p2106_s1056_strict_cmp1_disc_dressed_row_or_nonavailability.py",
    ROOT / "p2107_s1057_strict_cmp1_first_residual_execution.py",
    ROOT / "p2108_s1058_strict_cmp1_uncertainty_aware_residual_bound.py",
    ROOT / "p2109_s1059_strict_cmp1_covariance_row_link_and_d3_bound.py",
    ROOT / "p2110_s1060_strict_contract_level_cmp_table_d3_rowmap_calibration.py",
    ROOT / "p2111_s1061_strict_cmp2_covariance_row_export_and_residual_execution.py",
    ROOT / "p2112_s1062_strict_cmp2_multibin_covariance_residual_table.py",
    ROOT / "p2113_s1063_strict_cmp2_binwise_d3_covariance_transport_table.py",
    ROOT / "p2114_s1064_strict_cmp2_backend_binwise_covariance_and_rowmap_robustness_audit.py",
    ROOT / "p2115_s1065_strict_cmp2_full_binwise_channel_covariance_slices_and_coupled_robustness.py",
    ROOT / "p2116_s1066_strict_cmp2_full_slice_provenance_and_constructor_sensitivity_audit.py",
    ROOT / "p2117_s1067_strict_cmp2_binalignment_uncertainty_propagated_constructor_audit.py",
    ROOT / "p2118_s1068_strict_cmp2_backend_alignment_covariance_object_audit.py",
    ROOT / "p2119_s1069_strict_cmp2_alignment_covariance_crosscheck_intersection_envelope_audit.py",
    ROOT / "p2120_s1070_strict_cmp2_triple_alignment_estimator_intersection_consistency_audit.py",
    ROOT / "p2121_s1071_strict_cmp2_binlocal_eigenspectrum_triple_intersection_audit.py",
    ROOT / "p2122_s1072_strict_cmp2_index_to_s_provenance_uncertainty_triple_intersection_audit.py",
]
OUT = ROOT / "generated" / "p2122_s1072_strict_cmp2_index_to_s_provenance_uncertainty_triple_intersection_audit.json"


class TestP2122StrictCmp2IndexToSProvenanceUncertaintyTripleIntersectionAudit(unittest.TestCase):
    def test_p2122_exports_index_to_s_mapping_uncertainty_audit(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2122_s1072_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_STRICT_CMP2_INDEX_TO_S_PROVENANCE_UNCERTAINTY_TRIPLE_INTERSECTION_AUDIT_WITH_TRACE",
        )

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["index_to_s_rows_exported"])
        self.assertTrue(checks["mapping_uncertainty_propagated"])
        self.assertFalse(checks["full_d3_covariance_transport_proven"])
        self.assertFalse(checks["full_cutkosky_closure_proven"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
