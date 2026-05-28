from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


def sha256_json(payload: object) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


class TestP2316S1266StrictCurrentBestLagrangianEomCoverageAuditProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2316_s1266_strict_current_best_lagrangian_eom_coverage_audit_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2316_s1266_strict_current_best_lagrangian_eom_coverage_audit_probe.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(data["schema_version"], "p2316_s1266_v1")
        self.assertEqual(data["packet_id"], "P2316")
        self.assertEqual(
            data["result_kind"],
            "STRICT_REPO_GREP_AND_COMPUTATIONAL_LAGRANGIAN_EOM_COVERAGE_AUDIT_NO_G1_G3_UPDATE",
        )
        probe = data["strict_current_best_lagrangian_eom_coverage_audit_probe"]
        grep_audit = probe["repo_grep_audit"]
        self.assertGreaterEqual(grep_audit["hit_count"], 5)
        self.assertGreaterEqual(len(grep_audit["top_hits"]), 5)
        ranked = probe["candidate_lagrangian_sources_ranked"]
        self.assertGreaterEqual(len(ranked), 10)
        self.assertEqual(ranked[0]["source_id"], "P2086_FULL_LTOTAL_TERMWISE_EOM_AUDIT")
        selected = probe["strongest_current_lagrangian_form"]
        self.assertEqual(selected["selected_working_form_id"], "P2086_FULL_LTOTAL_TERMWISE_EOM_AUDIT")
        self.assertEqual(selected["selected_structural_anchor_id"], "P1653_NONSKELETON_FULL_LAGRANGIAN")
        self.assertIn("L_strict_scalar", selected["canonical_nonskeleton_decomposition"])
        coverage = probe["current_limitations_and_eom_coverage"]
        self.assertEqual(
            coverage["best_current_status"],
            "BEST_WORKING_LTOTAL_IDENTIFIED__FULL_TASK3_THEOREM_STILL_OPEN",
        )
        self.assertEqual(coverage["full_eom_answer"], "NO_FULL_TENSOR_RESOLVED_EOM_FOR_TASK3_YET")
        self.assertEqual(coverage["g1_g3_update"], "NONE__P2282_G1_G3_REMAIN_OPEN_UNCHANGED")
        p2086 = coverage["computational_coverage"]
        self.assertTrue(p2086["loaded"])
        self.assertEqual(p2086["term_count"], 11)
        self.assertEqual(p2086["field_count"], 3)
        self.assertTrue(p2086["all_symbolic_residual_zero"])
        self.assertTrue(p2086["all_numeric_residual_zero"])
        self.assertEqual(p2086["termwise_execution_status"], "TERMWISE_EXECUTED_NON_THEOREM_GRADE")
        self.assertEqual(p2086["theorem_grade_tensor_status"], "OPEN_C3_STILL_OPEN__NON_THEOREM_GRADE")
        gaps = coverage["tensor_and_metric_gaps"]
        self.assertTrue(gaps["all_gap_packets_loaded"])
        self.assertEqual(gaps["tensor_gap_packet_count"], 4)
        theorem = probe["theorem_export"]
        self.assertEqual(sha256_json(theorem), probe["theorem_fingerprint_sha256"])
        self.assertIn("does not export the full tensor-resolved", theorem["claim"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["repo_lagrangian_grep_hits_found"])
        self.assertTrue(g["candidate_sources_loaded"])
        self.assertTrue(g["p1653_nonskeleton_anchor_loaded"])
        self.assertTrue(g["p1693_multisector_bridge_loaded"])
        self.assertTrue(g["p2086_termwise_audit_loaded"])
        self.assertTrue(g["p2086_symbolic_recomposition_zero"])
        self.assertTrue(g["p2086_numeric_probe_zero"])
        self.assertTrue(g["p2086_non_theorem_grade_preserved"])
        self.assertTrue(g["tensor_gap_packets_loaded"])
        self.assertTrue(g["p2314_inventory_loaded"])
        self.assertTrue(g["p2315_schematic_spectrum_loaded"])
        self.assertTrue(g["best_working_ltotal_identified"])
        self.assertTrue(g["full_task3_theorem_not_claimed"])
        self.assertTrue(g["strict_g1_g3_not_updated"])
        self.assertTrue(g["no_legacy_kernel_role_transfer"])
        self.assertTrue(g["no_selector_premise_added"])
        self.assertTrue(g["no_qw2191_discharge_claimed"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
