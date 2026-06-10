from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2638_s1588_metric_pushforward_source_viability_and_neural_closure_audit.py"
OUT = ROOT / "generated" / "p2638_s1588_metric_pushforward_source_viability_and_neural_closure_audit.json"
MD = ROOT / "generated" / "p2638_s1588_metric_pushforward_source_viability_and_neural_closure_audit.md"


class P2638MetricPushforwardSourceViabilityAndNeuralClosureAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present_and_nonempty(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        self.assertIn("metric_pushforward_source_content", audit["patterns"])
        self.assertIn("neural_universe_empirical_closure_content", audit["patterns"])

    def test_phase_nodes_are_repaired_exactly(self) -> None:
        repair = self.payload["phase_node_repair_check"]
        self.assertEqual(repair["candidate_pushforward"], "r(d)=(4/3)*d+(-4/3)")
        self.assertTrue(repair["all_legacy_nodes_become_strict_phase_zeros"])
        self.assertTrue(all(row["is_exact_zero_with_tolerance_1e_12"] for row in repair["node_rows"]))

    def test_domain_admissibility_fails_as_global_uv_metric(self) -> None:
        domain = self.payload["domain_admissibility_check"]
        self.assertFalse(domain["is_real_on_open_uv_interval_0_to_1"])
        self.assertTrue(domain["collapses_d_equals_1_to_strict_zero_distance"])
        self.assertTrue(domain["requires_domain_cut_or_shifted_uv_origin"])
        self.assertIn("fails_as_global_positive_distance_metric", domain["admissibility_verdict"])

    def test_damping_distortion_does_not_preserve_inverse_hierarchy(self) -> None:
        distortion = self.payload["damping_and_attention_distortion"]
        ratios = distortion["ratio_abs_k7_over_abs_k1"]
        self.assertGreater(ratios["legacy_amplitude_normalized"], 1.0)
        self.assertLess(ratios["strict_after_metric_pushforward"], 1.0)
        self.assertFalse(ratios["pushforward_preserves_legacy_inverse_hierarchy_above_one"])
        self.assertEqual(distortion["normalized_attention_proxy_after_pushforward"]["argmax_d"], 1)

    def test_closure_decision_keeps_kernel_nonfull_and_nonpromoting(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertTrue(decision["gates"]["phase_nodes_repaired_exactly"])
        self.assertFalse(decision["gates"]["global_positive_distance_metric_on_uv_domain"])
        self.assertFalse(decision["gates"]["inverse_hierarchy_role_preserved_after_pushforward"])
        self.assertFalse(decision["metric_pushforward_promoted_to_bridge_theorem"])
        self.assertFalse(decision["full_kernel_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_professorial_path_and_docs_updated(self) -> None:
        path = self.payload["professorial_neural_toe_path"]
        self.assertIn("strict remains", path["full_kernel_answer"])
        self.assertIn("topology/selector certificates", path["next_honest_step"])
        self.assertEqual(len(path["strongest_toe_symptoms_now"]), 3)
        self.assertIn("P2638/S1588", MD.read_text(encoding="utf-8"))
        self.assertIn("P2638/S1588", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2638/S1588", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
