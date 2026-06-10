from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2637_s1587_phase_node_quotient_exhaustion_toe_closure_path_audit.py"
OUT = ROOT / "generated" / "p2637_s1587_phase_node_quotient_exhaustion_toe_closure_path_audit.json"
MD = ROOT / "generated" / "p2637_s1587_phase_node_quotient_exhaustion_toe_closure_path_audit.md"


class P2637PhaseNodeQuotientExhaustionToeClosurePathAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_is_present_and_nonempty(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        self.assertIn("phase_node_zero_lattice_content", audit["patterns"])
        self.assertIn("quotient_reindexing_metric_warp_content", audit["patterns"])

    def test_exact_zero_lattice_and_identity_failure(self) -> None:
        cert = self.payload["phase_node_quotient_exhaustion_certificate"]
        self.assertEqual(cert["exact_zero_lattice_formula"], "d_k = 4/3 + 4*k for k in Z")
        self.assertEqual(cert["declared_legacy_nodes"], ["2", "5", "8", "11"])
        self.assertEqual(cert["audited_zero_lattice"], ["4/3", "16/3", "28/3", "40/3"])
        self.assertFalse(all(row["is_exact_formula_zero"] for row in cert["identity_phase_rows"]))
        self.assertFalse(cert["silent_identity_or_simple_reindexing_passes"])

    def test_affine_metric_pushforward_is_constructive_but_sourced_guarded(self) -> None:
        cert = self.payload["phase_node_quotient_exhaustion_certificate"]
        affine = next(row for row in cert["exhausted_map_classes"] if row["map_class"] == "monotone_affine_metric_pushforward")
        self.assertEqual(affine["definition"], "r(d)=(4/3)*d+(-4/3)")
        self.assertTrue(affine["passes_exact_node_certificate"])
        self.assertEqual(affine["induced_phase_in_legacy_coordinate"], "cos((1/3)*pi*d + (-1/6)*pi)")
        self.assertTrue(cert["constructive_map_requires_new_source_theorem"])

    def test_acceptance_does_not_export_certificate_or_full_kernel(self) -> None:
        acceptance = self.payload["source_acceptance"]
        self.assertTrue(acceptance["gates"]["constructive_affine_pushforward_exists"])
        self.assertFalse(acceptance["gates"]["metric_pushforward_source_theorem_present"])
        self.assertFalse(acceptance["phase_frequency_node_gauge_certificate_exported"])
        decision = self.payload["professorial_toe_closure_path"]["strict_kernel_full_kernel_decision"]
        self.assertTrue(decision["toe_symptoms_present"])
        self.assertTrue(decision["strict_stability_positive"])
        self.assertFalse(decision["is_full_kernel_now"])

    def test_professorial_closure_path_and_negative_exports(self) -> None:
        closure = self.payload["professorial_toe_closure_path"]
        self.assertEqual(len(closure["professorial_closure_path"]), 6)
        self.assertIn("metric-pushforward source", closure["professorial_closure_path"][0]["closure_task"])
        self.assertIn("frozen-kernel", closure["professorial_closure_path"][-1]["computable_exit_condition"])
        self.assertIn("derive or falsify", closure["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))

    def test_docs_updated(self) -> None:
        self.assertIn("P2637/S1587", MD.read_text(encoding="utf-8"))
        self.assertIn("P2637/S1587", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2637/S1587", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
