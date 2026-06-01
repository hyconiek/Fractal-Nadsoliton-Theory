from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.py"
OUT = ROOT / "generated" / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.json"
MD = ROOT / "generated" / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.md"

PREREQ_SCRIPTS = [
    ROOT / "p2073_s1023_strict_full_ltotal_eom_derivation_scaffold_audit.py",
    ROOT / "p2086_s1036_strict_full_ltotal_eom_termwise_execution_audit.py",
    ROOT / "p2087_s1037_strict_full_ltotal_eom_normal_form_extraction_audit.py",
    ROOT / "p2088_s1038_strict_full_ltotal_eom_theorem_readiness_gap_audit.py",
    ROOT / "p2316_s1266_strict_current_best_lagrangian_eom_coverage_audit_probe.py",
    ROOT / "p2329_s1279_selector_independence_lagrangian_eom_audit_probe.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2362StrictLagrangianEomParallelCompletionProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_lagrangian_eom_parallel_completion_probe"]
        cls.summary = cls.probe["parallel_completion_summary"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2362_s1312_v1")
        self.assertEqual(self.payload["packet_id"], "P2362")
        self.assertEqual(self.payload["stage_id"], "S1312")
        self.assertEqual(
            self.payload["result_kind"],
            "STRICT_LAGRANGIAN_EOM_PARALLEL_COMPLETION_AUDIT_SELECTOR_INDEPENDENT_NO_QW2191_DISCHARGE",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_covariant_and_reduced_eom_coverage(self) -> None:
        self.assertEqual(len(self.probe["covariant_sector_eom_rows"]), 5)
        sectors = {row["sector"] for row in self.probe["covariant_sector_eom_rows"]}
        self.assertEqual(sectors, {"scalar_phi", "higgs_H", "gauge_A", "fermion_psi", "metric_g"})
        self.assertTrue(all(row["lagrangian_anchor"] for row in self.probe["covariant_sector_eom_rows"]))
        self.assertTrue(all(row["eom_row"] for row in self.probe["covariant_sector_eom_rows"]))
        self.assertTrue(all(not row["selector_required_for_eom_export"] for row in self.probe["covariant_sector_eom_rows"]))

        self.assertEqual(len(self.probe["reduced_termwise_incidence_rows"]), 11)
        self.assertEqual(self.probe["reduced_termwise_incidence_rank"], 3)
        self.assertEqual(self.probe["reduced_field_column_sums"], {"psi": 5, "A": 5, "h": 4})
        self.assertTrue(all(row["nonzero_fields"] for row in self.probe["reduced_termwise_incidence_rows"]))
        self.assertEqual(len(self.probe["reduced_eom_rows"]), 3)
        for row in self.probe["reduced_eom_rows"]:
            self.assertEqual(row["symbolic_recomposition_residual"], "Integer(0)")
            self.assertEqual(row["numeric_probe_residual"], "0")
            self.assertFalse(row["selector_required_for_reduced_eom"])

    def test_parallel_summary_and_gatekeepers(self) -> None:
        self.assertTrue(self.summary["eom_lagrangian_track_is_selector_independent"])
        self.assertTrue(self.summary["selector_closure_is_parallel_problem"])
        self.assertFalse(self.summary["selector_closure_required_before_eom_work"])
        self.assertEqual(self.summary["covariant_sector_eom_row_count"], 5)
        self.assertEqual(self.summary["reduced_term_count"], 11)
        self.assertEqual(self.summary["reduced_variation_field_count"], 3)
        self.assertEqual(self.summary["reduced_incidence_rank"], 3)
        self.assertEqual(
            self.summary["selector_qw2191_status"],
            "OPEN_PARALLEL_TRACK_NOT_USED_AS_EOM_PREREQUISITE",
        )
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_theorem_fingerprint_and_limits(self) -> None:
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("independently of selector/QW-2191 closure", theorem["claim"])
        self.assertIn("full tensor-resolved metric/background theorem", theorem["not_licensed"])
        self.assertIn("selector premise or QW-2191 selector discharge", theorem["not_licensed"])
        self.assertIn("ToE closure", theorem["not_licensed"])
        self.assertIn("not a prerequisite", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
