from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
P1848 = ROOT / "p1848_s798_strict_gravity_componentwise_variation_and_counterterm_witness_checkpoint.py"
P2030 = ROOT / "p2030_s980_strict_tensor_resolved_projection_source_audit.py"
P2031 = ROOT / "p2031_s981_strict_b1_tensor_component_profile_table_scaffold.py"
OUT2031 = ROOT / "generated" / "p2031_s981_strict_b1_tensor_component_profile_table_scaffold.json"


class TestStrictB1TensorComponentProfileTableScaffold(unittest.TestCase):
    def test_p2031_exports_4x4_scaffold_without_false_tensor_fill(self) -> None:
        subprocess.run([sys.executable, str(P1848)], check=True)
        subprocess.run([sys.executable, str(P2030)], check=True)
        subprocess.run([sys.executable, str(P2031)], check=True)
        data = json.loads(OUT2031.read_text(encoding="utf-8"))

        self.assertEqual(data["schema_version"], "p2031_s981_v1")
        self.assertEqual(data["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "OPEN_TENSOR_COMPONENT_PROFILE_TABLE_SCAFFOLD_EXPORTED_WITH_MISSING_ENTRIES",
        )
        self.assertEqual(data["table_object"]["object_id"], "strict_B1_tensor_component_profile_table_v1")

        entries = data["table_object"]["entries"]
        self.assertEqual(len(entries), 16)
        self.assertEqual(data["table_summary"]["total_required_entries"], 16)
        self.assertEqual(data["table_summary"]["derived_entry_count"], 0)
        self.assertEqual(data["table_summary"]["missing_entry_count"], 16)
        self.assertEqual(data["table_summary"]["conditional_gb_identity_row_count"], 4)
        self.assertFalse(data["table_summary"]["tensor_component_profile_table_ready"])

        self.assertTrue(all(row["profile_status"] == "OPEN_MISSING_TENSOR_COMPONENT_PROFILE" for row in entries))
        self.assertTrue(all(row["fill_status"] == "UNFILLED_BY_DESIGN" for row in entries))
        self.assertTrue(all(row["scalar_shadow"]["available"] for row in entries))
        self.assertTrue(all(row["covariant_template_available"] for row in entries))

        gb_entries = [row for row in entries if row["channel"] == "GB"]
        self.assertEqual(len(gb_entries), 4)
        self.assertTrue(
            all(
                row["gb_identity_condition"]["status"] == "CONDITIONAL_ON_NON_GB_COMPONENT_PROFILES"
                for row in gb_entries
            )
        )
        self.assertTrue(all(row["gb_identity_condition"]["can_fill_now"] is False for row in gb_entries))

        gb_00 = next(row for row in gb_entries if row["component"] == "00")
        self.assertEqual(gb_00["adm_lapse_relation"]["status"], "RELATED_REDUCED_LAPSE_EQUATION_NOT_H_00_PROFILE")
        self.assertFalse(gb_00["adm_lapse_relation"]["can_fill_from_adm_lapse"])
        self.assertTrue(data["adm_lapse_evidence_by_channel"]["GB"]["gb_lapse_cancels"])

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["scaffold_exported"])
        self.assertTrue(checks["entry_count_4x4"])
        self.assertTrue(checks["no_component_profile_marked_derived"])
        self.assertTrue(checks["all_entries_marked_missing"])
        self.assertTrue(checks["conditional_gb_identity_rows_only"])
        self.assertTrue(checks["scalar_shadows_available_but_not_promoted"])
        self.assertTrue(checks["adm_lapse_available_but_not_promoted"])
        self.assertTrue(checks["gb_lapse_cancellation_not_promoted_to_H00_zero"])
        self.assertTrue(checks["no_tensor_projection_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])

        self.assertIn("tensor-component Gram rule", data["blocking_requirements"])
        self.assertIn("using ADM lapse cancellation as H_00 cancellation", data["theorem_scope"]["not_licensed"])
        self.assertEqual(data["component_gram_rule_stub"]["status"], "OPEN_MISSING_RULE")
        self.assertEqual(data["divergence_tensor_target_stub"]["status"], "OPEN_MISSING_TARGET")


if __name__ == "__main__":
    unittest.main()
