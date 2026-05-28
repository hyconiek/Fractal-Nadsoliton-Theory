from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2331_s1281_strict_second_background_curvature_witness_gap_audit.py"
OUT = ROOT / "generated" / "p2331_s1281_strict_second_background_curvature_witness_gap_audit.json"
MD = ROOT / "generated" / "p2331_s1281_strict_second_background_curvature_witness_gap_audit.md"


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2331SecondBackgroundCurvatureWitnessGapAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["strict_second_background_curvature_witness_gap_audit"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2331_s1281_v1")
        self.assertEqual(self.payload["packet_id"], "P2331")
        self.assertEqual(self.payload["stage_id"], "S1281")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_capability_matrix_blocks_global_lift(self) -> None:
        matrix = self.probe["lift_capability_matrix"]
        self.assertFalse(matrix["all_fields_met_by_current_union"])
        self.assertEqual(matrix["admissible_single_source_rows"], [])
        self.assertLess(matrix["covered_field_count"], matrix["required_field_count"])
        self.assertEqual(self.probe["p2031_tensor_table_summary"]["missing_entry_count"], 16)
        self.assertEqual(self.probe["p2031_tensor_table_summary"]["total_required_entries"], 16)

    def test_existing_witnesses_preserved_but_not_promoted(self) -> None:
        rows = {row["source_id"]: row for row in self.probe["second_background_candidate_rows"]}
        self.assertTrue(rows["P1985_ADM_BIANCHI_NON_GB_LAPSE"]["independent_background_family"])
        self.assertFalse(rows["P1985_ADM_BIANCHI_NON_GB_LAPSE"]["gb_channel_in_same_basis"])
        self.assertTrue(rows["P2297_NON_GB_SPATIAL_EOM_PROVIDER"]["independent_background_family"])
        self.assertFalse(rows["P2297_NON_GB_SPATIAL_EOM_PROVIDER"]["transport_or_globalization_theorem"])
        self.assertEqual(self.probe["p2330_gb_dependence_summary"]["numeric_rank"], 3)
        self.assertEqual(self.probe["p2330_gb_dependence_summary"]["numeric_nullity"], 1)

    def test_gatekeepers_and_fingerprint(self) -> None:
        checks = self.payload["gatekeeper_checks"]
        for key in [
            "grep_hits_found",
            "p1985_loaded",
            "p2297_loaded",
            "p2030_tensor_projection_not_ready",
            "p2031_all_tensor_entries_missing",
            "p2033_ansatz_nonavailable",
            "p2330_rank3_nullity1_preserved",
            "no_admissible_second_background_lift_source",
            "current_union_does_not_cover_lift_fields",
            "current_exports_only_quotient_scope",
            "no_full_global_renormalization_claimed",
            "no_qw2191_discharge_claimed",
            "no_g1_g3_update_claimed",
            "no_toe_closure_claimed",
        ]:
            self.assertTrue(checks[key], key)
        theorem = self.probe["current_backend_scope_theorem"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertTrue(theorem["formal_current_export_condition_evaluates_true"])


if __name__ == "__main__":
    unittest.main()
