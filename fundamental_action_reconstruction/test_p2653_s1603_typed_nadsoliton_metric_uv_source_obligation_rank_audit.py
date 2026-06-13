from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2653_s1603_typed_nadsoliton_metric_uv_source_obligation_rank_audit.py"
OUT = ROOT / "generated" / "p2653_s1603_typed_nadsoliton_metric_uv_source_obligation_rank_audit.json"
MD = ROOT / "generated" / "p2653_s1603_typed_nadsoliton_metric_uv_source_obligation_rank_audit.md"


class P2653TypedNadsolitonMetricUVSourceObligationRankAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_and_upstream_consistency(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("typed_metric_uv_content", audit["patterns"])
        self.assertIn("scale_orbit_obstruction_content", audit["patterns"])
        self.assertIn("operator_identity_content", audit["patterns"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2649_beta_routes_blocked"])
        self.assertTrue(upstream["p2650_no_canonical_unit_candidates_pass"])
        self.assertTrue(upstream["p2651_beta_one_only_working_gauge"])
        self.assertTrue(upstream["p2652_unit_map_validator_ready_not_source"])
        self.assertTrue(upstream["p2652_real_payload_not_loaded"])

    def test_scale_orbit_indeterminacy_witness(self) -> None:
        orbit = self.payload["scale_orbit_indeterminacy_witness"]
        self.assertEqual(orbit["orbit_size_audited"], 25)
        self.assertTrue(orbit["all_representatives_equivalent_to_roundoff"])
        self.assertLess(orbit["max_error"], 1e-15)

    def test_obligation_matrix_blocks_current_routes(self) -> None:
        rank = self.payload["obligation_rank_audit"]
        self.assertEqual(rank["current_routes_passing"], [])
        self.assertIn("uv_unit_selected_by_nadsoliton_dynamics", rank["currently_missing_atoms_union"])
        self.assertIn("scale_orbit_quotient_discharge", rank["currently_missing_atoms_union"])
        self.assertIn("dimensionless_conservation_or_operator_identity", rank["currently_missing_atoms_union"])
        matrix = {row["route"]: row for row in self.payload["route_atom_matrix"]}
        self.assertFalse(matrix["p2652_precommitted_unit_map_validator"]["passes_typed_metric_uv_source"])
        self.assertTrue(matrix["hypothetical_valid_metric_uv_theorem"]["passes_typed_metric_uv_source"])

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertTrue(decision["scale_orbit_equivalence_verified"])
        self.assertEqual(decision["current_routes_passing_typed_metric_uv_source"], [])
        self.assertFalse(decision["typed_metric_uv_source_theorem_exported_now"])
        self.assertFalse(decision["canonical_unit_exported_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2653/S1603", MD.read_text(encoding="utf-8"))
        self.assertIn("P2653/S1603", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2653/S1603", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
