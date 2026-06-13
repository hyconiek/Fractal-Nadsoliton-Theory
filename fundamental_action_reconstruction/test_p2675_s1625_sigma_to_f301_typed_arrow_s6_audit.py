from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2675_s1625_sigma_to_f301_typed_arrow_s6_audit.py"
OUT = ROOT / "generated" / "p2675_s1625_sigma_to_f301_typed_arrow_s6_audit.json"
MD = ROOT / "generated" / "p2675_s1625_sigma_to_f301_typed_arrow_s6_audit.md"


class P2675SigmaToF301TypedArrowS6AuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_covers_s6_not_numbers_only(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "sigma_codomain_content",
            "f301_carrier_content",
            "typed_arrow_target_content",
            "precollapse_constraint_content",
            "projector_collapse_constraint_content",
            "fiat_blocker_content",
            "actual_export_blocker_content",
            "closure_guard_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_upstream_p2674_consistency(self) -> None:
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2674_names_s6_missing"])
        self.assertTrue(upstream["p2674_o3_not_passed"])
        self.assertTrue(upstream["p2674_blocks_o4"])
        self.assertTrue(upstream["p2674_blocks_ltotal"])

    def test_s6_acceptance_vector_separates_inputs_from_arrow(self) -> None:
        vector = {item["id"]: item for item in self.payload["s6_acceptance_vector"]}
        self.assertTrue(vector["A1_sigma_sel_src_target_present"]["satisfied_now"])
        self.assertTrue(vector["A2_surviving_f301_pair12_carrier_present"]["satisfied_now"])
        self.assertTrue(vector["A3_same_tau_packet_link_present"]["satisfied_now"])
        self.assertFalse(vector["A4_chart_label_retaining_arrow_exported"]["satisfied_now"])
        self.assertFalse(vector["A5_pre_collapse_nonquotient_descent_exported"]["satisfied_now"])
        self.assertFalse(vector["A6_nonprojector_local_atlas_descent_exported"]["satisfied_now"])
        self.assertFalse(vector["A7_no_fiat_identification_proof_exported"]["satisfied_now"])
        self.assertTrue(all(item["content_hits"] > 0 for item in vector.values()))

    def test_finite_s6_lattice_blocks_false_pass(self) -> None:
        lattice = self.payload["finite_s6_lattice"]
        self.assertEqual(lattice["total_states_checked"], 128)
        self.assertEqual(lattice["passing_s6_states_count"], 1)
        self.assertFalse(lattice["current_state_passes_s6"])
        self.assertEqual(lattice["hamming_distance_from_s6_pass"], 4)
        self.assertEqual(
            set(lattice["missing_s6_obligations_now"]),
            {
                "A4_chart_label_retaining_arrow_exported",
                "A5_pre_collapse_nonquotient_descent_exported",
                "A6_nonprojector_local_atlas_descent_exported",
                "A7_no_fiat_identification_proof_exported",
            },
        )

    def test_obstruction_table_and_no_closure_docs_updated(self) -> None:
        lanes = {row["lane"] for row in self.payload["obstruction_table"]}
        self.assertIn("basis_free_Q_basis_continuation", lanes)
        self.assertIn("local_pair12_atlas_projector_lane", lanes)
        self.assertIn("route_local_T220_T222_seed_target_family", lanes)
        self.assertIn("declaration_or_identification_shortcut", lanes)
        decision = self.payload["closure_decision"]
        self.assertFalse(decision["s6_exported_now"])
        self.assertFalse(decision["o3_exported_now"])
        self.assertFalse(decision["boundary_square_arrow_allowed_next"])
        self.assertFalse(decision["sector_swap_invariant_allowed_next"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2675/S1625", MD.read_text(encoding="utf-8"))
        self.assertIn("P2675/S1625", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2675/S1625", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
