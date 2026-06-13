from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2673_s1623_tau_src_pair12_boundary_square_subinterface_audit.py"
OUT = ROOT / "generated" / "p2673_s1623_tau_src_pair12_boundary_square_subinterface_audit.json"
MD = ROOT / "generated" / "p2673_s1623_tau_src_pair12_boundary_square_subinterface_audit.md"


class P2673TauSrcPair12BoundarySquareSubinterfaceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_covers_exact_subinterface(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "tau_src_sign_content",
            "pair12_same_packet_content",
            "chart_sensitive_descent_content",
            "boundary_square_cycle_content",
            "sourced_invariant_content",
            "nonclosure_guard_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_upstream_p2672_consistency(self) -> None:
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2672_no_typed_descent_gate"])
        self.assertTrue(upstream["p2672_no_passing_sourced_descent"])
        self.assertTrue(upstream["p2672_no_boundary_phase_bit_target"])
        self.assertTrue(upstream["p2672_no_ltotal_reopening"])

    def test_obligation_vector_has_two_real_inputs_and_three_missing_arrows(self) -> None:
        vector = {item["id"]: item for item in self.payload["obligation_vector"]}
        self.assertTrue(vector["O1_tau_src_sign_exists"]["satisfied_now"])
        self.assertTrue(vector["O2_pair12_carrier_same_tau_packet"]["satisfied_now"])
        self.assertFalse(vector["O3_chart_sensitive_pair12_typed_subinterface"]["satisfied_now"])
        self.assertFalse(vector["O4_boundary_square_cycle_arrow"]["satisfied_now"])
        self.assertFalse(vector["O5_sector_swap_sourced_invariant"]["satisfied_now"])
        self.assertTrue(all(item["content_hits"] > 0 for item in vector.values()))

    def test_finite_closure_lattice_blocks_false_pass(self) -> None:
        lattice = self.payload["finite_closure_lattice"]
        self.assertEqual(lattice["total_states_checked"], 32)
        self.assertEqual(lattice["closure_states_count"], 1)
        self.assertFalse(lattice["current_state_closes"])
        self.assertEqual(lattice["hamming_distance_from_closure"], 3)
        self.assertEqual(
            set(lattice["missing_obligations_now"]),
            {
                "O3_chart_sensitive_pair12_typed_subinterface",
                "O4_boundary_square_cycle_arrow",
                "O5_sector_swap_sourced_invariant",
            },
        )

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertFalse(decision["current_state_closes"])
        self.assertFalse(decision["boundary_phase_bit_target_exported_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["qw2191_discharged_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2673/S1623", MD.read_text(encoding="utf-8"))
        self.assertIn("P2673/S1623", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2673/S1623", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
