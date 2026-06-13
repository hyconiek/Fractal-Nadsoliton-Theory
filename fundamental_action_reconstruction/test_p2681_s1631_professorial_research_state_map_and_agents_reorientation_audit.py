from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2681_s1631_professorial_research_state_map_and_agents_reorientation_audit.py"
OUT = ROOT / "generated" / "p2681_s1631_professorial_research_state_map_and_agents_reorientation_audit.json"
MD = ROOT / "generated" / "p2681_s1631_professorial_research_state_map_and_agents_reorientation_audit.md"


class P2681ProfessorialResearchStateMapAndAgentsReorientationAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_broad_content_first_grep_covers_major_lanes(self) -> None:
        audit = self.payload["broad_state_grep"]
        self.assertIn("broad content-first", audit["mode"])
        for key in (
            "legacy_strict_bridge_and_role_transfer",
            "selector_orientation_qw2191",
            "tau_pair12_boundary_square",
            "beta_tors_chi11_and_damping",
            "strict_lagrangian_eom_qg",
            "direct_route_mass_defects_p46",
            "alpha_s_and_couplings",
            "noncyclic_new_provider_constraints",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_prior_ledger_blocks_repetition_lanes(self) -> None:
        ledger = {row["lane"]: row for row in self.payload["prior_lane_ledger"]}
        self.assertFalse(ledger["selector/orientation/QW-2191"]["repeat_now_admissible"])
        self.assertFalse(ledger["tau_src -> pair12 -> boundary-square / S6-O3 typed seed"]["repeat_now_admissible"])
        self.assertFalse(ledger["beta_tors -> chi11"]["repeat_now_admissible"])
        self.assertFalse(ledger["legacy -> strict bridge-source atoms"]["repeat_now_admissible"])
        self.assertTrue(ledger["strict Lagrangian/EOM/QG nonproxy closure"]["repeat_now_admissible"])
        self.assertTrue(ledger["kernel-split-robust direct route P46/N49 mass-defect witnesses"]["repeat_now_admissible"])

    def test_opportunity_matrix_selects_finite_live_frontier_not_generic_bridge(self) -> None:
        matrix = sorted(self.payload["opportunity_matrix"], key=lambda row: row["net_score"], reverse=True)
        self.assertEqual(matrix[0]["candidate"], "finite P46/N49 zero-witness obstruction matrix")
        generic = next(row for row in matrix if row["candidate"] == "generic legacy->strict bridge loop")
        replay = next(row for row in matrix if row["candidate"] == "selector/tau_pair12/beta_tors_chi11 replay")
        self.assertLess(generic["net_score"], matrix[0]["net_score"])
        self.assertLess(replay["net_score"], matrix[0]["net_score"])

    def test_priority_lattice_and_p2680_consistency(self) -> None:
        self.assertTrue(self.payload["p2680_consistency"]["p2680_exists"])
        self.assertTrue(self.payload["p2680_consistency"]["p2680_no_selector_replay"])
        self.assertTrue(self.payload["p2680_consistency"]["p2680_bridge_not_completed"])
        lattice = self.payload["finite_priority_lattice"]
        self.assertEqual(lattice["total_states"], 128)
        self.assertEqual(lattice["passing_states"], 1)
        self.assertEqual(lattice["selected_next_lane"], "finite_p46_n49_zero_witness_matrix")

    def test_agents_docs_updated_and_no_closure(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertTrue(decision["agents_reoriented"])
        self.assertFalse(decision["toe_closed_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2681/S1631", (REPO / "AGENTS.md").read_text(encoding="utf-8"))
        self.assertIn("P2681/S1631", MD.read_text(encoding="utf-8"))
        self.assertIn("P2681/S1631", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2681/S1631", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
