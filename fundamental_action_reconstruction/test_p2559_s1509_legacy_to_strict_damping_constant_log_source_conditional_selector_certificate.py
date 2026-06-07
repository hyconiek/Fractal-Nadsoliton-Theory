from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2559_s1509_legacy_to_strict_damping_constant_log_source_conditional_selector_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2559ConstantLogSourceConditionalSelectorTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["legacy_to_strict_damping_constant_log_source_conditional_selector_certificate"]["theorem_export"]

    def test_constant_log_source_selects_geometric(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2559")
        self.assertEqual(self.payload["stage_id"], "S1509")
        self.assertTrue(self.theorem["p2558_power_mean_continuum_inherited"])
        self.assertEqual(self.theorem["row_count"], 22)
        self.assertEqual(self.theorem["power_mean_q_values_audited"], [-2.0, -1.0, 0.0, 1.0, 2.0])
        self.assertEqual(self.theorem["source_density_samples"], [0.0, 0.5, 1.0])
        self.assertTrue(self.theorem["constant_log_source_selects_geometric_q0_for_all_rows"])
        self.assertTrue(self.theorem["conditional_geometric_selector_exported"])
        self.assertTrue(self.theorem["constant_log_source_premise_is_unsourced_bridge_obligation"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["constant_log_source_law_exported"])
        self.assertFalse(self.theorem["geometric_homotopy_source_exported"])
        self.assertFalse(self.theorem["unique_damping_homotopy_source_exported"])
        self.assertFalse(self.theorem["damping_compression_bridge_component_ready"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertIn("constant in log-denominator time", self.theorem["recommended_next_honest_step"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2559/S1509", MD.read_text(encoding="utf-8"))
        self.assertIn("P2559/S1509", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2559/S1509", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
