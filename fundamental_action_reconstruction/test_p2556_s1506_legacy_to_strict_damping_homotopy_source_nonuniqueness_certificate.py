from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2556_s1506_legacy_to_strict_damping_homotopy_source_nonuniqueness_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2556DampingHomotopySourceNonuniquenessTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["legacy_to_strict_damping_homotopy_source_nonuniqueness_certificate"]["theorem_export"]

    def test_identity_and_nonuniqueness(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2556")
        self.assertEqual(self.payload["stage_id"], "S1506")
        self.assertTrue(self.theorem["p2555_denominator_nonrenormalization_inherited"])
        self.assertEqual(self.theorem["row_count"], 22)
        self.assertTrue(self.theorem["both_homotopies_have_same_endpoints_and_transport_primitive"])
        self.assertTrue(self.theorem["instantaneous_sources_differ_for_all_rows"])
        self.assertTrue(self.theorem["homotopy_endpoint_data_do_not_select_unique_source_dynamics"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["unique_damping_homotopy_source_exported"])
        self.assertFalse(self.theorem["damping_compression_bridge_component_ready"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertIn("homotopy/source density", self.theorem["recommended_next_honest_step"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2556/S1506", MD.read_text(encoding="utf-8"))
        self.assertIn("P2556/S1506", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2556/S1506", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
