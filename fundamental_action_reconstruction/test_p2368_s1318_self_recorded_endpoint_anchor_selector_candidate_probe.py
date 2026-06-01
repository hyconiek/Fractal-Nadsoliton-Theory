from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2368_s1318_self_recorded_endpoint_anchor_selector_candidate_probe.py"
OUT = ROOT / "generated" / "p2368_s1318_self_recorded_endpoint_anchor_selector_candidate_probe.json"
MD = ROOT / "generated" / "p2368_s1318_self_recorded_endpoint_anchor_selector_candidate_probe.md"

PREREQ_SCRIPTS = [
    ROOT / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.py",
    ROOT / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.py",
    ROOT / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.py",
    ROOT / "p2366_s1316_selector_candidate_phase_reference_chi11_audit_probe.py",
    ROOT / "p2367_s1317_selector_phase_origin_admissibility_no_go_probe.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2368SelfRecordedEndpointAnchorSelectorCandidateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["self_recorded_endpoint_anchor_selector_candidate_probe"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2368_s1318_v1")
        self.assertEqual(self.payload["packet_id"], "P2368")
        self.assertEqual(self.payload["stage_id"], "S1318")
        self.assertEqual(
            self.payload["result_kind"],
            "SELF_RECORDED_ENDPOINT_ANCHOR_SELECTOR_CANDIDATE",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_ledger_action_certificate(self) -> None:
        ledger = self.probe["ledger_action_certificate"]
        self.assertEqual(ledger["composition_count"], 35)
        self.assertTrue(ledger["unique_lexicographic_winner"])
        self.assertEqual(ledger["winners"], [[2, 2, 2, 1, 1]])
        self.assertEqual(ledger["minimum_pair"], [14, 0])
        self.assertGreater(ledger["equal_ripple_competitor_count"], 0)
        self.assertTrue(ledger["all_equal_ripple_competitors_have_positive_arrow_penalty"])

    def test_anchor_reconstruction_and_equivariance(self) -> None:
        anchor = self.probe["anchor_reconstruction_certificate"]
        self.assertEqual(anchor["row_count"], 24)
        self.assertTrue(anchor["all_sources_recovered"])
        self.assertTrue(anchor["all_orientations_recovered"])
        self.assertTrue(anchor["all_ordered_values_match_balanced_ledger"])

        equivariance = self.probe["d12_equivariance_certificate"]
        self.assertEqual(equivariance["checked_cases"], 576)
        self.assertEqual(equivariance["mismatch_count"], 0)
        self.assertTrue(equivariance["all_cases_equivariant"])

    def test_negative_controls_and_candidate_boundary(self) -> None:
        controls = self.probe["negative_controls"]
        self.assertEqual(controls["value_multiset_signature_count"], 1)
        self.assertTrue(controls["value_multiset_is_source_and_orientation_blind"])
        self.assertTrue(controls["reversed_ledger_moves_source_to_opposite_endpoint_and_flips_orientation"])
        self.assertTrue(controls["absolute_origin_not_selected_by_equivariant_anchor"])

        boundary = self.probe["candidate_boundary"]
        self.assertIn("external Fourier phase reference", boundary["stronger_than_p2366_phase_origin_candidate_because"])
        self.assertIn("does not contradict P2367", boundary["relation_to_p2367_no_go"])

    def test_gatekeepers_fingerprint_and_limits(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("D12-equivariant", theorem["claim"])
        self.assertIn("QW-2191 discharge", theorem["not_licensed"])
        self.assertIn("self-recorded arrow action", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
