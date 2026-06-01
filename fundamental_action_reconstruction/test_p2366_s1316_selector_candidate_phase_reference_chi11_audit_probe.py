from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2366_s1316_selector_candidate_phase_reference_chi11_audit_probe.py"
OUT = ROOT / "generated" / "p2366_s1316_selector_candidate_phase_reference_chi11_audit_probe.json"
MD = ROOT / "generated" / "p2366_s1316_selector_candidate_phase_reference_chi11_audit_probe.md"

PREREQ_SCRIPTS = [
    ROOT / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.py",
    ROOT / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.py",
    ROOT / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2366SelectorCandidatePhaseReferenceChi11AuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["selector_candidate_phase_reference_chi11_audit_probe"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2366_s1316_v1")
        self.assertEqual(self.payload["packet_id"], "P2366")
        self.assertEqual(self.payload["stage_id"], "S1316")
        self.assertEqual(
            self.payload["result_kind"],
            "SELECTOR_CANDIDATE_PHASE_REFERENCE_CHI11_AUDIT",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_candidate_ranking_selects_phase_origin(self) -> None:
        ranking = self.probe["candidate_ranking"]
        self.assertEqual(ranking[0]["rank"], 1)
        self.assertEqual(
            ranking[0]["candidate"],
            "phase_origin_plus_chiral_bispectrum_for_chi11_selector_source",
        )
        self.assertIn("BEST_OPERATIONAL_CANDIDATE", ranking[0]["status"])
        self.assertEqual(self.probe["selected_candidate"], ranking[0])
        self.assertIn("candidate-not-theorem", ranking[1]["status"].lower())

    def test_phase_origin_recovery_scan(self) -> None:
        research = self.probe["phase_origin_selector_research"]
        self.assertEqual(research["row_count"], 24)
        self.assertEqual(research["expected_pair_count"], 24)
        self.assertEqual(research["unique_predicted_orientation_source_pairs_using_mode_1"], 24)
        self.assertTrue(research["all_orientations_recovered_by_chiral_bispectrum"])
        self.assertTrue(research["all_rows_recovered_by_all_coprime_modes"])
        self.assertEqual(set(research["source_recovery_success_by_mode"]), {"1", "5", "7", "11"})
        self.assertTrue(all(research["source_recovery_success_by_mode"].values()))

    def test_negative_controls_and_grep_synthesis(self) -> None:
        controls = self.probe["negative_controls"]
        self.assertTrue(controls["translation_invariant_magnitudes_cannot_select_source"])
        self.assertTrue(controls["noncoprime_modes_fail_full_source_resolution"])
        self.assertTrue(controls["without_chiral_marker_orientation_not_unique"])
        self.assertTrue(all(value == 1 for value in controls["translation_invariant_magnitude_signature_counts_by_orientation"].values()))
        self.assertTrue(all(value < 12 for value in controls["noncoprime_mode_alias_class_counts"].values()))

        synthesis = self.probe["grep_synthesis"]
        self.assertEqual(synthesis["top_logical_bottleneck"], ["chi11_selector_source"])
        self.assertEqual(synthesis["chi11_total_critical_count"], 73)
        self.assertTrue(synthesis["chi11_is_only_singleton_unlock"])
        self.assertEqual(synthesis["strict_full_aut_internal_chi11_polarity_source_antichain"], [])

    def test_bridge_eom_separation_limits_and_fingerprint(self) -> None:
        separation = self.probe["bridge_and_eom_separation"]
        self.assertTrue(separation["selector_work_is_parallel_to_eom_lagrangian"])
        self.assertFalse(separation["selector_required_for_p2365_eom_lift"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("phase-origin source-localizing selector", theorem["claim"])
        self.assertIn("selector/QW-2191 discharge", theorem["not_licensed"])
        self.assertIn("beta_tors -> chi11 theorem", theorem["not_licensed"])
        self.assertIn("phase-origin reference", self.payload["recommended_next_honest_step"])


if __name__ == "__main__":
    unittest.main()
