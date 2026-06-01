from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2367_s1317_selector_phase_origin_admissibility_no_go_probe.py"
OUT = ROOT / "generated" / "p2367_s1317_selector_phase_origin_admissibility_no_go_probe.json"
MD = ROOT / "generated" / "p2367_s1317_selector_phase_origin_admissibility_no_go_probe.md"

PREREQ_SCRIPTS = [
    ROOT / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.py",
    ROOT / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.py",
    ROOT / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.py",
    ROOT / "p2366_s1316_selector_candidate_phase_reference_chi11_audit_probe.py",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


class P2367SelectorPhaseOriginAdmissibilityNoGoProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        for script in PREREQ_SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT.parent)
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT.parent)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.probe = cls.payload["selector_phase_origin_admissibility_no_go_probe"]

    def test_packet_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["schema_version"], "p2367_s1317_v1")
        self.assertEqual(self.payload["packet_id"], "P2367")
        self.assertEqual(self.payload["stage_id"], "S1317")
        self.assertEqual(
            self.payload["result_kind"],
            "SELECTOR_PHASE_ORIGIN_ADMISSIBILITY_NO_GO_BOUNDARY",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_free_transitive_orbit_certificate(self) -> None:
        rows = self.probe["orbit_certificate"]
        self.assertEqual({row["orientation"] for row in rows}, {-1, 1})
        for row in rows:
            self.assertEqual(row["orbit_size_from_source_rows"], 12)
            self.assertEqual(row["translated_sources_from_source_0"], list(range(12)))
            self.assertTrue(row["translation_action_is_transitive_on_sources"])
            self.assertTrue(row["translation_action_is_free_on_sources"])
            self.assertEqual(row["stabilizer_of_source_0"], [0])

    def test_admissibility_boundary(self) -> None:
        for row in self.probe["admissibility_rows"]:
            self.assertEqual(row["unique_power_signatures"], 1)
            self.assertEqual(row["unique_complete_bispectrum_signatures"], 1)
            self.assertEqual(row["unique_chiral_bispectrum_orientations"], [row["orientation"]])
            self.assertEqual(row["unique_raw_coprime_phase_signatures"], 12)
            self.assertTrue(all(count < 12 for count in row["noncoprime_phase_alias_class_counts"].values()))
            self.assertTrue(all(row["calibrated_coprime_phase_recovers_all_sources_by_mode"].values()))

    def test_candidate_boundary_table_and_limits(self) -> None:
        table = self.probe["candidate_boundary_table"]
        verdicts = {row["candidate_class"]: row["verdict"] for row in table}
        self.assertEqual(
            verdicts["translation_invariant_power_or_histogram_score"],
            "REJECT_AS_SOURCE_LOCALIZER",
        )
        self.assertEqual(
            verdicts["coprime_phase_with_phase_origin_reference"],
            "BEST_OPERATIONAL_CANDIDATE_BUT_PREMISE_BASED",
        )
        theorem = self.probe["theorem_export"]
        self.assertEqual(self.probe["theorem_fingerprint_sha256"], sha256_json(theorem))
        self.assertIn("translation-invariant strict-core score is source-blind", theorem["claim"])
        self.assertIn("QW-2191 discharge", theorem["not_licensed"])
        self.assertIn("internal strict phase-origin", self.payload["recommended_next_honest_step"])

    def test_gatekeepers(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(
            self.payload["global_status"],
            "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
        )


if __name__ == "__main__":
    unittest.main()
