from __future__ import annotations
import hashlib
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


def sha256_json(payload: object) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


class TestP2298S1248StrictNonGbProviderSourceNonadmissibilityTheoremProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2298_s1248_strict_non_gb_provider_source_nonadmissibility_theorem_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2298_s1248_strict_non_gb_provider_source_nonadmissibility_theorem_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2298_s1248_v1")
        probe = data["strict_non_gb_provider_source_nonadmissibility_theorem_probe"]
        result = probe["nonadmissibility_result"]
        self.assertEqual(result["legal_strict_capable_sources"], [])
        self.assertEqual(
            result["theorem_grade_current_inventory_verdict"],
            "NO_CURRENT_LEGAL_STRICT_PROVIDER_SOURCE_FOR_FULL_P2297_RESIDUAL_FAMILY",
        )
        self.assertIn("formal_full_residual_basis_provider_from_P2297", result["capable_but_non_admissible_sources"])
        self.assertTrue(probe["admissibility_matrix"]["strict_legal_polynomial_is_zero"])
        self.assertEqual(sha256_json(probe["theorem_export"]), probe["theorem_fingerprint_sha256"])
        g = data["gatekeeper_checks"]
        self.assertTrue(g["p2297_provider_matrix_loaded"])
        self.assertTrue(g["full_residual_family_detected"])
        self.assertTrue(g["formal_basis_rank_positive"])
        self.assertTrue(g["strict_core_matrix_inconsistent"])
        self.assertTrue(g["augmented_non_strict_matrix_inconsistent"])
        self.assertTrue(g["formal_full_basis_algebraically_capable"])
        self.assertTrue(g["formal_full_basis_marked_non_admissible"])
        self.assertTrue(g["no_current_legal_strict_capable_provider_source"])
        self.assertTrue(g["all_capable_sources_require_new_source_or_non_strict_selector"])
        self.assertTrue(g["qw2191_replay_certified"])
        self.assertTrue(g["qw2191_not_discharged"])
        self.assertTrue(g["no_selector_closure_claimed"])
        self.assertTrue(g["no_legacy_transfer_used"])
        self.assertTrue(g["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
