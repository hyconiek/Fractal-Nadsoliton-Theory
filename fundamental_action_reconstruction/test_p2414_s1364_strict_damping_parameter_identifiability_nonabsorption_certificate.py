#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2414_s1364_strict_damping_parameter_identifiability_nonabsorption_certificate.py"
OUT = ROOT / "generated" / "p2414_s1364_strict_damping_parameter_identifiability_nonabsorption_certificate.json"
MD = ROOT / "generated" / "p2414_s1364_strict_damping_parameter_identifiability_nonabsorption_certificate.md"
P2413 = ROOT / "generated" / "p2413_s1363_amplitude_scalar_normalization_bridge_witness_certificate.json"


class P2414StrictDampingParameterIdentifiabilityNonabsorptionCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2413.exists():
            subprocess.run([sys.executable, str(ROOT / "p2413_s1363_amplitude_scalar_normalization_bridge_witness_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_damping_parameter_identifiability_nonabsorption_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.finite = cls.cert["finite_witness_certificate"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2414")
        self.assertEqual(self.payload["stage_id"], "S1364")
        self.assertEqual(
            self.payload["status"],
            "PASS_STRICT_DAMPING_PARAMETERS_IDENTIFIED_FROM_ACCEPTED_SAMPLES_NO_SOURCE_NO_ABSORPTION",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_beta_eta_identifiability(self) -> None:
        self.assertEqual(self.theorem["domain_size"], 11)
        self.assertEqual(self.theorem["accepted_strict_denominator_samples"], "S(d)=1+d^(9/5) on d=1..11")
        self.assertTrue(self.theorem["candidate_grid_unique_eta_match"])
        self.assertEqual(self.theorem["matching_rational_candidates"], ["9/5"])
        self.assertTrue(self.theorem["strict_beta_eta_identified_from_samples"])
        self.assertGreater(self.theorem["candidate_grid_size"], 100)

    def test_linear_nonabsorption_rows_and_proof(self) -> None:
        rows = self.finite["finite_rows"]
        self.assertEqual([row["d"] for row in rows], list(range(1, 12)))
        self.assertTrue(self.theorem["required_gamma_values_strictly_increase"])
        self.assertTrue(self.theorem["no_single_linear_gamma_matches_two_distinct_positive_nodes"])
        self.assertTrue(self.theorem["legacy_beta_tors_matches_no_positive_strict_node"])
        self.assertTrue(self.theorem["legacy_denominator_residuals_all_positive"])
        self.assertIn("gamma=d^(4/5)", self.finite["nonabsorption_proof"]["linear_gamma_no_go"])
        self.assertIn("(S(2)-1)^5 = 2^9", rows[1]["fifth_power_cover_identity"])

    def test_hard_limits_and_inheritance(self) -> None:
        self.assertFalse(self.theorem["strict_beta_eta_source_exported"])
        self.assertFalse(self.theorem["legacy_beta_tors_to_beta_eta_theorem_exported"])
        self.assertFalse(self.theorem["damping_compression_bridge_component_ready"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertTrue(self.theorem["p2411_full_bridge_still_requires_all_obligations"])
        self.assertTrue(self.theorem["p2413_amplitude_witness_inherited"])
        self.assertTrue(self.theorem["scratch_damping_separation_inherited"])
        self.assertTrue(self.theorem["scratch_damping_identifiability_inherited"])
        self.assertIn("No strict dynamic source", "\n".join(self.theorem["not_licensed"]))

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("damping parameter identifiability", MD.read_text(encoding="utf-8"))
        self.assertIn("P2414/S1364 strict damping", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2414/S1364 strict damping", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
