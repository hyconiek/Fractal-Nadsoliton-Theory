#!/usr/bin/env python3
"""Acceptance tests for FIN P527--P536 generated artifacts."""

from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parent


class Programs527536Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = json.loads((ROOT / "FIN_Programs_527_536_Results.json").read_text(encoding="utf-8"))
        cls.results = cls.payload["results"]

    def test_p527_minimal_mechanism(self) -> None:
        row = self.results["P527"]
        self.assertIn("-g**2*rho**2/(2*b)", row["effective_density_energy"])
        self.assertFalse(row["magnitude_identified"])
        self.assertEqual(len(row["axiom_removal_ledger"]), 4)

    def test_p528_replay(self) -> None:
        row = self.results["P528"]
        self.assertTrue(row["all_acceptance_inequalities_pass"])
        self.assertEqual(row["chart_count"] + row["bridge_count"], 415)

    def test_p529_dynamical_falsification(self) -> None:
        row = self.results["P529"]
        self.assertTrue(row["stable_sector_remains_orbitally_bounded"])
        self.assertTrue(row["negative_VK_case_exhibits_growth"])
        for case in row["cases"]:
            self.assertLess(case["runs"][-1]["power_relative_drift"], 1e-8)

    def test_p530_realization_dependence(self) -> None:
        row = self.results["P530"]
        self.assertEqual(row["relaxation_stable_speed_count"], 0)
        self.assertGreater(row["hamiltonian_stable_speed_count"], 0)

    def test_p531_no_candidate(self) -> None:
        row = self.results["P531"]
        self.assertEqual(row["linear_localization_threshold"], 0.125)
        self.assertEqual(row["localized_relative_periodic_candidates"], [])

    def test_p532_recurrence(self) -> None:
        row = self.results["P532"]
        self.assertLessEqual(row["best_integer_time"], row["integer_time_search_limit"])
        self.assertLess(row["best_phase_sup_error"], 0.15)

    def test_p533_certified_bound(self) -> None:
        row = self.results["P533"]
        self.assertGreater(row["minimum_response_denominator"], 0)
        self.assertLess(row["maximum_combined_bound"], 0.07)

    def test_p534_independent_agreement(self) -> None:
        row = self.results["P534"]
        self.assertEqual(row["exact_class_agreement_count"], row["terminal_box_count"])
        self.assertEqual(row["contradiction_indices"], [])

    def test_p535_nonempty(self) -> None:
        row = self.results["P535"]
        self.assertTrue(row["optimizer_success"])
        self.assertGreater(row["minimum_recommended_pairwise_TV_lower"], 0)

    def test_p536_prior_dependence(self) -> None:
        row = self.results["P536"]
        self.assertTrue(row["selection_varies_across_declared_priors"])

    def test_standard_library_checkers(self) -> None:
        completed = subprocess.run(
            ["python3", str(ROOT / "check_fin_p528_p533_p534_certificates.py")],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        replay = json.loads(completed.stdout)
        self.assertTrue(all(block["accepted"] for block in replay.values()))


if __name__ == "__main__":
    unittest.main(verbosity=2)
