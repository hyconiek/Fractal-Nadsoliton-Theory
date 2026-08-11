#!/usr/bin/env python3
import json
import unittest
from pathlib import Path

import numpy as np

from fin_st01_st15_research import N, strict_operator
from fin_st28_st45_research import dyadic_lift


ROOT = Path(__file__).resolve().parent


class ST94ST105Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((ROOT / "FIN_ST94_ST105_Results.json").read_text(encoding="utf-8"))
        _, cls.a, _ = strict_operator()

    def test_programs_present(self):
        self.assertTrue(all(f"ST{k}" in self.data for k in range(94,106)))

    def test_st94_unital(self):
        self.assertLess(self.data["ST94"]["maximum_live_residuals"]["unital"],1e-12)

    def test_st94_bimodule(self):
        self.assertLess(self.data["ST94"]["maximum_live_residuals"]["bimodule"],3e-12)

    def test_st94_relative_not_absolute(self):
        self.assertIn("after the observable algebra",self.data["ST94"]["projection_uniqueness_verdict"])

    def test_st95_dimension(self):
        self.assertEqual(self.data["ST95"]["fine_block_parameter_space_dimension"],78)

    def test_st95_interior(self):
        self.assertGreater(self.data["ST95"]["interior_graph_offdiagonal_margin"],0)

    def test_st95_live_fiber(self):
        plus=np.zeros((24,12))
        for x in range(12): plus[x,x]=plus[x+12,x]=1/np.sqrt(2)
        for q in [0,.41,1.1]: self.assertLess(np.linalg.norm(dyadic_lift(self.a,q)@plus-plus@self.a),1e-12)

    def test_st96_pythagorean(self):
        self.assertLess(self.data["ST96"]["maximum_pythagorean_residual"],1e-12)

    def test_st96_pinsker(self):
        self.assertGreater(self.data["ST96"]["minimum_pinsker_margin"],0)

    def test_st97_all_branches(self):
        self.assertEqual(self.data["ST97"]["branches_reached"],list(range(12)))

    def test_st97_symmetric_state_not_selected(self):
        self.assertTrue(self.data["ST97"]["zero_state_remains_zero"])

    def test_st97_qw_open(self):
        self.assertIn("QW2191_remains_open",self.data["ST97"]["status"])

    def test_st98_passivity(self):
        self.assertGreater(self.data["ST98"]["minimum_sampled_real_transfer"],0)

    def test_st99_inclusion(self):
        self.assertIsNotNone(self.data["ST99"]["accepted"])
        self.assertGreater(self.data["ST99"]["accepted"]["minimum_strict_inclusion_margin"],0)

    def test_st99_scope(self):
        self.assertIn("frozen binary64",self.data["ST99"]["operator_scope"])

    def test_st100_pullback(self):
        self.assertEqual(self.data["ST100"]["compatible_pullback_sequences"],[[1,1,1],[-1,1,1]])

    def test_st100_no_stationary_antiperiodic(self):
        self.assertFalse(self.data["ST100"]["stationary_antiperiodic_sequence_exists"])

    def test_st101_shift(self):
        self.assertLess(max(row["fine_scalar_shift_residual"] for row in self.data["ST101"]["rows"]),1e-12)

    def test_st101_recovery(self):
        self.assertLess(max(abs(row["q"]-row["q_recovered_from_probability"]) for row in self.data["ST101"]["rows"]),1e-12)

    def test_st102_commuting_square(self):
        self.assertLess(self.data["ST102"]["maximum_conditional_expectation_commutator_residual"],1e-12)

    def test_st103_common_monotone(self):
        self.assertTrue(self.data["ST103"]["common_monotone"])

    def test_st103_full_goal_open(self):
        self.assertFalse(self.data["ST103"]["independent_two_parameter_interval_certificate_obtained"])

    def test_st104_zero_mode_dimension(self):
        self.assertEqual(self.data["ST104"]["zero_Bohr_frequency_operator_dimension"],22)

    def test_st104_left_open(self):
        self.assertIn("open",self.data["ST104"]["left_inclusion_strictness_for_this_spectrum"])

    def test_st105_coarse_impossibility(self):
        self.assertEqual(self.data["ST105"]["coarse_conditioned_channel_TV"],0)

    def test_st105_fine_protocol(self):
        self.assertEqual(self.data["ST105"]["first_event_count_with_Bayes_error_at_most_0_01"],154)

    def test_recommendations(self):
        self.assertEqual([row["id"] for row in self.data["recommended_next_programs"]],[f"ST{k}" for k in range(106,118)])

    def test_global_boundary(self):
        text=self.data["epistemic_boundary"]
        for token in ["QW-2191","laboratory","Standard Model","ToE"]: self.assertIn(token,text)


if __name__=="__main__":
    unittest.main(verbosity=2)
