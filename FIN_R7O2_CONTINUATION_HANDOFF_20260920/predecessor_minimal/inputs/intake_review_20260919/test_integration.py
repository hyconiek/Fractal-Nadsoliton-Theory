"""Negative controls for the new consolidation checks; no archive writes."""
from copy import deepcopy
from fractions import Fraction as F
import unittest
from fin_rank7_intake_review import finalize as audit
from fin_rank7_intake_review import scientific_rechecks as science

class IntegrationTests(unittest.TestCase):
    def test_all_repair_trees_complete(self):
        self.assertEqual(sum(audit.validate_repair(r) for r in audit.load('FR223_subdivision_repairs.json')),18)

    def test_removed_leaf_rejected(self):
        r=deepcopy(audit.load('FR223_subdivision_repairs.json')[0]);r['leaves'].pop()
        with self.assertRaises(AssertionError):audit.validate_repair(r)

    def test_false_sign_flag_rejected(self):
        r=deepcopy(audit.load('FR223_subdivision_repairs.json')[0])
        r['leaves'][0]['bounds']['P']=['-1','1'];r['leaves'][0]['bounds']['P1']=['-2','-1']
        with self.assertRaises(AssertionError):audit.validate_repair(r)

    def test_false_coverage_box_rejected(self):
        r=deepcopy(audit.load('FR223_subdivision_repairs.json')[0])
        r['leaves'][0]['box'][0][0]='0'
        with self.assertRaises(AssertionError):audit.validate_repair(r)

    def test_true_phase_symmetry_gradient(self):
        for kind in ['quartic','full']:
            g,_=science.phase_FH([science.I(0)]*3,kind)
            self.assertTrue(all(science.bounds(v)[0]<=0<=science.bounds(v)[1] for v in g))

    def test_bad_phase_locator_rejected(self):
        row=audit.load('phase_recertification.json')['full']['roots'][0]
        phase=[float(F(x)) for x in row['center']];phase[0]+=.03
        with self.assertRaises(AssertionError):
            science.phase_root({'phase':phase,'negative_index':row['negative_index']},'full')

    def test_navigation_buffer_not_a_proof(self):
        data=audit.load('FR223_union_replay.json')
        self.assertIs(data['navigation_buffer_used'],False)
        self.assertEqual(data['count'],106)
        self.assertEqual(data['passed'],99)

if __name__=='__main__':unittest.main()
