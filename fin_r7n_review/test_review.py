"""Adversarial controls for the new intake checker; no source-archive writes."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import unittest
from fin_r7n_review import review as r
from fin_r7n_review import finalize as final
from fin_r7n_review import partial_cover as partial

class ReviewTests(unittest.TestCase):
    def test_complete_small_partition(self):
        root=((F(0),F(1)),)*3
        children=[tuple((F(i,2),F(i+1,2)) for i in idx) for idx in product(range(2),repeat=3)]
        r.partition(root,children)

    def test_missing_cell_rejected(self):
        root=((F(0),F(1)),)*3
        children=[tuple((F(i,2),F(i+1,2)) for i in idx) for idx in product(range(2),repeat=3)]
        with self.assertRaises(AssertionError):r.partition(root,children[:-1])

    def test_duplicate_cell_rejected(self):
        root=((F(0),F(1)),)*3
        with self.assertRaises(AssertionError):r.partition(root,[root,root])

    def test_incomplete_phase_replay_rejected(self):
        records={k:final.load('leaf_replay_'+k+'.json') for k in ['quartic','16','20']}
        records['20']['complete']=False
        with self.assertRaises(AssertionError):final.phase_gate(records,final.load('phase_geometry.json'),final.load('collar_replay.json'))

    def test_nonunique_collar_rejected(self):
        records={k:final.load('leaf_replay_'+k+'.json') for k in ['quartic','16','20']}
        collars=deepcopy(final.load('collar_replay.json'));collars['full']['roots'][0]['q_upper']='1'
        with self.assertRaises(AssertionError):final.phase_gate(records,final.load('phase_geometry.json'),collars)

    def test_quartic_formula_agrees_with_direct_field(self):
        phi=[r.sr.I(F(x,10)) for x in [4,7,11]]
        H=r.quartet_hessian(phi);direct=r.sr.phase_FH(phi,'quartic')[1]
        for i in range(3):
            for j in range(3):
                a,b=r.sr.bounds(H[i][j]);c,d=r.sr.bounds(direct[i][j])
                self.assertLessEqual(max(a,c),min(b,d))

    def test_rank_deficient_compression_rejected(self):
        row={'basis_den':1,'basis_num':[[0,0,0] for _ in range(4)]}
        with self.assertRaises(AssertionError):partial.compression(row,[(F(1),F(1))]*7,[[r.sr.I(0)]*4 for _ in range(7)])

    def test_probability_product_bounds(self):
        w=[(F(1),F(2))]*7
        pair=partial.pair_upper(w,0,1);triple=partial.triple_upper(w,(0,1,2))
        for a in product([F(1),F(3,2),F(2)],repeat=3):
            denominator=sum(a)+4
            self.assertGreaterEqual(pair,a[0]*a[1]/denominator**2)
            self.assertGreaterEqual(triple,a[0]*a[1]*a[2]/denominator**3)

if __name__=='__main__':unittest.main()
