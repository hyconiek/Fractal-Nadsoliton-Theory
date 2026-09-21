"""New intake controls. The six inherited interval-backend tests run separately."""
import copy
from fractions import Fraction as F
import unittest

from fin_r7o3_review import review as r
from fin_r7o3_review.finalize import proof_gate, completion_gate


class ReviewTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certs=r.lines(r.SOURCE/'certificates/active_leaf_certificates.jsonl')
        cls.pilot=r.load(r.HERE/'pilot_replay.json')['certificates']
        cls.ck=r.backend(True)

    def fixture(self):
        parent=((F(0),F(1)),)*4
        L,R=list(parent),list(parent)
        L[0]=(F(0),F(1,2));R[0]=(F(1,2),F(1))
        node={'kind':'SPLIT','cell':parent,'axis':0,'split':'1/2',
              'left':{'kind':'SAFE','cell':L,'leaf_id':0},
              'right':{'kind':'SAFE','cell':R,'leaf_id':1}}
        leaves={0:{'cell':L,'original_index':5},1:{'cell':R,'original_index':5}}
        return parent,node,leaves

    def test_geometry_complete(self):
        parent,node,leaves=self.fixture();used=set()
        r.tree(node,parent,leaves,5,used)
        self.assertEqual(used,{0,1})

    def test_deleted_leaf(self):
        parent,node,leaves=self.fixture();del leaves[0]
        with self.assertRaises(KeyError):r.tree(node,parent,leaves,5,set())

    def test_equal_volume_overlap_gap(self):
        parent,node,leaves=self.fixture()
        node['left']['cell'][0]=(F(0),F(3,5))
        node['right']['cell'][0]=(F(2,5),F(4,5))
        with self.assertRaises(AssertionError):r.tree(node,parent,leaves,5,set())

    def test_wrong_parent(self):
        parent,node,leaves=self.fixture()
        with self.assertRaises(AssertionError):r.tree(node,parent,leaves,6,set())

    def test_duplicate_active_leaf(self):
        parent,node,leaves=self.fixture();used=set()
        r.tree(node,parent,leaves,5,used)
        with self.assertRaises(AssertionError):r.tree(node,parent,leaves,5,used)

    def test_bad_split_axis(self):
        parent,node,leaves=self.fixture();node['axis']=-1
        with self.assertRaises(AssertionError):r.tree(node,parent,leaves,5,set())

    def test_compact_ancestor_leaf(self):
        parent,node,leaves=self.fixture()
        rows=[{'path':'0L','cell':parent},{'path':'0L1L','cell':parent}]
        with self.assertRaises(AssertionError):r.compact_partition(rows)

    def test_exact_pd_recheck(self):
        for c,out in zip(self.certs,self.pilot):
            self.assertIn(proof_gate(c,out),['SYLVESTER','GERSHGORIN'])

    def test_forged_pass_does_not_override_moment(self):
        c=self.certs[0];out=copy.deepcopy(self.pilot[0])
        out['moment_entry_enclosures'][0][0]=['100','100']
        out['ok']=True
        with self.assertRaises(AssertionError):proof_gate(c,out)

    def test_wrong_threshold(self):
        c=copy.deepcopy(self.certs[0]);c['threshold']='1/4'
        self.assertFalse(self.ck.certify_fixed(c)['ok'])
        with self.assertRaises(AssertionError):proof_gate(c,self.pilot[0])

    def test_rank_deficient_basis(self):
        c=copy.deepcopy(self.certs[0]);c['basis_num']=[[0]*3 for _ in range(4)]
        self.assertFalse(self.ck.certify_fixed(c)['ok'])

    def test_corrupt_center_inequality(self):
        c=copy.deepcopy(self.certs[0]);c['center_c']=['1000']*3
        self.assertFalse(self.ck.certify_fixed(c)['ok'])

    def test_false_full_replay_flag(self):
        fake={'complete':True,'failed':[],'total':12425,'processed':12425,
              'certificates':self.pilot}
        with self.assertRaises(AssertionError):completion_gate(fake,self.certs)

    def test_forged_gram(self):
        out=copy.deepcopy(self.pilot[0]);out['gram_matrix_exact'][0][0]='1'
        with self.assertRaises(AssertionError):proof_gate(self.certs[0],out)

    def test_true_sqrt3_exponent_bracket(self):
        self.assertLess(F(19,11)**2,3)
        self.assertGreater(F(26,15)**2,3)


if __name__=='__main__':unittest.main()
