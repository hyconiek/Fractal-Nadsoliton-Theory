import copy,json,sys,tempfile,unittest
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import boundary_checker as chk

SPEC=ROOT/'certificates/R7P-048_boundary_proof_spec.json'
TREE=ROOT/'results/R7P-053_refined_boundary_cover.json'

class CheckerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.good=chk.check(SPEC,TREE)
    def test_054_independent_checker_passes(self):
        self.assertEqual(self.good['status'],'CHECK_PASS')
        self.assertEqual(self.good['leaves'],436)
        self.assertEqual(self.good['reason_counts']['LOCAL_EQUALITY_CERTIFICATE'],1)
    def _write(self,d,spec,tree):
        sp=Path(d)/'spec.json';tp=Path(d)/'tree.json';sp.write_text(json.dumps(spec,indent=2)+'\n');tp.write_text(json.dumps(tree,indent=2)+'\n');return sp,tp
    def test_remove_leaf_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            spec=json.loads(SPEC.read_text());tree=json.loads(TREE.read_text());tree['leaves'].pop()
            # keep original hash valid by copying exact original spec bytes
            sp=Path(d)/'spec.json';sp.write_bytes(SPEC.read_bytes());tp=Path(d)/'tree.json';tp.write_text(json.dumps(tree))
            with self.assertRaises(ValueError):chk.check(sp,tp)
    def test_corrupt_split_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            tree=json.loads(TREE.read_text());leaf=tree['leaves'][0];p=leaf['path'];leaf['path']=str((int(p[0])+1)%3)+p[1:]
            sp=Path(d)/'spec.json';sp.write_bytes(SPEC.read_bytes());tp=Path(d)/'tree.json';tp.write_text(json.dumps(tree))
            with self.assertRaises(ValueError):chk.check(sp,tp)
    def test_flip_reason_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            tree=json.loads(TREE.read_text())
            leaf=next(x for x in tree['leaves'] if x['reason']=='SAFE_A_NONPOS' and F(x['B_bounds'][0])<0)
            leaf['reason']='SAFE_B_NONNEG'
            sp=Path(d)/'spec.json';sp.write_bytes(SPEC.read_bytes());tp=Path(d)/'tree.json';tp.write_text(json.dumps(tree))
            with self.assertRaises(ValueError):chk.check(sp,tp)
    def test_altered_spectral_interval_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            spec=json.loads(SPEC.read_text());tree=json.loads(TREE.read_text())
            spec['frozen_polynomials']['spectral_intervals']['l3'][1]='3'
            sp,tp=self._write(d,spec,tree)
            with self.assertRaisesRegex(ValueError,'hash mismatch'):chk.check(sp,tp)

if __name__=='__main__':unittest.main(verbosity=2)
