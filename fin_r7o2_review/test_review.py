from fractions import Fraction as F
from itertools import product
import math, random, unittest
from fin_r7o2_review.intervals_fast import Interval as I
from fin_r7o2_review.review import partition

class ReviewTests(unittest.TestCase):
    def test_rational_inputs_enclosed(self):
        for x in [F(1,3),F(-1,7),F(10**30+1,10**30),F(1,10**310)]:
            a=I(x);self.assertLessEqual(a.lo,x);self.assertGreaterEqual(a.hi,x)
    def test_arithmetic_encloses_exact_endpoints(self):
        rng=random.Random(20260920)
        for _ in range(300):
            a,b=sorted([F(rng.randint(-100,100),37),F(rng.randint(-100,100),41)])
            c,d=sorted([F(rng.randint(-100,100),43),F(rng.randint(-100,100),47)])
            A,B=I(a,b),I(c,d)
            for result,fn in [(A+B,lambda x,y:x+y),(A-B,lambda x,y:x-y),(A*B,lambda x,y:x*y)]:
                vals=[fn(x,y) for x,y in product([a,b],[c,d])]
                self.assertLessEqual(result.lo,min(vals));self.assertGreaterEqual(result.hi,max(vals))
            if not c<=0<=d:
                result=A/B;vals=[x/y for x,y in product([a,b],[c,d])]
                self.assertLessEqual(result.lo,min(vals));self.assertGreaterEqual(result.hi,max(vals))
    def test_zero_denominator_rejected(self):
        with self.assertRaises(ValueError):I(1)/I(-1,1)
    def test_nonfinite_rejected(self):
        with self.assertRaises(ValueError):I(F(10)**400)
    def test_exact_geometry_and_missing_leaf(self):
        parent=((F(0),F(1)),)*4
        left=list(parent);right=list(parent);left[0]=(F(0),F(1,3));right[0]=(F(1,3),F(1))
        rows=[{'path':'L','cell':left},{'path':'R','cell':right}]
        partition(parent,rows)
        with self.assertRaises(AssertionError):partition(parent,rows[:1])
    def test_equal_volume_gap_overlap_rejected(self):
        parent=((F(0),F(1)),)*4
        left=list(parent);left[0]=(F(0),F(1,2))
        with self.assertRaises(AssertionError):partition(parent,[{'path':'L','cell':left},{'path':'R','cell':left}])

if __name__=='__main__':unittest.main()
