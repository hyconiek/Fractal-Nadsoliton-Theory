import math,sys,unittest
from fractions import Fraction as F
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.model import *
from src.derivatives import dual_all,directional_34
from src.intervals import QI,sqrt_interval,exp_interval,log_interval,tanh_interval,sech_interval

class BTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.W,cls.A,cls.L,cls.X,cls.C,cls.A7=feature_spaces()
 def test_009_factorization_embedding_formula(self):
  self.assertLess(np.linalg.norm(self.A7-self.X@self.X.T),1e-14)
  self.assertLess(np.linalg.norm(self.X.sum(axis=0)),1e-14)
  np.testing.assert_allclose(self.C,self.X[:,[0,2,4,6]],atol=0,rtol=0)
  np.testing.assert_allclose(self.X,canonical_feature_formula(self.L),atol=1e-15)
  norms=np.diag(self.X.T@self.X)
  np.testing.assert_allclose(norms,[self.L[3],self.L[3],self.L[4],self.L[4],self.L[5],self.L[5],self.L[6]],rtol=2e-14)
 def test_009_mutation_detected(self):
  bad=self.X.copy();bad[:,3]*=-1
  self.assertGreater(np.linalg.norm(bad-canonical_feature_formula(self.L)),1.0)
 def test_009_field_routes(self):
  rng=np.random.default_rng(9009)
  for _ in range(20):
   th=rng.normal(size=7); s=rng.normal(size=4)
   np.testing.assert_allclose(self.X@th,sum(th[k]*self.X[:,k] for k in range(7)),atol=1e-14)
   np.testing.assert_allclose(self.C@s,self.X[:,[0,2,4,6]]@s,atol=1e-14)
 def test_010_group_and_equivariance(self):
  acts=d12_actions(self.X); self.assertEqual(len(acts),24)
  keys=list(acts)
  for k,(P,T) in acts.items():
   np.testing.assert_allclose(P@self.X,self.X@T,atol=2e-14)
   np.testing.assert_allclose(T.T@T,np.eye(7),atol=3e-14)
  mats=[v[0] for v in acts.values()]
  for P1 in mats:
   for P2 in mats:
    Q=P1@P2
    self.assertTrue(any(np.array_equal(Q,R) for R in mats))
  rng=np.random.default_rng(9010); th=rng.normal(size=7);g=3.8
  base=dual7(th,g,self.X)[0]
  for P,T in acts.values(): self.assertAlmostEqual(base,dual7(T@th,g,self.X)[0],places=12)
 def test_011_joint_min_identities(self):
  rng=np.random.default_rng(9011);g=3.7
  for _ in range(10):
   th=rng.normal(size=7);phi,grad,H,p=dual7(th,g,self.X)
   # joint F evaluated at softmax p equals dual Phi
   D=np.sum(p*np.log(12*p)); joint=D+th@th/(2*g)-th@(self.X.T@p)
   self.assertAlmostEqual(phi,joint,places=13)
   # theta minimizer at fixed p reproduces primal
   thp=g*self.X.T@p
   joint2=D+thp@thp/(2*g)-thp@(self.X.T@p)
   self.assertAlmostEqual(joint2,primal(p,g,self.X),places=13)
 def test_012_derivatives(self):
  rng=np.random.default_rng(9012);th=rng.normal(size=7)*.3;v=rng.normal(size=7);v/=np.linalg.norm(v);g=3.8
  val,grad,H,T3,T4,p=dual_all(th,g,self.X)
  eps=1e-2
  def f(t): return dual_all(th+t*v,g,self.X)[0]
  # 3rd central stencil and 4th central stencil diagnostics
  d3=(f(2*eps)-2*f(eps)+2*f(-eps)-f(-2*eps))/(2*eps**3)
  d4=(f(2*eps)-4*f(eps)+6*f(0)-4*f(-eps)+f(-2*eps))/eps**4
  a3,a4=directional_34(th,g,self.X,v)
  self.assertLess(abs(d3-a3),3e-5)
  self.assertLess(abs(d4-a4),2e-5)
  np.testing.assert_allclose(T3,np.transpose(T3,(1,0,2)),atol=1e-14)
 def test_013_interval_primitives(self):
  import mpmath as mp; mp.mp.dps=60
  for a,b in [(F(-2),F(3)),(F(0),F(10)),(F(1,10),F(2))]:
   E=exp_interval(QI(a,b)); self.assertLessEqual(float(E.lo),float(mp.e**mp.mpf(str(float(a)))));self.assertGreaterEqual(float(E.hi),float(mp.e**mp.mpf(str(float(b)))))
  for a,b in [(F(1,10),F(10)),(F(1),F(3))]:
   L=log_interval(QI(a,b)); self.assertLessEqual(float(L.lo),math.log(float(a)));self.assertGreaterEqual(float(L.hi),math.log(float(b)))
  T=tanh_interval(QI(-2,3)); self.assertLessEqual(float(T.lo),math.tanh(-2));self.assertGreaterEqual(float(T.hi),math.tanh(3))
  S=sech_interval(QI(-2,3)); self.assertLessEqual(float(S.lo),1/math.cosh(3));self.assertGreaterEqual(float(S.hi),1.0)
  R=sqrt_interval(QI(F(2),F(3))); self.assertLessEqual(float(R.lo),math.sqrt(2));self.assertGreaterEqual(float(R.hi),math.sqrt(3))
  with self.assertRaises(ValueError): log_interval(QI(-1,2))
 def test_014_inertia_stationary_example(self):
  # solve known C4 root embedded in full X only to test Schur inertia identity algebraically
  from scipy.optimize import root
  g=3.71834489812038
  s0=np.array([1.8199035813,1.9139895547,1.9145691325,1.3672032802])
  # full theta with sine zeros is stationary because reflection-fixed representative
  th=np.zeros(7);th[[0,2,4,6]]=s0
  phi,grad,H,p=dual7(th,g,self.X)
  self.assertLess(np.linalg.norm(grad),2e-8)
  # Tangent basis 12x11 orthogonal to ones
  Q=np.linalg.qr(np.column_stack([np.ones(12),np.eye(12)[:,1:]]))[0][:,1:]
  Hp=Q.T@(np.diag(1/p)-g*self.X@self.X.T)@Q
  ni7=sum(np.linalg.eigvalsh(H)<-1e-8); nip=sum(np.linalg.eigvalsh(Hp)<-1e-8)
  nz7=sum(abs(np.linalg.eigvalsh(H))<=1e-8); nzp=sum(abs(np.linalg.eigvalsh(Hp))<=1e-8)
  self.assertEqual((ni7,nz7),(nip,nzp))

if __name__=='__main__':unittest.main()
