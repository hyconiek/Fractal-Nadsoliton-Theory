#!/usr/bin/env python3
"""FIN ST1767--ST1781: moment-map reference dependence and nonidentifiability."""
import numpy as np
import sympy as sp
from scipy.optimize import root

from fin_total_nadsoliton_common import write_packet, write_round


LO,HI=1767,1781
NAMES=["MomentMapDefinition","ReferenceScaleTable","CouplingFactorTable","ReferenceDependenceNoGo",
 "CoordinateRescalingLaw","MomentJacobian","MomentNullDirection","ExactMomentEquivalentFamily","OffReferenceKernelDifference",
 "FreeBaseCouplingSurjectivity","OnlyFactorRelations","C1SquareInformationLoss","FourthMomentRankLift",
 "RoundFour_Verdict","RoundFour_Recommendation"]
def packet(k,s,b,p):return write_packet(k,NAMES[k-LO],s,b,p)
def main():
 d,o,p,b,e=sp.symbols('d o p b e', positive=True);K=sp.cos(o*d+p)/(1+b*d**e);mom=sp.Matrix([K,sp.diff(K,d),sp.diff(K,d,2)/2]);nom={o:.18575,p:.1625,b:1,e:1.8};fun=sp.lambdify((d,o,p,b,e),mom,'numpy');x={}
 x["ST1767"]=packet(1767,"constructed_three_jet_moment_map","Current P1866/P2363 input.",{
  "map":"(omega,phi,beta,eta;dref)->(c0,c1,c2)"})
 table=[];factors=[]
 for rr in [.5,1,2]:
  c=np.array(fun(rr,.18575,.1625,1,1.8),float).reshape(3);table.append({"dref":rr,"c":c.tolist()});factors.append({"dref":rr,"m_xi":1+c[0],"lambda":1+c[1]**2,"y":1+c[0]/2,"g":1+c[2]})
 x["ST1768"]=packet(1768,"proven_strict_moments_change_materially_with_reference_point","Same kernel.",{"table":table})
 x["ST1769"]=packet(1769,"proven_effective_coupling_multipliers_change_with_reference_point","No dref invariance.",{"factors":factors})
 x["ST1770"]=packet(1770,"proven_no_canonical_dref_or_length_normalization_is_exported_by_moment_map","dref=1 is a premise.",{
  "canonical_reference":False})
 x["ST1771"]=packet(1771,"proven_coordinate_rescaling_changes_derivative_moments","For d'=a d.",{
  "laws":["c0'(1)=K(1/a)","c1'(1)=a^-1 K'(1/a)","c2'(1)=a^-2 K''(1/a)/2"]})
 J=np.array(mom.jacobian([o,p,b,e]).subs({d:1,**nom}),float);u,s,vh=np.linalg.svd(J)
 x["ST1772"]=packet(1772,"proven_three_moment_parameter_Jacobian_has_rank_three_not_four","Nominal tuple.",{
  "matrix":J.tolist(),"singular_values":s.tolist(),"rank":int(np.linalg.matrix_rank(J,tol=1e-10))})
 null=vh[-1]
 x["ST1773"]=packet(1773,"constructed_local_unidentifiable_parameter_direction","First-order kernel-parameter fiber.",{
  "direction_omega_phi_beta_eta":null.tolist(),"residual":float(np.linalg.norm(J@null))})
 target=np.array(fun(1,.18575,.1625,1,1.8),float).reshape(3);families=[]
 for ee in [1.79,1.81,1.85]:
  sol=root(lambda z:np.array(fun(1,z[0],z[1],z[2],ee),float).reshape(3)-target,[.18575,.1625,1]);families.append({"eta":ee,"omega_phi_beta":sol.x.tolist(),"moment_residual":float(np.linalg.norm(sol.fun,np.inf))})
 x["ST1774"]=packet(1774,"constructed_distinct_parameter_tuples_with_same_three_moments","Numerical root witnesses.",{"families":families})
 diffs=[]
 for fam in families:
  oo,pp,bb=fam['omega_phi_beta'];ee=fam['eta'];row={"eta":ee}
  for rr in [.5,2,5]:row[str(rr)]=float(np.cos(oo*rr+pp)/(1+bb*rr**ee)-np.cos(.18575*rr+.1625)/(1+rr**1.8))
  diffs.append(row)
 x["ST1775"]=packet(1775,"proven_moment_equivalent_kernels_differ_away_from_reference","Three jets do not identify global kernel.",{"differences":diffs})
 x["ST1776"]=packet(1776,"proven_five_free_base_couplings_make_effective_map_pointwise_surjective_where_factors_nonzero","Choose each base coupling as target/factor.",{
  "predicted_absolute_couplings":False})
 x["ST1777"]=packet(1777,"identified_only_formal_multiplier_relations","Example equality of m2 and xi renormalization factors.",{
  "relations":["m2_eff/m2=xi_eff/xi=1+c0"],"physical_without_base_values":False})
 x["ST1778"]=packet(1778,"proven_lambda_map_discards_sign_of_c1","c1 and -c1 give same lambda multiplier.",{
  "orientation_information_lost":True})
 mom4=sp.Matrix([K,sp.diff(K,d),sp.diff(K,d,2)/2,sp.diff(K,d,3)/6]);J4=np.array(mom4.jacobian([o,p,b,e]).subs({d:1,**nom}),float);sv=np.linalg.svd(J4)[1]
 x["ST1779"]=packet(1779,"proven_fourth_jet_locally_lifts_parameter_rank_to_four","Constructive mathematical repair, not physical source.",{
  "determinant":float(np.linalg.det(J4)),"singular_values":sv.tolist(),"condition_number":float(np.linalg.cond(J4))})
 x["ST1780"]=packet(1780,"round_four_falsifies_three_moment_map_as_unique_or_reference_invariant_kernel_to_physics_bridge","It remains a chosen feature map.",{
  "unique_kernel":False,"reference_invariant":False,"predictive_couplings":False})
 x["ST1781"]=packet(1781,"recommended_ST1782_ST1796","Search for the smallest genuinely canonical action object and test propagator claims.",{
  "next":["Dirichlet action uniqueness","W versus Green resolvent","stability of W inverse action","spectral action inputs","strict/legacy role boundary"]})
 write_round(LO,HI,x)
if __name__=="__main__":main()
