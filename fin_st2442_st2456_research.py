#!/usr/bin/env python3
"""Round I: inventory and strict-shape ternary Gibbs candidate."""
import numpy as np
from fin_st2442_st2531_common import data,law,write_packet,write_round
LO=2442;N=['NonlinearInventory','DNLSBoundary','NegativeInformationBoundary','HebbBoundary','STDPBoundary','AdaptivePotentialBoundary','StrictLoopTensor','TernaryHamiltonian','FinitePartition','GlobalFlipDuality','LinearResponse','NonGaussianWitness','PairLeakage','ThetaSourceGate','RoundOneRecommendation']
def p(k,s,b,x):return write_packet(k,N[k-LO],s,b,x)
def main():
 tau,X,Y,S=data();r01=law(.1);r001=law(.001);x={}
 x['ST2442']=p(2442,'Classified','Current exported nonlinear families.',{'families':['focusing/saturating DNLS','negative-information Landau','Hebb/Oja/BCM/STDP','operator potential','ternary Gibbs candidate']})
 x['ST2443']=p(2443,'Conditional','P507/P517 boundary.',{'explicit_three_site_statistic':False,'sign_coefficient_strict':False})
 x['ST2444']=p(2444,'Conditional','Active gain g is supplied.',{'strict_active_gain':False,'first_angular_term_order':12})
 x['ST2445']=p(2445,'No-go','Covariance rules are second-order.',{'connected_three_body_source':False})
 x['ST2446']=p(2446,'Conditional','Delay can source direction but is external and generally non-Hermitian.',{'strict_delay':False})
 x['ST2447']=p(2447,'Underdetermined','V is a family until sourced.',{'unique_third_jet':False})
 x['ST2448']=p(2448,'Proven strict shape','Loop tensor from W.',{'faces':220,'positive':int((tau>0).sum()),'D12_orbits':12})
 x['ST2449']=p(2449,'Constructed conditional','Minimal all-triangle non-Gaussian state law.',{'H':'-theta sum_f tau_f product_(i in f) x_i','free':['theta']})
 x['ST2450']=p(2450,'Proven exact finite','All 4096 configurations.',{'states':len(X),'theta01':r01})
 x['ST2451']=p(2451,'Proven','Global bit flip sends S to -S and theta to -theta.',{'Z_even':True,'S_mean_odd':True,'selects_sign':False})
 x['ST2452']=p(2452,'Proven','At theta=0, d<E[Y_f]>/dtheta=tau_f by cube-character orthogonality.',{'small_theta_residual':r001['linear_response_residual']})
 x['ST2453']=p(2453,'Proven witness','Nonzero theta produces nonzero triangle moments.',{'theta':.1,'tri_min':r01['tri_min'],'tri_max':r01['tri_max']})
 x['ST2454']=p(2454,'Computed boundary','Overlapping triangles induce pair correlations beyond first order.',{'theta':.1,'pair_max_abs':r01['pair_max_abs'],'one_point_max_abs':r01['mean_max_abs']})
 x['ST2455']=p(2455,'Gate fail','Strict W fixes tensor shape, not theta sign/magnitude or state preparation.',{'strict_theta':False,'strict_nonGaussian_law':False})
 x['ST2456']=p(2456,'Recommendation','Test symmetry/polarity and whether theta can arise from existing scalar strict data.',{'next':['D12 covariance','reflection','global flip','entropy/Fisher','source no-go']})
 write_round(2442,2456,x)
if __name__=='__main__':main()
