#!/usr/bin/env python3
"""Prospective composite FIN/flat test (DESIGN ONLY). No physical observations.
Conservative exclusion using a certified-by-formula gain mesh Lipschitz error.
Only applicable if actual calibration licenses these hypothesis sets.
"""
import math,json
from pathlib import Path
import numpy as np
from scipy.special import softmax
from source_model import all_models
ROOT=Path(__file__).resolve().parent
MODELS=('FIN_A7','FLAT_P7_TRACE')

def tv(a,b):return np.abs(np.array(a)-np.array(b)).sum()/2

def empirical_radius(n,k=13,alpha=.025):
    if n<=0:raise ValueError('positive sample size required')
    return math.sqrt(math.log((2**k-2)/alpha)/(2*n))

def lower_distance_to_hypothesis(emp,A,detector_tv=.005,lo=2.97,hi=3.03,intervals=240):
    """Lower bound on inf_{g in [lo,hi],TV(q,p_g)<=detector_tv} TV(emp,q).
    Arbitrary output TV-ball is an ENLARGEMENT of physically possible channels;
    hence the rejection result is conservative for the stated hypothetical set.
    """
    assert len(emp)==13 and min(emp)>=-1e-14 and abs(sum(emp)-1)<1e-9
    assert 0<=detector_tv<=1 and lo<=hi and intervals>0
    gg=np.linspace(lo,hi,intervals+1)
    row=A[0]
    dist=[]
    for g in gg:
        p=np.r_[softmax(g*row/2),0.]
        dist.append(tv(emp,p))
    delta=(hi-lo)/intervals
    # For this softmax: ||d p_g / d g||_TV <= range(row)/8.
    # Use twice this conservative derivative bound to avoid fragile factors.
    mesh_tv=float(np.ptp(row)/4*delta/2)
    bound=max(0.,float(min(dist))-mesh_tv-detector_tv)
    return {'rejection_distance_lower':bound,'min_grid_distance':float(min(dist)),
       'mesh_Lipschitz_TV_allowance':mesh_tv,'grid_points':len(gg),'detector_TV_ball_radius':detector_tv,
       'gain_box':[lo,hi]}

def classify(emp,n,alpha=.025,detector_tv=.005):
    eps=empirical_radius(n,13,alpha)
    stats={m:lower_distance_to_hypothesis(emp,all_models()[m],detector_tv) for m in MODELS}
    excluded={m:stats[m]['rejection_distance_lower']>eps for m in MODELS}
    if all(excluded.values()):out='NEITHER_MODEL_FITS'
    elif excluded['FIN_A7']:out='FLAT_ONLY'
    elif excluded['FLAT_P7_TRACE']:out='FIN_ONLY'
    else:out='INCONCLUSIVE'
    return {'outcome':out,'confidence_radius_TV':eps,'lower_distances':stats,
       'excluded_models':excluded,'conditional_sampling':'iid independent draws of complete 13-outcome record',
       'conditional_calibration':'independent certified gain [2.97,3.03], true observed law within TV .005 from at least one corresponding p_g; NO SUCH RECORDS PROVIDED'}

def main():
    A=all_models();n=40278676
    from scipy.special import softmax
    fin=np.r_[softmax(3*A['FIN_A7'][0]/2),0.]
    flat=np.r_[softmax(3*A['FLAT_P7_TRACE'][0]/2),0.]
    erasure=np.r_[np.ones(12)/12,0.]
    invalid=.03*erasure+.97*fin;invalid[-1]+=.03;invalid[:-1]*=.97 # discard this crafted input; use below instead
    invalid=np.r_[.97*fin[:12],.03]
    cases={'FIN_law':fin,'FLAT_law':flat,'g0_erasure':erasure,'invalid_3pct':invalid,
           'midpoint':.5*(fin+flat)}
    scenarios={k:classify(v,n) for k,v in cases.items()}
    low_n=classify(fin,1000)
    result={'scenario_tag':'SYNTHETIC_SOFTWARE_DESIGN_TESTS_ONLY','sample_n_conditional':n,
            'hypothetical_gain_box':[2.97,3.03],'hypothetical_detector_TV_ball_each':.005,
            'scenario_results':scenarios,'low_n_fin_1000':low_n,
            'NOTE':'These input probability vectors are theoretical and not empirical samples. The sample n is a design requirement; no observations were collected.'}
    (ROOT/'RG-F04/SYNTHETIC_COMPOSITE_TESTS.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print({k:v['outcome'] for k,v in scenarios.items()},'n=1000',low_n['outcome'])
if __name__=='__main__':main()
