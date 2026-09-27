#!/usr/bin/env python3
import json, math, hashlib
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar

# Secondary analysis of accepted exact N=6 microscopic correlation artifacts.
# Provenance (Git blob SHAs recorded through the connected GitHub source):
SOURCES = {
    "micro_to_12_semigroup_error_N6.json": "fd9425f16a0aa20e293a231ea5f6a171dcac8c7e",
    "localized12_MZ_fourier_N6.json": "0fb7a9e9daad2ea1e2f2cc4d58dc2c247ad297cb",
    "localized12_exact_mode_validation_N6.json": "d82d9cd827a39b952c45c0837454728d60c17094",
}

TIMES = np.array([0.0,0.1,0.25,0.5,1.0,2.0,4.0,8.0,16.0,32.0,64.0])
# Columns k=1,...,6. Values are exact equilibrium-projected microscopic correlations
# from micro_to_12_semigroup_error_N6.json.
C = np.array([
[1.0000000000000002,1.0000000000000002,1.0,1.0000000000000002,1.0000000000000002,1.0000000000000009],
[0.9878316176490285,0.9866222549141779,0.9908131124194399,0.9904179277265379,0.9893165895579813,0.9889253434676482],
[0.9749743290220367,0.972424859926432,0.9814284988451932,0.980556178154703,0.9782132101205531,0.97752438310305],
[0.9603148408644194,0.9561684581690341,0.9711339419904296,0.9696589549079745,0.96578389832034,0.9649361899311594],
[0.9395152159818849,0.9330606230192104,0.9568645920176717,0.9545259688246804,0.9483119292025027,0.9474926760685389],
[0.9066399831631925,0.8966035918913011,0.93414399073212,0.9305069191904872,0.9204605490601216,0.9198838231777777],
[0.850005122644495,0.8341723747248251,0.8940543634391346,0.8882328301466839,0.8716793516435882,0.8715550477947529],
[0.7498710684112068,0.7250218354181417,0.820750854490928,0.8112150235257789,0.783863725521548,0.7844312254079722],
[0.5840247316805685,0.5481241133035883,0.6920156721650234,0.6769713666629259,0.6342522389315838,0.6357773486434393],
[0.35426752466068967,0.31328867218741147,0.4919681386605228,0.47146671550097324,0.4152584279704751,0.41765430217485205],
[0.13035616288123322,0.10234727211141975,0.24864465508253827,0.22867178011010456,0.1780046762906335,0.18023570494002314],
])

LAMBDA_EXACT = np.array([
-0.03124318526101317,
-0.034961043176981066,
-0.02132466144971339,
-0.02261127518711381,
-0.026471600724339933,
-0.026262143621502046,
])
Q_EXACT = np.array([
0.000470645731792,
0.00135372597139,
0.00424680663344,
0.00386463648571,
0.00184808354015,
0.0013893428767,
])
RHO = 0.022618912154471746
TAU_STAR_297 = 0.5427059873

# Shell q_d -> Fourier lambda_k matrix.
MSHELL = np.zeros((6,6))
for k in range(1,7):
    for d in range(1,6):
        MSHELL[k-1,d-1] = 2.0*(math.cos(2*math.pi*k*d/12.0)-1.0)
    MSHELL[k-1,5] = (-1)**k - 1.0

idx = {float(t):i for i,t in enumerate(TIMES)}
intervals = [(0.1,0.25),(0.25,0.5),(0.5,1.0),(1.0,2.0),(2.0,4.0),(4.0,8.0),(8.0,16.0),(16.0,32.0),(32.0,64.0)]
slopes = {}
for k in range(6):
    arr=[]
    for a,b in intervals:
        s=(math.log(C[idx[b],k])-math.log(C[idx[a],k]))/(b-a)
        arr.append({"interval":[a,b],"slope":s})
    slopes[str(k+1)] = arr
    # log C convex => adjacent interval slopes are nondecreasing.
    ss=np.array([x["slope"] for x in arr])
    assert np.all(np.diff(ss) >= -1e-12)

# Three-point divided-slope drift; exact exponential (including arbitrary fixed residue Z) gives zero.
triples = [(0.1,0.25,0.5),(0.25,0.5,1.0),(0.5,1.0,2.0),(1.0,2.0,4.0),(2.0,4.0,8.0),(4.0,8.0,16.0),(8.0,16.0,32.0),(16.0,32.0,64.0)]
drift=[]
for a,b,c in triples:
    vals=[]
    for k in range(6):
        s1=(math.log(C[idx[b],k])-math.log(C[idx[a],k]))/(b-a)
        s2=(math.log(C[idx[c],k])-math.log(C[idx[b],k]))/(c-b)
        vals.append(s2-s1)
        assert s2-s1 >= -1e-12
    drift.append({"triple":[a,b,c],"by_mode":vals,"max":max(vals),"max_mode":int(np.argmax(vals)+1)})

late_slopes = np.array([(math.log(C[idx[64.0],k])-math.log(C[idx[32.0],k]))/32.0 for k in range(6)])
late_rel = np.abs((late_slopes-LAMBDA_EXACT)/LAMBDA_EXACT)
assert float(late_rel.max()) < 3e-12
q_late = np.linalg.solve(MSHELL, late_slopes)
q_late_rel = np.abs((q_late-Q_EXACT)/Q_EXACT)
assert float(q_late_rel.max()) < 4e-11

# Reconstruct 12-label reflection-symmetric projected transition row from C_k.
def p_from_modes(vals):
    hat=np.concatenate(([1.0],np.asarray(vals),np.asarray(vals)[-2::-1]))
    p=np.fft.ifft(hat).real
    assert abs(float(p.sum())-1.0)<1e-12
    assert float(p.min())>-1e-12
    return p

def p_markov(q,t):
    lam=MSHELL@np.asarray(q)
    hat=np.concatenate(([1.0],np.exp(lam*t),np.exp(lam[-2::-1]*t)))
    p=np.fft.ifft(hat).real
    assert abs(float(p.sum())-1.0)<1e-12
    return p

def noisy(p,eta):
    return (1.0-12.0*eta/11.0)*p + eta/11.0

def chernoff(p,q):
    lp=np.log(p); lq=np.log(q)
    f=lambda s: math.log(float(np.sum(np.exp(s*lp+(1.0-s)*lq))))
    rr=minimize_scalar(f,bounds=(0.0,1.0),method='bounded',options={'xatol':1e-14})
    return -float(rr.fun), float(rr.x)

P={float(t):p_from_modes(C[i]) for i,t in enumerate(TIMES)}
# Late incremental rates 32->64 cancel the initial-slip residue and define the late Markov generator.
# Use t=1 as an early held-out validation point from the pre-existing grid.
noise_rows=[]
for eta in [0.0,0.05,0.10,0.20]:
    p_true=noisy(P[1.0],eta)
    p_pred=noisy(p_markov(q_late,1.0),eta)
    ci,s=chernoff(p_true,p_pred)
    noise_rows.append({
        "eta":eta,
        "chernoff":ci,
        "chernoff_s":s,
        "shots_bound_5pct":math.ceil(math.log(10.0)/ci),
        "shots_bound_1pct":math.ceil(math.log(50.0)/ci),
        "tv":0.5*float(np.abs(p_true-p_pred).sum()),
    })

# 7-bin Y1 collapse at eta=5%, true vs late-Markov extrapolation at t=1.
orbits=[[0],[1,11],[2,10],[3,9],[4,8],[5,7],[6]]
p_true_05=noisy(P[1.0],0.05)
p_pred_05=noisy(p_markov(q_late,1.0),0.05)
seven_true=[float(p_true_05[o].sum()) for o in orbits]
seven_pred=[float(p_pred_05[o].sum()) for o in orbits]

# Illustrative delta-method sample size for a direct three-time rate-drift statistic.
# One-sided alpha=0.05, power=0.80 normal approximation, eta=5%.
# Triple 0,1,32; mode k=2. t=0 is measured with same readout channel, so attenuation cancels.
from scipy.stats import norm
eta=0.05
a_noise=1.0-12.0*eta/11.0
P0=p_from_modes(C[idx[0.0]])
Pobs={t:noisy(P[t],eta) for t in [0.0,1.0,32.0]}
def logmean_var_unit(t,k):
    j=np.arange(12)
    x=np.cos(2*np.pi*k*j/12.0)
    po=Pobs[t]
    m=float(np.dot(po,x))
    v=float(np.dot(po,x*x)-m*m)
    return m, v/(m*m)
a,b,c=0.0,1.0,32.0; k=2
s1=(math.log(C[idx[b],k-1])-math.log(C[idx[a],k-1]))/(b-a)
s2=(math.log(C[idx[c],k-1])-math.log(C[idx[b],k-1]))/(c-b)
delta=s2-s1
coeff=np.array([1/(b-a),-(1/(b-a)+1/(c-b)),1/(c-b)])
vars_u=np.array([logmean_var_unit(t,k)[1] for t in [a,b,c]])
var_equal=float(np.sum(coeff*coeff*vars_u))
zsum=float(norm.ppf(0.95)+norm.ppf(0.80))
M_equal=zsum*zsum*var_equal/(delta*delta)
w=np.abs(coeff)*np.sqrt(vars_u)
alloc=w/w.sum()
Mtot_opt=zsum*zsum*(float(w.sum())**2)/(delta*delta)

audit={
    "report":298,
    "N":6,
    "source_git_blob_shas":SOURCES,
    "theorem":{
        "correlation_form":"C(t)=sum_a w_a exp(-r_a t), w_a>=0",
        "log_convex":True,
        "effective_rate_derivative":"d[-d log C/dt]/dt = -Var_t(r) <= 0",
        "three_time_single_exponential_null":"divided-slope drift = 0",
        "time_independent_symmetric_readout_attenuation_cancels":True,
    },
    "shell_matrix":{
        "rank":int(np.linalg.matrix_rank(MSHELL)),
        "det":float(np.linalg.det(MSHELL)),
        "det_exact":"-3456*sqrt(3)",
    },
    "rho":RHO,
    "tau_star_297":TAU_STAR_297,
    "t_star_297_microscopic":TAU_STAR_297/RHO,
    "interval_slopes":slopes,
    "divided_slope_drift":drift,
    "strongest_mode_on_early_grid":2,
    "late_32_64":{
        "slopes":late_slopes.tolist(),
        "exact_eigenvalues":LAMBDA_EXACT.tolist(),
        "max_relative_eigenvalue_error":float(late_rel.max()),
        "q_from_incremental_slopes":q_late.tolist(),
        "q_exact_spectral":Q_EXACT.tolist(),
        "max_relative_q_error":float(q_late_rel.max()),
    },
    "late_generator_to_early_t1":{
        "seven_bin_true_eta_005":seven_true,
        "seven_bin_late_markov_prediction_eta_005":seven_pred,
        "noise_robustness":noise_rows,
        "note":"Chernoff shot bounds condition on an independently calibrated late generator; calibration uncertainty is not included."
    },
    "direct_rate_drift_normal_approx":{
        "triple":[0.0,1.0,32.0],
        "mode":2,
        "eta":0.05,
        "slope_0_1":s1,
        "slope_1_32":s2,
        "drift":delta,
        "equal_shots_per_time_for_alpha005_power080":M_equal,
        "equal_total_shots":3*M_equal,
        "optimal_total_shots_same_approx":Mtot_opt,
        "optimal_fraction_by_time":alloc.tolist(),
        "warning":"Delta-method normal approximation, not an exact finite-sample guarantee."
    },
    "scope_boundary":[
        "N=6 exact microscopic projected correlation artifact only",
        "does not prove N-uniform memory decay",
        "does not turn localized basins into physical spatial sites",
        "does not derive SI time or laboratory noise",
        "late-generator Chernoff counts assume independent shots and independent calibration"
    ]
}

out=Path('/mnt/data/fin298/multitime_semigroup_fingerprint_298.json')
out.write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps({
    "PASS":True,
    "max_late_lambda_relerr":float(late_rel.max()),
    "max_late_q_relerr":float(q_late_rel.max()),
    "tstar297":TAU_STAR_297/RHO,
    "k2_slopes":[x["slope"] for x in slopes["2"]],
    "eta005_chernoff":noise_rows[1]["chernoff"],
    "eta005_shots5":noise_rows[1]["shots_bound_5pct"],
    "direct_drift":delta,
    "direct_equal_M_per_time":M_equal,
    "direct_opt_total":Mtot_opt,
},indent=2))
