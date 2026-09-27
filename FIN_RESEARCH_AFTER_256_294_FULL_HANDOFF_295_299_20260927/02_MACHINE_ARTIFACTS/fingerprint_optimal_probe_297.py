#!/usr/bin/env python3
"""Replay for FIN report 297: fingerprint-optimal probe.

Scope:
- effective 12-state D12-circulant chain from report 293;
- same-rho/same-total-exit equalized comparator;
- localized state-0 preparation;
- symmetric categorical readout error eta;
- independent shots.

The script deliberately uses N=3..7 for probe design and keeps N=8 held out.
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar

ROWS = [
    (3,[0.00179034869,0.00932054908,0.0227890128,0.0211310797,0.0117600081,0.00874797691]),
    (4,[0.00126917356,0.00479252953,0.0126549862,0.0116826528,0.0062071051,0.0046699612]),
    (5,[0.00079921456,0.00254190843,0.00729500082,0.00668854178,0.00337858745,0.00254628983]),
    (6,[0.000468550925,0.00135413049,0.00425133082,0.00386837877,0.0018485772,0.00138883579]),
    (7,[0.000258887038,0.000715830917,0.00247712206,0.00223477624,0.00100447616,0.000748823886]),
    (8,[0.000136873225,0.000373516559,0.0014327756,0.00128061241,0.00053897181,0.000397420823]),
]
ETA_DESIGN = 0.05


def comparator(q):
    q=np.asarray(q,float)
    kz3=q[0]+q[1]+q[3]+q[4]
    exit_rate=2*np.sum(q[:5])+q[5]
    a=kz3/4.0
    b=(exit_rate-8*a)/3.0
    qc=np.array([a,a,b,a,a,b],float)
    assert np.all(qc>0)
    assert abs((qc[0]+qc[1]+qc[3]+qc[4])-kz3)<1e-15
    assert abs((2*np.sum(qc[:5])+qc[5])-exit_rate)<1e-15
    return qc


def generator_first_row(q):
    q=np.asarray(q,float)
    c=np.zeros(12,float)
    for d in range(1,6):
        c[d]+=q[d-1]; c[-d]+=q[d-1]
    c[6]+=q[5]
    c[0]=-np.sum(c[1:])
    return c


def transition_row(q,t):
    # Exact spectral evaluation of exp(Qt)[0,:] for the real symmetric circulant Q.
    lam=np.fft.fft(generator_first_row(q)).real
    p=np.fft.ifft(np.exp(lam*t)).real
    p[p<0]=np.maximum(p[p<0],-1e-15)
    p=np.maximum(p,0.0); p/=p.sum()
    return p


def rho(q):
    q=np.asarray(q,float)
    return 3.0*(q[0]+q[1]+q[3]+q[4])


def readout_noise(p,eta):
    # P(report=j|true=i)=1-eta if j=i, eta/11 otherwise.
    p=np.asarray(p,float)
    return (1-eta)*p + (eta/11.0)*(1-p)


def mode_pushforward(p,k):
    j=np.arange(12)
    vals=np.round(np.cos(2*np.pi*k*j/12),12)
    uniq=np.array(sorted(set(vals)))
    probs=np.array([p[vals==u].sum() for u in uniq])
    return uniq, probs


def chernoff(p,q):
    p=np.asarray(p,float); q=np.asarray(q,float)
    assert np.all(p>0) and np.all(q>0)
    def overlap(s):
        return float(np.sum(np.exp(s*np.log(p)+(1-s)*np.log(q))))
    r=minimize_scalar(overlap,bounds=(1e-10,1-1e-10),method='bounded',options={'xatol':1e-13})
    return -math.log(r.fun), float(r.x)


def distributions(N,q,tau,eta):
    qc=comparator(q)
    rr=rho(q)
    t=tau/rr
    pa=readout_noise(transition_row(q,t),eta)
    pc=readout_noise(transition_row(qc,t),eta)
    return pa,pc


def C_mode(N,q,tau,k,eta):
    pa,pc=distributions(N,q,tau,eta)
    _,a=mode_pushforward(pa,k); _,c=mode_pushforward(pc,k)
    return chernoff(a,c)[0]


def C_full(N,q,tau,eta):
    pa,pc=distributions(N,q,tau,eta)
    return chernoff(pa,pc)[0]


def optimize_mode(k, train, eta):
    grid=np.linspace(0.10,1.50,701)
    scores=[]
    for tau in grid:
        scores.append(min(C_mode(N,q,tau,k,eta) for N,q in train))
    i=int(np.argmax(scores)); tc=float(grid[i])
    lo=max(0.02,tc-0.03); hi=tc+0.03
    f=lambda x: -min(C_mode(N,q,x,k,eta) for N,q in train)
    r=minimize_scalar(f,bounds=(lo,hi),method='bounded',options={'xatol':1e-12})
    tau=float(r.x)
    cs=[C_mode(N,q,tau,k,eta) for N,q in train]
    return tau,min(cs),cs


def shell_map_matrix():
    M=np.zeros((6,6),float)
    for k in range(1,7):
        for d in range(1,6):
            M[k-1,d-1]=2*(math.cos(2*math.pi*k*d/12)-1)
        M[k-1,5]=((-1)**k-1)
    return M


def main():
    train=ROWS[:-1]
    test=ROWS[-1]

    # exact algebraic shell-map determinant is -3456*sqrt(3)
    M=shell_map_matrix()
    det=float(np.linalg.det(M))
    det_exact=-3456*math.sqrt(3)
    assert np.linalg.matrix_rank(M)==6
    assert abs(det-det_exact)<1e-9

    # Scan one-cosine observables. k=1 and k=5 are reflection-orbit encodings.
    mode_scan={}
    for k in range(1,7):
        tau,score,cs=optimize_mode(k,train,ETA_DESIGN)
        mode_scan[str(k)]={"tau_star":tau,"worst_train_C":score,"train_C":cs}

    tau_star=mode_scan['1']['tau_star']
    train_rows=[]
    for N,q in train:
        C=C_mode(N,q,tau_star,1,ETA_DESIGN)
        Cfull=C_full(N,q,tau_star,ETA_DESIGN)
        assert abs(C-Cfull)<2e-13
        train_rows.append({
            "N":N,
            "C":C,
            "chernoff_s":chernoff(mode_pushforward(distributions(N,q,tau_star,ETA_DESIGN)[0],1)[1],
                                  mode_pushforward(distributions(N,q,tau_star,ETA_DESIGN)[1],1)[1])[1],
            "M_for_Pe_le_0p05_bound":math.ceil(math.log(10)/C),
            "M_for_Pe_le_0p01_bound":math.ceil(math.log(50)/C),
        })

    N8,q8=test
    pa8,pc8=distributions(N8,q8,tau_star,ETA_DESIGN)
    vals8,da8=mode_pushforward(pa8,1); _,dc8=mode_pushforward(pc8,1)
    C8,s8=chernoff(da8,dc8)
    assert abs(C8-C_full(N8,q8,tau_star,ETA_DESIGN))<2e-13

    noise_robust=[]
    for eta in (0.0,0.05,0.10,0.20):
        ctr=[C_mode(N,q,tau_star,1,eta) for N,q in train]
        c8=C_mode(N8,q8,tau_star,1,eta)
        noise_robust.append({
            "eta":eta,
            "worst_train_C":min(ctr),
            "heldout_N8_C":c8,
            "worst_train_M_for_Pe_le_0p05_bound":math.ceil(math.log(10)/min(ctr)),
            "heldout_N8_M_for_Pe_le_0p05_bound":math.ceil(math.log(10)/c8),
        })

    timing_robust=[]
    for factor in (0.8,0.9,1.0,1.1,1.2):
        tau=tau_star*factor
        ctr=[C_mode(N,q,tau,1,ETA_DESIGN) for N,q in train]
        c8=C_mode(N8,q8,tau,1,ETA_DESIGN)
        timing_robust.append({"tau_factor":factor,"tau":tau,"worst_train_C":min(ctr),"heldout_N8_C":c8})

    result={
        "status":"P1-297 positive for the declared report-293 comparator; universal minimax separation over all same-rho/same-exit alternatives is impossible without a minimum alternative distance.",
        "noise_model":{
            "description":"Independent categorical readout error: correct label with probability 1-eta; otherwise uniformly one of the other 11 labels.",
            "eta_design":ETA_DESIGN,
        },
        "design_split":{"train_N":[3,4,5,6,7],"heldout_N":[8]},
        "shell_map":{"rank":6,"determinant_numeric":det,"determinant_exact":"-3456*sqrt(3)"},
        "mode_scan":mode_scan,
        "selected_probe":{
            "preparation":"localized effective state J=0",
            "observable":"Y=cos(2*pi*J/12), retain its seven-bin outcome histogram",
            "tau_star":tau_star,
            "time":"t=tau_star/rho calibrated on an independent dataset",
            "reason":"Y labels the seven reflection orbits exactly; for reflection-symmetric transition rows it is a sufficient statistic and preserves the full-label Chernoff information.",
        },
        "training_results":train_rows,
        "heldout_N8":{
            "C":C8,"chernoff_s":s8,
            "M_for_Pe_le_0p05_bound":math.ceil(math.log(10)/C8),
            "M_for_Pe_le_0p01_bound":math.ceil(math.log(50)/C8),
            "Y_values":vals8.tolist(),
            "FIN_probs":da8.tolist(),
            "comparator_probs":dc8.tolist(),
        },
        "noise_robustness_fixed_probe":noise_robust,
        "timing_robustness_fixed_probe":timing_robust,
        "chernoff_bound":"For equal priors and M independent shots, optimal Bayes error satisfies the standard Chernoff upper bound Pe <= 0.5*exp(-M*C).",
    }
    out=Path(__file__).with_name('fingerprint_optimal_probe_297.json')
    out.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print('PASS')
    print('tau_star',tau_star)
    print('worst_train_C',min(x['C'] for x in train_rows))
    print('heldout_N8_C',C8)
    print('heldout_N8_M_5pct_bound',math.ceil(math.log(10)/C8))
    print('shell_map_det',det)
    print('output',out)

if __name__=='__main__':
    main()
