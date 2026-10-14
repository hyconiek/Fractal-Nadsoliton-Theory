#!/usr/bin/env python3
"""Read-only source intake and independent, bounded three-copy verification.

The new campaign's model/propagation code is NOT imported. The versioned
safe-seed defines the numerical proposal. Exact-decimal moment certificates
use mpmath.iv, dyadic candidate functions, and the analytic hidden gap >= 1.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import importlib.util
import itertools
import json
import math
import re
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
from scipy.linalg import expm
from scipy.sparse import coo_matrix, diags, eye
from scipy.special import logsumexp, softmax

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "fin_generative_memory_campaign_20261014"
Q, N, Y, G = 12, 1728, 144, 5.145228719489142


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(base, name):
    return json.loads((base / name).read_text())


def endpoints(value):
    out = []
    for sign, mantissa, exponent, _ in value._mpi_:
        n = -mantissa if sign else mantissa
        out.append(Fraction(n << exponent, 1) if exponent >= 0 else Fraction(n, 1 << -exponent))
    return out


def decimal(value, digits=15, upper=True):
    n, d = value.numerator * 10**digits, value.denominator
    r = -((-n) // d) if upper else n // d
    sign = "-" if r < 0 else ""
    a, b = divmod(abs(r), 10**digits)
    return f"{sign}{a}.{b:0{digits}d}"


def box(value, digits=15):
    a, b = endpoints(value)
    return [decimal(a, digits, False), decimal(b, digits, True)]


def dyadic(value):
    f = Fraction.from_float(float(value))
    return mp.iv.mpf(f.numerator) / mp.iv.mpf(f.denominator)


def build(A):
    states = np.array(list(itertools.product(range(Q), repeat=3)))
    count = np.eye(Q, dtype=int)[states].sum(axis=1)
    logpi = (G / 6) * np.einsum("ni,ij,nj->n", count, A, count)
    pi = np.exp(logpi - logsumexp(logpi))
    y = states[:, 0] * Q + states[:, 1]
    py = np.bincount(y, weights=pi, minlength=Y)
    rows, cols, values = [], [], []
    third_rows, third_cols, third_values = [], [], []
    for axis, stride in enumerate((Q**2, Q, 1)):
        other = count - np.eye(Q, dtype=int)[states[:, axis]]
        rate = softmax((G / 3) * (other @ A), axis=1)
        for new in range(Q):
            dest = np.arange(N) + (new - states[:, axis]) * stride
            rows.extend(range(N)); cols.extend(dest); values.extend(rate[:, new])
            if axis == 2:
                third_rows.extend(range(N)); third_cols.extend(dest); third_values.extend(rate[:, new])
    L = coo_matrix((values, (rows, cols)), shape=(N, N)).tocsr() - 3 * eye(N, format="csr")
    P3 = coo_matrix((third_values, (third_rows, third_cols)), shape=(N, N)).tocsr()
    U = coo_matrix((np.sqrt(pi / py[y]), (np.arange(N), y)), shape=(N, Y)).tocsr()
    S = diags(np.sqrt(pi)) @ L @ diags(1 / np.sqrt(pi))
    Acoarse = (U.T @ S @ U).toarray()
    C = S @ U - U @ Acoarse
    return states, pi, py, y, L, S.tocsr(), U, Acoarse, C, P3


def propagate(L, p, t):
    # Different uniformization rate and cutoff from the source's lambda=3,m=65.
    P = eye(N, format="csr") + L / 4
    w = math.exp(-4 * t)
    term, total = p.copy(), w * p
    for k in range(1, 91):
        term = P.T @ term
        w *= 4 * t / k
        total += w * term
    return total


def exact_inputs():
    mp.iv.dps = 55
    I = mp.iv
    gain, omega, phi, eta = map(I.mpf, ("5.145228719489142", "0.18575", "0.1625", "1.8"))
    weights = [I.mpf(0)] + [I.cos(omega * min(d, 12-d) + phi) / (1 + I.power(I.mpf(min(d, 12-d)), eta)) for d in range(1, Q)]
    lam = {k: sum((weights[d] * (1-I.cos(2*I.pi*k*d/12)) for d in range(1, Q)), I.mpf(0)) for k in (3, 4, 5, 6)}
    row = [sum((lam[k] * I.cos(2*I.pi*k*d/12)/6 for k in (3, 4, 5)), I.mpf(0)) + lam[6]*(-1)**d/12 for d in range(Q)]
    rates = {}
    for a, b in itertools.product(range(Q), repeat=2):
        rr = [I.exp(gain/3 * (row[(j-a)%Q] + row[(j-b)%Q])) for j in range(Q)]
        denom = sum(rr, I.mpf(0))
        rates[a, b] = [x/denom for x in rr]
    ww = []
    for state in itertools.product(range(Q), repeat=3):
        ww.append(I.exp(gain/6 * sum((row[(a-b)%Q] for a in state for b in state), I.mpf(0))))
    total = sum(ww, I.mpf(0))
    pi = [w/total for w in ww]
    py = [sum(pi[12*y:12*y+12], I.mpf(0)) for y in range(Y)]
    return row, rates, pi, py


def moment_certificate(solution_function, exact):
    """Validate a dyadic function via an interval residual, not CG convergence.

    Use indicator f=1{Y=00}; divide weighted quadratic channels by piY00
    afterwards. This avoids treating nominal sqrt(pi) as exact input.
    """
    _, rates, pi, py = exact
    I = mp.iv
    zero = I.mpf(0)
    states = list(itertools.product(range(Q), repeat=3))

    def project(v):
        means = [sum((pi[12*y+h]*v[12*y+h] for h in range(Q)), zero)/py[y] for y in range(Y)]
        return [v[x] - means[x//12] for x in range(N)]

    def backward(v):
        result = []
        for x, (a, b, c) in enumerate(states):
            targets = ((rates[b,c], 12*b+c, 144),
                       (rates[a,c], 144*a+c, 12),
                       (rates[a,b], 144*a+12*b, 1))
            result.append(sum((rr[j]*(v[start+stride*j]-v[x]) for rr,start,stride in targets for j in range(Q)), zero))
        return result

    f = [I.mpf(int(x < 12)) for x in range(N)]
    h = project(backward(f))
    psi = project([dyadic(v) for v in solution_function])
    Lpsi = project(backward(psi))
    r = [h[x]+Lpsi[x] for x in range(N)]  # r=h-H psi
    dot = lambda a,b: sum((pi[x]*a[x]*b[x] for x in range(N)), zero)/py[0]
    # Interval squaring must not be multiplication of two independent boxes.
    norm2 = lambda a: sum((pi[x]*I.power(a[x],2) for x in range(N)), zero)/py[0]
    k0 = norm2(h)
    rho = I.sqrt(norm2(r))
    m0hat, m1hat = dot(h, psi), norm2(psi)
    err0 = I.sqrt(k0) * rho
    err1 = 2*I.sqrt(m1hat)*rho + I.power(rho,2)
    er0, er1 = endpoints(err0)[1], endpoints(err1)[1]
    to_iv = lambda v: I.mpf(v.numerator)/I.mpf(v.denominator)
    m0 = m0hat + I.mpf([-1,1])*to_iv(er0)
    m1 = m1hat + I.mpf([-1,1])*to_iv(er1)
    assert endpoints(m0)[0] > 0 and endpoints(m1)[0] > 0
    return {
        "arithmetic": "mpmath.iv at 55 digits, exact-decimal model, exact dyadic candidate, public outward decimal endpoints",
        "backend_assumption": "correct mpmath.iv elementary/transcendental enclosures; not proof-assistant formalization",
        "analytic_gap_used": "1 <= H <= 3 on the exact hidden subspace",
        "K0_normalized_interval": box(k0),
        "residual_norm_upper": decimal(endpoints(rho)[1], 25),
        "M0_integral_interval": box(m0),
        "M1_first_integral_moment_interval": box(m1),
        "M0_residual_error_upper": decimal(er0, 25),
        "M1_residual_error_upper": decimal(er1, 25),
        "candidate_function_sha256": hashlib.sha256(np.asarray(solution_function,dtype="<f8").tobytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replay", type=Path)
    parser.add_argument("--witness-input", type=Path,
                        help="recheck a saved exact dyadic moment candidate")
    args = parser.parse_args()
    manifest = []
    for line in (SOURCE/"MANIFEST.sha256").read_text().splitlines():
        expected, name = line.split(None, 1)
        assert sha(SOURCE/name) == expected, name
        manifest.append(name)
    uncovered = sorted(str(p.relative_to(SOURCE)) for p in SOURCE.rglob("*") if p.is_file() and str(p.relative_to(SOURCE)) not in manifest and p.name != "MANIFEST.sha256")
    assert not uncovered, uncovered
    stage_input_count = 0
    for stage in sorted(SOURCE.glob("RG-D0[1-5]")):
        for line in (stage/"INPUTS.sha256").read_text().splitlines():
            expected, relative = line.split(None,1)
            assert sha(stage/relative) == expected, (stage,relative)
            stage_input_count += 1
    lock = read(SOURCE,"RG-D04/LOCK.json")
    assert all(sha(SOURCE/f)==h for f,h in lock["pretest_files_sha256"].items())
    frozen = read(SOURCE,"RG-D04/FROZEN_PREDICTOR.json")
    assert all(sha(SOURCE/f)==h for f,h in frozen["input_hashes"].items())
    core = ROOT/"fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"
    spec = importlib.util.spec_from_file_location("versioned_core",core)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    _,_,_,_,A = module.build_rank7()
    proposals = module.build_countermodels(A)
    kernels = {name:proposals[name] for name in ("FIN_A7","FULL_POTTS_TRACE","FLAT_P7_TRACE")}
    pred = read(SOURCE,"RG-D04/PREDICTIONS_PREOPEN.json")
    supplied = read(SOURCE,"RG-D05/RAW_MODEL_OUTPUTS.json")
    envelopes = read(SOURCE,"RG-D04/MODEL_ERROR_ENVELOPES_PREOPEN.json")
    rows, diagnostics, moments = [], {}, {}
    solution = None
    for name,kernel in kernels.items():
        states,pi,py,y,L,S,U,Acoarse,C,P3 = build(kernel)
        sqrtpi = np.sqrt(pi)
        qproj = lambda v: v-U@(U.T@v)
        H = lambda v: -qproj(S@qproj(v))
        c = np.asarray(C[:,0]).ravel()
        x = np.zeros(N)
        for _ in range(90):
            x = qproj(x + (c-H(x))/2)
        moments[name] = {"K0":float(c@c), "Kprime0":float(c@(S@c)), "M0_Richardson":float(c@x), "M1_Richardson":float(x@x), "residual_norm_numeric":float(np.linalg.norm(c-H(x)))}
        for field,mkey in (("K0","K0_00"),("Kprime0","Kprime0_00"),("M0_Richardson","M0_integral_00_numeric"),("M1_Richardson","M1_first_integral_moment_00_numeric")):
            assert abs(moments[name][field]-read(SOURCE,"RG-D03/MEMORY_REDUCTION.json")[name][mkey]) < 1e-12
        E3 = diags(sqrtpi)@P3@diags(1/sqrtpi)
        defect = E3-U@U.T
        identity = (U.T@(S@(S@U))).toarray()-(Acoarse@Acoarse)-C.T@C
        diagnostics[name] = {"P3_equals_UUt_max":float(np.max(np.abs(defect.data))), "semigroup_identity_max":float(np.max(np.abs(identity))), "L_rowsum_max":float(np.max(np.abs(L.sum(axis=1)))), "S_symmetry_max":float(np.max(np.abs((S-S.T).data)))}
        assert max(diagnostics[name].values()) < 1e-12
        if name == "FIN_A7":
            solution = x/sqrtpi*np.sqrt(py[0])
        observe = lambda p: np.bincount(y,weights=p,minlength=Y)
        a0 = np.zeros(Y); a0[0]=1/np.sqrt(py[0])
        btime = {t:np.sqrt(py)*(expm(t*Acoarse)@a0) for t in (0.35,0.8)}
        bfirst = expm(.35*Acoarse)@a0
        event = np.zeros(Y);event[0]=bfirst[0]
        bjoint = np.r_[np.sqrt(py)*(expm(.45*Acoarse)@event),np.sqrt(py)*(expm(.45*Acoarse)@(bfirst-event))]
        for prep in ("conditional_00","delta_0","delta_1","delta_2"):
            p0 = np.zeros(N)
            if prep == "conditional_00": p0[:12]=pi[:12]/py[0]
            else: p0[int(prep.split("_")[-1])]=1
            vals = pred[name][prep]
            source = supplied[name][prep]
            pfirst = propagate(L,p0,.35)
            yes = pfirst.copy();yes[12:]=0
            joint = np.r_[observe(propagate(L,yes,.45)),observe(propagate(L,pfirst-yes,.45))]
            cases = [("one_time_"+str(t), observe(propagate(L,p0,t)), source["one_time_reference_full144"][str(t)], vals["single_time_full144"][str(t)], btime[t], envelopes[name][prep]["single_time"][str(t)]["TV_upper_formula"]) for t in (.35,.8)]
            cases.append(("two_time_0.35_0.8", joint,source["joint_reference_288"],vals["joint_time_first_binary_second144"],bjoint,envelopes[name][prep]["two_time"]["TV_upper_formula"]))
            for scope,ref,saved,forecast,baseline,limit in cases:
                gap = float(np.abs(ref-np.asarray(forecast)).sum()/2)
                agreement = float(np.max(np.abs(ref-np.asarray(saved))))
                assert agreement < 1e-12 and gap < .01 and gap < limit+1.1e-7
                assert min(forecast) >= 0 and abs(math.fsum(forecast)-1) < 1e-11
                rows.append({"model":name,"preparation":prep,"scope":scope,"TV_forecast":gap,"TV_markov":float(np.abs(ref-baseline).sum()/2),"saved_reference_max_discrepancy":agreement,"formula_plus_assumed_allowance":limit+1.1e-7})
    if args.witness_input:
        witness = read(args.witness_input.parent,args.witness_input.name)
        solution = np.array([float.fromhex(v) for v in witness["candidate_function_hex"]])
        assert len(solution) == N
        assert hashlib.sha256(np.asarray(solution,dtype="<f8").tobytes()).hexdigest() == witness["candidate_function_sha256"]
    certificate = moment_certificate(solution,exact_inputs())
    assert Fraction(certificate["residual_norm_upper"]) < Fraction(1,10**12)
    witness = {"candidate_function_hex":[float(v).hex() for v in solution],
               "candidate_function_sha256":certificate["candidate_function_sha256"],
               "meaning":"exact dyadic backward function before exact conditional Q projection; unnormalized indicator channel"}
    args.output.with_name("moment_witness.json").write_text(json.dumps(witness,indent=2)+"\n")
    I = mp.iv
    tail = I.exp(-I.mpf("4.5"))*I.mpf("4.5")**66/math.factorial(66)/(1-I.mpf("4.5")/67)
    assert endpoints(tail)[1] < Fraction(2832,10**55)
    normalized = lambda text: re.sub(r"\s+","",html.unescape(text).replace("\\", "").replace("*", ""))
    replay = None
    if args.replay:
        rs = read(args.replay,"RG-D05/VALIDATION.json")
        rp = read(args.replay,"RG-D04/PREDICTIONS_PREOPEN.json")
        delta = max(abs(a-b) for model in pred for prep in pred[model] for t in ("0.35","0.8") for a,b in zip(pred[model][prep]["single_time_full144"][t],rp[model][prep]["single_time_full144"][t]))
        replay = {"directory":str(args.replay),"checks":read(args.replay,"INDEPENDENT_VERIFICATION.json"),"laws_pass":rs["n_tests_meeting_bound_and_0.01"],"FIN_max_TV":rs["FIN_max_TV_rom"],"max_one_time_prediction_coordinate_change":delta,"byte_identical_lock":sha(args.replay/"RG-D04/LOCK.json")==sha(SOURCE/"RG-D04/LOCK.json"),"context":"different local numpy/scipy/mpmath versions; byte identity not required"}
    out = {
        "date":"2026-10-14","manifest":{"entries":len(manifest),"mismatches":0,"uncovered_files":uncovered,"stage_input_hashes_valid":stage_input_count},
        "source_hashes":{str(core.relative_to(ROOT)):sha(core),"AGENTS.md":sha(ROOT/"AGENTS.md"),"FIN_GENERATIVE_MEMORY_NEXT_CHAIN_20261014.md":sha(ROOT/"FIN_GENERATIVE_MEMORY_NEXT_CHAIN_20261014.md"),"campaign_manifest":sha(SOURCE/"MANIFEST.sha256"),"checker":sha(Path(__file__))},
        "plan_byte_identical":sha(ROOT/"FIN_GENERATIVE_MEMORY_NEXT_CHAIN_20261014.md")==sha(SOURCE/"source_context/FIN_GENERATIVE_MEMORY_NEXT_CHAIN_20261014.md"),
        "plan_format_normalized_identical":normalized((ROOT/"FIN_GENERATIVE_MEMORY_NEXT_CHAIN_20261014.md").read_text())==normalized((SOURCE/"source_context/FIN_GENERATIVE_MEMORY_NEXT_CHAIN_20261014.md").read_text()),
        "predecessor_intake_byte_identical":sha(ROOT/"fin_generative_review_20261014/INTAKE_20261014.md")==sha(SOURCE/"source_context/INTAKE_20261014.md"),
        "predecessor_intake_format_normalized_identical":normalized((ROOT/"fin_generative_review_20261014/INTAKE_20261014.md").read_text())==normalized((SOURCE/"source_context/INTAKE_20261014.md").read_text()),
        "freeze_hashes_valid":True,"frozen_input_hashes_valid":True,
        "Poisson_tail_public_validated_ceiling":"2.832e-52 at cutoff65,z<=4.5",
        "runtime":{"numpy":np.__version__,"scipy":scipy.__version__,"mpmath":mp.__version__},
        "independent_numeric":{"laws":len(rows),"rows":rows,"diagnostics":diagnostics,"moments":moments,"reference_method":"versioned safe-seed -> independently assembled labelled generator -> lambda4/cutoff90 uniformization; no new producer model imports"},
        "FIN_memory_moment_certificate":certificate,"source_replay":replay,
        "limits":["No historical independent custody/timestamp established by hash checks", "Source ROM expm and polynomial-norm numeric margins are assumptions, not outward certificates", "Interval moment certificate assumes mpmath.iv correctness", "No physical source, physical units, experiment, or general probability-process closure"]}
    args.output.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"manifest":out["manifest"],"laws":len(rows),"moment_certificate":certificate,"replay":None if replay is None else {k:v for k,v in replay.items() if k!="checks"}},indent=2))


if __name__ == "__main__":
    main()
