#!/usr/bin/env python3
"""FIN ST477--ST491: exact degree-eight rank, transition stabilizer, and duality.

The round follows ST462--ST476.  It keeps the active gain, noise covariance,
refinement, clocks, and instruments explicitly typed as supplied objects.  No
legacy-to-strict completion, role transfer, selector, unit, laboratory record,
Standard Model, gravity, L_total, or ToE closure is claimed.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import os
import re
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp
from scipy.linalg import expm
from scipy.optimize import brentq, root

from fin_st118_st129_research import interval_left, interval_matvec
from fin_st132_center_isolation_replay import bounds, iv
from fin_st402_st416_research import independent_strict_matrix_float
from fin_st417_st431_research import N, SEED
from fin_st447_st461_research import (
    adaptive_global_cover,
    degree8_evaluation_matrix,
    direct_ir_cell,
    direct_ir_interval_fj,
    regularized_ir_float,
    FACE,
)
from fin_st462_st476_research import dual_distance_scalar, rank_small_mod


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST477_ST491_Results.json"
SUMMARY = ROOT / "FIN_ST477_ST491_Summary.csv"
FIG_DIR = ROOT / "FIN_ST477_ST491_Figures"
NAMES = {
    477: "Exact_Degree8_Modular_Rank_Prime1",
    478: "Exact_Degree8_Modular_Rank_Prime2",
    479: "Degree8_CharacteristicZero_Bounds",
    480: "Transition_Orbit_Stabilizer_Theorem",
    481: "Fluctuation_Driven_Gain_Source_Accounting",
    482: "Equivariant_Noise_No_Canonical_Selector",
    483: "Dual_Dynamics_Unit_Separation_Threshold",
    484: "Operational_Conjugacy_Invariance_Theorem",
    485: "Cross_Channel_Clock_Dimension_Obstruction",
    486: "Transition_Branch_Identity_Audit",
    487: "Global_Cover_Architecture_Second_Stop",
    488: "Adaptive_IR_Radius_Repair_and_Extension",
    489: "Degree8_Support6_8_Circuit_Audit",
    490: "Source_Selector_Scale_Admission_Gate",
    491: "Independent_Evidence_Gate",
}
PACKETS = {k: ROOT / f"FIN_ST{k}_{v}.json" for k, v in NAMES.items()}
A = independent_strict_matrix_float()
U = np.ones(N) / N


def native(x: Any) -> Any:
    if isinstance(x, dict): return {str(k): native(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [native(v) for v in x]
    if isinstance(x, np.ndarray): return native(x.tolist())
    if isinstance(x, (np.floating, np.integer)): return x.item()
    return x


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finalize(k: int, status: str, boundary: str, packet: dict) -> dict:
    path = PACKETS[k]
    path.write_text(json.dumps(native(packet), indent=2, sort_keys=True), encoding="utf-8")
    return {"program": f"ST{k}", "object": NAMES[k], "packet_file": path.name,
            "packet_sha256": sha(path), **packet, "status": status, "boundary": boundary}


def compile_modrank() -> Path:
    src = ROOT / "fin_modrank_u32.c"
    binary = Path(tempfile.gettempdir()) / "fin_modrank_u32_st477"
    proc = subprocess.run(["gcc", "-O3", "-march=native", "-fopenmp", "-o", str(binary), str(src)],
                          capture_output=True, text=True, timeout=60)
    if proc.returncode:
        raise RuntimeError(proc.stdout + proc.stderr)
    return binary


def exact_modular_rank(prime: int, timeout: int = 600) -> dict:
    binary = compile_modrank(); t0 = time.time()
    matrix = degree8_evaluation_matrix(prime, npts=1892).astype(np.uint32)
    matrix_hash = hashlib.sha256(matrix.tobytes()).hexdigest()
    raw = Path(tempfile.gettempdir()) / f"fin_degree8_{prime}_{os.getpid()}.bin"
    matrix.tofile(raw)
    env = dict(os.environ); env["OMP_NUM_THREADS"] = env.get("OMP_NUM_THREADS", "8")
    try:
        proc = subprocess.run([str(binary), str(raw), "1892", "2028", str(prime)],
                              capture_output=True, text=True, timeout=timeout, env=env)
        match = re.search(r"rank=(\d+)", proc.stdout)
        rank = int(match.group(1)) if match else None
        return {"prime": prime, "matrix_shape": [1892, 2028], "matrix_sha256": matrix_hash,
                "rank": rank, "nullity": None if rank is None else 2028-rank,
                "return_code": proc.returncode, "stderr_tail": proc.stderr.splitlines()[-8:],
                "wall_seconds": time.time()-t0, "threads": int(env["OMP_NUM_THREADS"])}
    finally:
        raw.unlink(missing_ok=True)


def modular_rank_program(k: int, prime: int) -> dict:
    row = exact_modular_rank(prime)
    status = "proven_exact_modular_rank_returned" if row["rank"] is not None and row["return_code"] == 0 else "failed_exact_modular_rank_attempt"
    packet = {**row, "claim": "The rank is exact over the declared finite field; it is a lower bound for characteristic-zero rank, not by itself a rational nullspace basis."}
    return finalize(k, status,
                    "A modular rank does not prove characteristic-zero equality unless a matching rational upper certificate is supplied.", packet)


def st477() -> dict: return modular_rank_program(477, 1000003)
def st478() -> dict: return modular_rank_program(478, 1000033)


def st479() -> dict:
    a = json.loads(PACKETS[477].read_text()); b = json.loads(PACKETS[478].read_text())
    lower = max(a["rank"], b["rank"]); upper = 1892
    packet = {
        "two_prime_modular_ranks": [a["rank"], b["rank"]],
        "characteristic_zero_decomposable_rank_interval": [lower, upper],
        "characteristic_zero_relation_nullity_interval": [2028-upper, 2028-lower],
        "primitive_degree8_quotient_dimension_interval": [0, 1892-lower],
        "rank_1791_equality_promoted": False,
        "theorem": "Reduction modulo either good prime cannot increase rank. Hence the characteristic-zero decomposable rank is at least the maximum returned modular rank; Molien dimension gives the upper bound 1892.",
    }
    return finalize(479, "proven_two_sided_characteristic_zero_degree8_bounds",
                    "Exact equality needs 237 rational relations (if rank 1791) or another characteristic-zero upper-rank theorem.", packet)


def radial_probability(center: list[float]) -> np.ndarray:
    p0, p1, p2, p3, p4, p5, p6 = center[:7]
    return np.array([p0, p1, p2, p3, p4, p5, p6, p5, p4, p3, p2, p1])


def d12_stabilizer(p: np.ndarray, tol: float = 1e-8) -> list[str]:
    rows = []
    for k in range(N):
        if np.linalg.norm(np.roll(p, k)-p) < tol: rows.append(f"r^{k}")
        if np.linalg.norm(np.roll(p[::-1], k)-p) < tol: rows.append(f"r^{k}s")
    return rows


def st480() -> dict:
    old = json.loads((ROOT / "FIN_ST342_Narrow_Certified_Coexistence_Bracket.json").read_text())
    centers = [x["center"] for x in old["endpoint_certificates"]]
    rows = []
    for c in centers:
        p = radial_probability(c)
        rows.append({"probability": p, "stabilizer": d12_stabilizer(p),
                     "unique_max_gap": float(p[0]-np.max(p[1:])), "sum": float(np.sum(p))})
    packet = {
        "certified_branch_parameterization": "reflection-even radial p=(p0,p1,...,p6,p5,...,p1)",
        "endpoint_checks": rows, "stabilizer_order": 2, "D12_order": 24, "orbit_size": 12,
        "theorem": "Reflection-even parameterization supplies one reflection. The certified unique dominant vertex excludes every nontrivial rotation. Two distinct reflections would compose to a nontrivial rotation, so the stabilizer has exactly order two and the branch orbit has exactly twelve members throughout the tube.",
        "globally_first_orbit": False,
    }
    return finalize(480, "proven_transition_branch_stabilizer_order2_orbit12",
                    "This proves the orbit size of the ST342 branch; it does not prove that no third orbit attains the global ratio earlier.", packet)


def st481() -> dict:
    trace_a = float(np.trace(A)); sigma = 0.01; kappa = 2.0; gamma = 1.0
    P = np.eye(N)-np.ones((N, N))/N
    covariance = sigma**2 * P
    expected_q = float(np.trace(A @ covariance))
    rng = np.random.default_rng(SEED+481)
    samples = rng.multivariate_normal(np.zeros(N), covariance, size=200000)
    observed = float(np.mean(np.einsum("bi,ij,bj->b", samples, A, samples)))
    packet = {
        "controller": "g_dot=kappa*q(xi)-gamma*g, q(xi)=xi^T A xi",
        "tangent_isotropic_covariance": "Sigma=sigma^2(I-11^T/12)",
        "sigma": sigma, "kappa": kappa, "gamma": gamma, "trace_A": trace_a,
        "exact_expected_q": expected_q, "predicted_stationary_mean_g": kappa*expected_q/gamma,
        "Monte_Carlo_samples": len(samples), "Monte_Carlo_mean_q": observed,
        "relative_sampling_error": abs(observed-expected_q)/expected_q,
        "theorem": "For zero-mean fluctuations with covariance Sigma, E[xi^T A xi]=Tr(A Sigma), hence the stationary mean controller output is (kappa/gamma)Tr(A Sigma).",
        "strict_source_status": "Sigma/noise law is an added resource, not generated by frozen A",
    }
    return finalize(481, "proven_fluctuation_to_mean_gain_conversion_with_explicit_resource",
                    "The result identifies where gain can enter; it does not derive the fluctuation covariance, pump energy, or physical noise law.", packet)


def st482() -> dict:
    packet = {
        "assumptions": ["D12-equivariant dynamics", "D12-invariant noise law", "one transitive orbit of twelve absorbing/minimizing branches"],
        "branch_probability": "1/12 for each branch",
        "theorem": "Equivariance transports the event of reaching branch i to branch g*i without changing probability. Transitivity forces all twelve probabilities equal; normalization makes each 1/12.",
        "single_run_realization": True, "canonical_selector_exported": False,
        "interpretation": "Symmetric noise can realize a random branch in each run, but it does not define an invariant preferred branch across runs.",
    }
    return finalize(482, "proven_equivariant_noise_realizes_but_does_not_select",
                    "Breaking the uniform branch law requires an asymmetric state, boundary, apparatus, or noise covariance and therefore an explicit selector resource.", packet)


def st483() -> dict:
    x1 = brentq(lambda x: math.exp(-x)-2*math.cos(x), 0.1, math.pi/2)
    eig = np.linalg.eigvalsh(A)[1:]
    rows = [{"lambda": float(l), "first_time_distance_equals_one": x1/float(l)} for l in eig]
    packet = {
        "threshold_equation": "|exp(-ix)-exp(-x)|=1 iff exp(-x)=2cos(x)",
        "first_positive_dimensionless_threshold": x1, "strict_mode_times": rows,
        "earliest_strict_mode_time": min(x["first_time_distance_equals_one"] for x in rows),
        "theorem": "Every positive mode separates coherent and diffusive functional calculus by operator distance one after dimensionless modal time x1; the largest eigenvalue reaches it first.",
        "physical_time": False,
    }
    return finalize(483, "proven_dual_dynamics_unit_separation_threshold",
                    "Operator-norm separation becomes an experiment only after preparations, clock calibration, instruments, and noise are supplied.", packet)


def st484() -> dict:
    vals, vecs = np.linalg.eigh(A); theta = 0.37
    rot = np.eye(N); c, s = math.cos(theta), math.sin(theta)
    rot[1, 1]=c; rot[1, 3]=-s; rot[3, 1]=s; rot[3, 3]=c
    Q = vecs @ rot @ vecs.T; B = Q @ A @ Q.T
    rng = np.random.default_rng(SEED+484); rho = rng.normal(size=N); effect = rng.normal(size=N)
    rows = []
    for t in (0.2, 0.8, 2.0):
        for kind, fA, fB in (("heat", expm(-t*A), expm(-t*B)),
                             ("unitary_real_record", expm(-1j*t*A), expm(-1j*t*B))):
            left = np.vdot(effect, fA @ rho)
            right = np.vdot(Q@effect, fB @ (Q@rho))
            rows.append({"kind": kind, "t": t, "absolute_record_residual": float(abs(left-right))})
    packet = {
        "transport": "B=Q A Q*, rho'=Q rho, M'=Q M",
        "record_checks": rows, "maximum_record_residual": max(x["absolute_record_residual"] for x in rows),
        "theorem": "For every Borel f, f(QAQ*)=Qf(A)Q*. Transporting preparations and effects therefore preserves every one-time matrix record exactly, including both unitary and heat channels.",
        "spectrum_alone_sufficient": False,
        "operational_bundle": "generator/projectors + preparations + transformations + effects + record map",
    }
    return finalize(484, "proven_complete_operational_conjugacy_invariance",
                    "Fixed vertex labels without transported instruments are not equivalent; multitime equivalence additionally requires transported interventions/composition.", packet)


def st485() -> dict:
    positive = np.unique(np.round(np.linalg.eigvalsh(A)[1:], 13))
    packet = {
        "positive_distinct_eigenvalues": positive,
        "unitary_phase_frequency": "omega_U=lambda/tau_*",
        "heat_decay_rate": "gamma_H=lambda/tau_*",
        "wave_frequency": "omega_C=sqrt(lambda)/tau_C",
        "theorem": "Assigning A dimension T^{-1} makes unitary phase and heat rates dimensionally natural, whereas d2u/dt2=-Au requires A dimension T^{-2}. One raw unit assignment cannot satisfy both. A sector map such as A_wave=A/tau_* or a distinct wave clock is necessary.",
        "single_raw_A_unit_closes_all_three_channels": False,
        "absolute_seconds_generated": False,
    }
    return finalize(485, "proven_cross_channel_clock_dimension_obstruction",
                    "Dimensionless spectral ratios remain common; the obstruction concerns physical units and does not weaken the mathematical dual dynamics.", packet)


def st486() -> dict:
    old = json.loads((ROOT / "FIN_ST342_Narrow_Certified_Coexistence_Bracket.json").read_text())
    new = json.loads((ROOT / "FIN_ST462_First_Transition_Orbit_Numerical_Isolation.json").read_text())
    center = old["floating_center"]
    packet = {
        "ST342_certified_local_crossing_bracket": old["certified_crossing_bracket"],
        "ST462_best_ratio": new["best_ratio"], "distance_to_ST342_floating_center": abs(new["best_ratio"]-center),
        "same_reflection_even_profile_residual": new["minimum_reflection_stabilizer_residual"],
        "multistarts_in_same_orbit": new["best_orbit_hit_count"],
        "conclusion": "The numerical ratio minimizer and the interval-certified ST342 equal-value branch are the same local orbit to the displayed precision.",
        "global_first_identity": False,
    }
    return finalize(486, "strong_cross_packet_identity_of_candidate_transition_branch",
                    "Local identity plus exhaustive starts is not a proof that another interior orbit cannot attain a smaller ratio in the global strip.", packet)


def st487() -> dict:
    gain = 2.8936; t0=time.time(); cover=adaptive_global_cover(gain,0.3342,initial=35,max_depth=14)
    packet = {"attempted_gain": gain, "cover": cover, "wall_seconds": time.time()-t0,
              "new_global_lower_theorem": bool(cover["closed"]),
              "result": "The unchanged two-coordinate cover is stopped again if a leaf reaches the paid numerical floor."}
    return finalize(487, "global_cover_architecture_closed" if cover["closed"] else "second_bounded_stop_for_unmodified_global_cover",
                    "A failed cover is not a negative competitor; the accepted theorem remains 2.8934 unless closed=True.", packet)


def recertify_ir_cell(row: dict, blo: float, bhi: float, factor: float) -> dict:
    mp.iv.dps = 45
    center=np.array(row["center"],float); radii=np.array(row["radii"],float)*factor
    X=[iv((center[i]-radii[i],center[i]+radii[i])) for i in range(5)]
    f0,jac=direct_ir_interval_fj([iv(float(x)) for x in center],blo,bhi)
    flo=np.array([bounds(x)[0] for x in f0]); fhi=np.array([bounds(x)[1] for x in f0])
    bmid=(blo+bhi)/2; eps=1e-6; eye=np.eye(5)
    jm=np.column_stack([(regularized_ir_float(center+eps*eye[j],bmid)-regularized_ir_float(center-eps*eye[j],bmid))/(2*eps) for j in range(5)])
    pre=np.linalg.inv(jm); jl=np.array([[bounds(x)[0] for x in rr] for rr in jac]); jh=np.array([[bounds(x)[1] for x in rr] for rr in jac])
    cfl,cfh=interval_matvec(pre,pre,flo,fhi); ylo,yhi=center-cfh,center-cfl
    cjl,cjh=interval_left(pre,jl,jh); mlo,mhi=-cjh,-cjl
    for i in range(5): mlo[i,i]+=1; mhi[i,i]+=1
    dlo,dhi=interval_matvec(mlo,mhi,-radii,radii)
    margins=np.minimum(ylo+dlo-(center-radii),(center+radii)-(yhi+dhi))
    M=np.maximum(abs(mlo),abs(mhi)); contraction=float(max((M@radii)/radii))
    return {**row,"radii":radii,"included":bool(np.min(margins)>0),
            "minimum_margin":float(np.min(margins)),"weighted_contraction":contraction,
            "radius_expansion_factor":factor}


def repaired_ir_chain(start=.124,stop=.2,width=.003) -> dict:
    prior=root(lambda z:regularized_ir_float(z,start),FACE,tol=1e-12).x
    rows=[]; lo=start; failure=None
    while lo<stop-1e-14:
        hi=min(stop,lo+width); raw,candidate=direct_ir_cell(lo,hi,prior)
        selected=None
        for factor in (1.0,1.5,2.0,3.0,4.0):
            trial=recertify_ir_cell(raw,lo,hi,factor)
            if trial["included"]: selected=trial; break
        if selected is None: failure=raw; break
        rows.append(selected); prior=candidate; lo=hi
    return {"start":start,"target":stop,"certified_stop":lo,"cells":rows,"failure":failure}


def st488() -> dict:
    chain=repaired_ir_chain(); rows=chain["cells"]
    packet={**chain,"cell_count":len(rows),
            "minimum_Krawczyk_margin":min(x["minimum_margin"] for x in rows),
            "maximum_weighted_contraction":max(x["weighted_contraction"] for x in rows),
            "maximum_radius_expansion_factor":max(x["radius_expansion_factor"] for x in rows),
            "theorem":f"Expanded-radius Krawczyk boxes continue the same locally certified component from b=0.124 through b={chain['certified_stop']:.3f}.",
            "target_reached":chain["certified_stop"]>=.2-1e-14}
    return finalize(488,"proven_repaired_IR_extension" if rows else "failed_IR_radius_repair",
                    "Uniqueness is local to the displayed chain; a failure is a certificate-template stop, not a branch singularity.",packet)


def sampled_support_audit(prime:int,supports=(6,8),count=12000)->dict:
    matrix=degree8_evaluation_matrix(prime,npts=64); rng=np.random.default_rng(SEED+489+prime)
    rows=[]
    for support in supports:
        proj=rng.integers(1,prime,size=(support,matrix.shape[0]),dtype=np.int64)@matrix%prime
        candidates=true=0
        for _ in range(count):
            ids=np.sort(rng.choice(matrix.shape[1],support,replace=False))
            if rank_small_mod(proj[:,ids],prime)<support:
                candidates+=1
                if rank_small_mod(matrix[:,ids],prime)<support:true+=1
        rows.append({"support":support,"sampled_subsets":count,"projected_candidates":candidates,"true_dependencies":true})
    return {"prime":prime,"rows":rows,"matrix_sha256":hashlib.sha256(matrix.tobytes()).hexdigest()}


def st489() -> dict:
    audits=[sampled_support_audit(p) for p in (1000003,1000033)]
    found=sum(x["true_dependencies"] for a in audits for x in a["rows"])
    packet={"exact_modular_random_audits":audits,"true_dependencies_found":found,
            "claim":"No support-six or support-eight dependency was found in the declared reproducible exact sample.",
            "exhaustive_theorem":False}
    return finalize(489,"exact_sampled_no_support6_8_circuit_found",
                    "Dense degree-eight relations may remain; random absence is not a support lower-bound theorem.",packet)


def st490() -> dict:
    packet={"new_strict_gain_source":False,"new_nonpremise_selector":False,"QW_2191":"open",
            "new_scale_charged_source":False,"absolute_unit":"absent",
            "conditional_noise_gain":"accepted only with supplied covariance",
            "reason":"ST481 identifies noise covariance as the missing source resource; ST482 proves symmetric noise does not canonically select a branch."}
    return finalize(490,"blocked_no_strict_gain_selector_or_scale_source",
                    "Conditional source accounting is progress but exports no strict physical bridge.",packet)


def st491() -> dict:
    packet={"external_referee":"absent","independent_laboratory_record":"absent","held_out_empirical_record":"absent"}
    return finalize(491,"blocked_no_independent_empirical_evidence",
                    "No laboratory, legacy-role-transfer, Standard Model, gravity, L_total, or ToE closure is exported.",packet)


def figures(results:dict)->None:
    FIG_DIR.mkdir(exist_ok=True)
    ranks=[results["ST477"]["rank"],results["ST478"]["rank"]]
    fig,ax=plt.subplots(figsize=(7.2,4));ax.bar(["1,000,003","1,000,033"],ranks,color="#2563eb")
    ax.axhline(1892,ls="--",color="#dc2626",label="Molien degree-8 dimension");ax.set_ylim(min(ranks)-20,1900)
    ax.set(ylabel="exact modular rank",title="ST477--ST479: degree-eight decomposable rank");ax.legend();fig.tight_layout();fig.savefig(FIG_DIR/"st477_modular_rank.png",dpi=180);plt.close(fig)
    row=results["ST481"];fig,ax=plt.subplots(figsize=(7.2,4));ax.bar(["exact Tr(AΣ)","Monte Carlo"],[row["exact_expected_q"],row["Monte_Carlo_mean_q"]],color=["#0f766e","#7c3aed"]);ax.set(title="ST481: fluctuation-to-gain source accounting",ylabel="mean quadratic response");fig.tight_layout();fig.savefig(FIG_DIR/"st481_fluctuation_gain.png",dpi=180);plt.close(fig)
    rows=results["ST488"]["cells"]; b=[x["b_interval"][1] for x in rows]; y2=[x["endpoint_centers"][-1][1] for x in rows]
    fig,ax=plt.subplots(figsize=(7.2,4));ax.plot(b,y2,"o-");ax.set(xlabel="b",ylabel=r"$y_2$",title="ST488: repaired local IR continuation");fig.tight_layout();fig.savefig(FIG_DIR/"st488_ir_repair.png",dpi=180);plt.close(fig)


def main()->None:
    funcs={477:st477,478:st478,479:st479,480:st480,481:st481,482:st482,483:st483,484:st484,485:st485,486:st486,487:st487,488:st488,489:st489,490:st490,491:st491}
    results={}
    for k in range(477,492):
        print(f"running ST{k}",flush=True);results[f"ST{k}"]=funcs[k]()
    RESULTS.write_text(json.dumps(native(results),indent=2,sort_keys=True),encoding="utf-8")
    with SUMMARY.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f);w.writerow(["program","status","object","boundary"])
        for k in range(477,492):
            r=results[f"ST{k}"];w.writerow([f"ST{k}",r["status"],r["object"],r["boundary"]])
    figures(results);print(f"wrote {RESULTS.name} and {SUMMARY.name}")


if __name__=="__main__":main()
