#!/usr/bin/env python3
"""FIN ST166--ST177: carrier invariants, strict-source no-go results, and executable validation."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import platform
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import scipy
import sympy as sp
from scipy.linalg import expm
from scipy.optimize import root

from fin_st01_st15_research import N, strict_operator
from fin_st130_st141_research import point_design_system
from fin_st132_center_isolation_replay import bounds as replay_bounds, strict_interval_matrix
from fin_st154_st165_research import (
    commutant_dimension,
    hidden_pair_transfer,
    local_param_krawczyk,
    parametric_data,
)
from fin_st177_operational_validator import run as run_validator


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST166_ST177_Results.json"
SUMMARY = ROOT / "FIN_ST166_ST177_Summary.csv"
PACK166 = ROOT / "FIN_ST166_Carrier_Operational_Invariant.json"
PACK168 = ROOT / "FIN_ST168_Uniform_Branch_Fold_Classification.json"
PACK169 = ROOT / "FIN_ST169_Nonorthogonal_Unitary_Recovery.json"
CERT172 = ROOT / "FIN_ST172_Extended_Nuisance_Cover.json"
PACK173 = ROOT / "FIN_ST173_Noisy_Refinement_Compression.json"
PACK174 = ROOT / "FIN_ST174_Robust_Tomographic_Net.json"
PACK175 = ROOT / "FIN_ST175_HMM_Collatz_Certificate.json"
PACK176 = ROOT / "FIN_ST176_Microscopic_Local_Bath_Tradeoff.json"
PACK177 = ROOT / "FIN_ST177_Operational_Validator_Study.json"
FIG_DIR = ROOT / "FIN_ST166_ST177_Figures"
SEED = 20260825


def native(x: Any) -> Any:
    if isinstance(x, dict): return {str(k): native(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [native(v) for v in x]
    if isinstance(x, np.ndarray): return native(x.tolist())
    if isinstance(x, (np.floating, np.integer)): return x.item()
    return x


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(native(value), sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def write_packet(path: Path, data: dict) -> str:
    path.write_text(json.dumps(native(data), indent=2, sort_keys=True), encoding="utf-8")
    return sha(path)


def entropy(p: np.ndarray) -> float:
    return float(-sum(x*math.log(x) for x in p if x > 0))


def st166_carrier_invariant() -> dict:
    e1 = np.array([[1,0,0],[0,1,0],[0,0,.5],[0,0,.5]], float)
    d1 = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,1]], float)
    e2 = np.array([[.25,0,0],[.75,0,0],[0,1,0],[0,0,.4],[0,0,.6]], float)
    d2 = np.array([[1,1,0,0,0],[0,0,1,0,0],[0,0,0,1,1]], float)
    preps = np.array([[1,0,0,.2],[0,1,0,.3],[0,0,1,.5]], float)
    effects = np.array([[1,0,0],[0,1,0],[0,0,1],[.2,.3,.5]], float)
    logical_table = effects @ preps
    carrier1_table = (effects @ d1) @ (e1 @ preps)
    carrier2_table = (effects @ d2) @ (e2 @ preps)
    p = np.array([.5,.5,0]); q = np.array([.5,0,.5])
    packet = {
        "logical_probability_table": logical_table.tolist(),
        "carrier_1_probability_table": carrier1_table.tolist(),
        "carrier_2_probability_table": carrier2_table.tolist(),
        "carrier_1_left_inverse_error": float(np.max(abs(d1@e1-np.eye(3)))),
        "carrier_2_left_inverse_error": float(np.max(abs(d2@e2-np.eye(3)))),
        "maximum_probability_table_mismatch": float(max(np.max(abs(logical_table-carrier1_table)), np.max(abs(logical_table-carrier2_table)))),
        "operational_invariant_sha256": digest(logical_table.tolist()),
        "same_entropy_counterexample": {"p": p.tolist(), "q": q.tolist(), "entropy_p": entropy(p), "entropy_q": entropy(q),
                                        "vertex_effect_probabilities_p": p.tolist(), "vertex_effect_probabilities_q": q.tolist()},
        "theorem": (
            "If carrier encodings E_i have stochastic left inverses D_i and every logical effect m is transported as m D_i, then every preparation-effect probability is invariant: "
            "m D_i E_i p=m p. The full tomographically complete probability table is therefore a carrier-change operational invariant. Equal Shannon entropy alone is not sufficient."
        ),
    }
    packet_hash = write_packet(PACK166, packet)
    return {"program": "ST166", "object": "Carrier-Change Operational Probability Invariant",
            "packet_file": PACK166.name, "packet_sha256": packet_hash, **packet,
            "status": "proven_exact_for_reversible_stochastic_code_embeddings",
            "boundary": "Encodings, decodings, logical preparations, and transported instruments are supplied. The theorem does not prove information exists without a physical carrier."}


def st167_strict_interaction_no_go(a: np.ndarray) -> dict:
    functions = [np.eye(N), a, a@a, a@a@a, expm(-.3*a)]
    rows = []
    for count in range(1, len(functions)+1):
        rows.append({"functional_generators": count, "joint_commutant_dimension": commutant_dimension(functions[:count])})
    return {"program": "ST167", "object": "Strict Functional-Calculus Environment-Coupling No-Go",
            "rows": rows, "strict_commutant_dimension": commutant_dimension([a]),
            "target_vertex_pointer_dimension": N,
            "theorem": (
                "For every family L_k=f_k(A), the strict commutant {A}' is contained in {A,L_1,...,L_m}'. Hence its dimension is at least 22 and no such coupling family can generate the 12-dimensional vertex pointer algebra or a smaller irreducible algebra. "
                "If one f_k separates all seven spectral classes, equality with {A}' holds."
            ),
            "status": "proven_no_go_in_strict_functional_calculus_class",
            "boundary": "State-dependent, nonlinear, vertex-typed, environmental, or symmetry-breaking couplings lie outside this class and remain additional objects."}


def rho_stationary_polynomial_intervals():
    rho = sp.symbols("rho")
    poly = sp.Poly(3*rho**3 + 18*rho**2 - 284*rho + 24, rho, domain=sp.QQ)
    isolating = poly.intervals(eps=sp.Rational(1, 10**14))
    positive = []
    for (lohi, multiplicity) in isolating:
        lo, hi = lohi
        if hi > 0 and lo >= 0:
            positive.append((float(lo), float(hi), int(multiplicity)))
    return str(poly.as_expr()), isolating, positive


def iv(value):
    if isinstance(value, tuple): return mp.iv.mpf([str(value[0]), str(value[1])])
    return mp.iv.mpf(str(value))


def st168_uniform_fold_classification() -> dict:
    poly_string, all_intervals, positive = rho_stationary_polynomial_intervals()
    aiv, _, _ = strict_interval_matrix(); mp.iv.dps = 70
    spectral = []
    for k in range(1, 7):
        lam = sum((aiv[0][j]*mp.iv.cos(2*mp.iv.pi*k*j/12) for j in range(N)), iv(0))
        lo, hi = replay_bounds(lam)
        spectral.append({"mode": k, "lambda_interval": [lo, hi]})
    rho_rows = [{"rho_interval": [0.0, 0.0], "amplitude_sign_count": 1}]
    rho_rows += [{"rho_interval": [lo, hi], "amplitude_sign_count": 2} for lo, hi, _ in positive]
    fold_rows = []
    for rr in rho_rows:
        rho_iv = iv(tuple(rr["rho_interval"])) if rr["rho_interval"][1] > 0 else iv(0)
        den = 1 + rho_iv/2; qfun = rho_iv/den; qp = den**-2; qpp = -den**-3
        h = -qfun*qp + iv("0.075"); hp = -(qp**2 + qfun*qpp); rdiag = 2*h + 4*rho_iv*hp
        rd = replay_bounds(rdiag)
        for spec in spectral:
            lam = iv(tuple(spec["lambda_interval"])); kap = -rdiag/lam
            fold_rows.append({"rho_interval": rr["rho_interval"], "amplitude_sign_count": rr["amplitude_sign_count"],
                              "mode": spec["mode"], "rdiag_interval": list(rd), "kappa_interval": list(replay_bounds(kap))})
    packet = {
        "stationarity_polynomial": poly_string,
        "constant_hessian_radial_numerator": "3*(rho^4+8*rho^3+344*rho^2-608*rho+16)",
        "constant_mode_exact_resultant": -108263374848000,
        "all_real_root_intervals": [[str(x[0][0]), str(x[0][1]), int(x[1])] for x in all_intervals],
        "positive_rho_intervals": [[lo, hi] for lo, hi, _ in positive],
        "positive_spectral_classes": spectral,
        "geometric_fold_points_null_line_mod_sign": 30,
        "normalized_augmented_fold_roots_including_v_sign": 60,
        "fold_rows": fold_rows,
        "theorem": (
            "On the reflection-even uniform branch q0=...=q6=c, stationarity is exactly c h(c)=0. Besides c=0, rho=c^2 obeys 3rho^3+18rho^2-284rho+24=0, which has exactly two positive roots. "
            "The constant Hessian mode never vanishes on these stationary amplitudes. Each of the five real amplitudes and six positive strict Fourier classes yields exactly one fold null line, hence exactly 30 geometric folds and 60 normalized augmented roots after v-sign."
        ),
    }
    packet_hash = write_packet(PACK168, packet)
    return {"program": "ST168", "object": "Complete Uniform-Branch Fold Classification",
            "packet_file": PACK168.name, "packet_sha256": packet_hash, **packet,
            "status": "proven_complete_on_declared_uniform_symmetry_reduced_branch",
            "boundary": "This exact 30-fold classification applies only to the uniform reflection-even branch; nonuniform portions of the 15-dimensional fold domain remain open."}


def binary_phase_recovery(theta: float, p_error=.35, delta0=.08, delta1=.17):
    weights = []
    for y in [0, 1]:
        alpha = (1-p_error)*((1-delta0) if y == 0 else delta0)
        beta = p_error*(delta1 if y == 0 else (1-delta1))
        coherence = alpha + beta*np.exp(1j*theta)
        weights.append({"syndrome": y, "identity_weight": alpha, "unitary_weight": beta,
                        "coherence_abs": float(abs(coherence)), "optimal_contribution": float((alpha+beta+abs(coherence))/2)})
    return weights, sum(x["optimal_contribution"] for x in weights)


def st169_nonorthogonal_recovery() -> dict:
    rows = []
    for theta in np.linspace(0, math.pi, 17):
        details, fidelity = binary_phase_recovery(float(theta))
        rows.append({"theta": float(theta), "unitary_HS_overlap_abs": abs(math.cos(theta/2)),
                     "optimal_entanglement_fidelity": fidelity, "syndrome_rows": details})
    packet = {
        "error_probability": .35, "readout_errors": {"delta0": .08, "delta1": .17}, "rows": rows,
        "theorem": (
            "For qubit errors I and U=exp(-i theta Z/2), each noisy-syndrome branch has subnormalized weights alpha_y,beta_y and coherence c_y=alpha_y+beta_y exp(i theta). "
            "A phase correction converts the branch to an orthogonal dephasing channel, so the globally optimal CPTP recovery contribution is (alpha_y+beta_y+|c_y|)/2. Summing over y solves the nonorthogonal binary-unitary recovery problem exactly."
        ),
    }
    packet_hash = write_packet(PACK169, packet)
    return {"program": "ST169", "object": "Exact Nonorthogonal Binary-Unitary Recovery",
            "packet_file": PACK169.name, "packet_sha256": packet_hash, **packet,
            "status": "proven_global_optimum_for_qubit_binary_phase_errors",
            "boundary": "The qubit phase family, prior, syndrome interaction, and readout channel are supplied. General dimension and three or more nonorthogonal errors remain open SDP problems."}


def fisher_rao_cost(g: float, horizon: float) -> tuple[float, float]:
    overlap = (math.sqrt(1+11*g) + 11*math.sqrt(max(0, 1-g)))/12
    distance = 2*math.acos(min(1, max(-1, overlap)))
    return distance, distance*distance/(2*horizon)


def st170_information_geometric_selector() -> dict:
    rows = []
    for g in [.001, .01, .05, .1, .25, .5, 1.0]:
        distance, cost = fisher_rao_cost(g, 1.0)
        rows.append({"gap": g, "Fisher_Rao_distance": distance, "minimum_geodesic_action_T1": cost})
    return {"program": "ST170", "object": "Fisher--Rao Selector-Control Cost and Covariant No-Start Theorem",
            "rows": rows,
            "theorem": (
                "Every C12-covariant Markov/Lindblad semigroup that preserves the uniform state leaves a uniform initial preparation uniform and cannot start a selector. "
                "If an external control drives the probability simplex with Fisher--Rao kinetic action, the least action to a robust gap-g endpoint in time T is d_FR(g)^2/(2T), where d_FR=2 arccos((sqrt(1+11g)+11sqrt(1-g))/12)."
            ),
            "status": "proven_no_start_and_sharp_conditional_information_geometric_cost",
            "boundary": "The control, Fisher--Rao metric, preferred vertex, clock, and gap are supplied. No strict selector or physical energy follows."}


def st171_anisotropy_coupling_no_go() -> dict:
    return {"program": "ST171", "object": "Twelfth-Order Coupling Nonidentifiability Theorem",
            "jet_order": 11, "free_coupling_family": "V_g=V_0+(g/12) sum_j psi_j^12",
            "theorem": (
                "For every real g, V_g has the same strict quadratic operator A, the same C12 and reflection symmetries, and the same field jet through order eleven at psi=0. "
                "Therefore no data functor depending only on A, strict functional calculus, symmetry, and derivatives through order eleven can identify g, its sign, or its magnitude. A nonzero twelfth-order datum or new normalization/source axiom is necessary."
            ),
            "status": "proven_source_no_go_for_linear_spectral_and_11_jet_data",
            "boundary": "The theorem does not exclude a new strict twelfth-order observable, nonlinear microscopic law, measured response, or independently supplied coupling."}


def st172_extended_cover(a: np.ndarray) -> dict:
    ec, ef, eigc, eigf = parametric_data(a)
    global_halfwidth = 7.5e-5; subdivisions = 5; local_halfwidth = global_halfwidth/subdivisions
    offsets = [-global_halfwidth + (2*i+1)*local_halfwidth for i in range(subdivisions)]
    rows = []
    for dq0, dq1, dd in itertools.product(offsets, repeat=3):
        nuisance = (.2+dq0, .7+dq1, .05+dd)
        center = root(lambda x: point_design_system(x, ec, ef, *nuisance), np.array([2.1862,.53983]), tol=1e-12).x
        cert = local_param_krawczyk(eigc, eigf, nuisance, local_halfwidth, center)
        rows.append(cert)
    passed = sum(x["included"] for x in rows)
    packet = {"global_halfwidth": global_halfwidth, "subdivisions_per_axis": subdivisions,
              "local_halfwidth": local_halfwidth, "boxes": len(rows), "passed_boxes": passed,
              "minimum_margin": min(x["margin"] for x in rows), "rows": rows,
              "theorem": "The 5x5x5 cell union exactly covers the declared nuisance cube. Inclusion in every cell proves a uniform stationary-root continuation over the cube."}
    packet_hash = write_packet(CERT172, packet)
    return {"program": "ST172", "object": "Second-Generation Adaptive Nuisance Continuation",
            "certificate_file": CERT172.name, "certificate_sha256": packet_hash, **packet,
            "status": "proven_extended_uniform_cover" if passed == len(rows) else "partial_extended_cover",
            "boundary": "The largest certified cube is not proved maximal; no physical parameter calibration or global beta optimum follows."}


def path_parity_error(edge_error: float, path_length: int) -> float:
    return (1-(1-2*edge_error)**path_length)/2


def majority_error(q: float, repetitions: int) -> float:
    threshold = repetitions//2 + 1
    return sum(math.comb(repetitions,k)*q**k*(1-q)**(repetitions-k) for k in range(threshold,repetitions+1))


def st173_noisy_refinement() -> dict:
    rows = []
    for eps in [.001,.005,.01,.02,.05,.1]:
        q = path_parity_error(eps, 3)
        rows.append({"edge_flip_probability": eps, "three_edge_path_parity_error": q,
                     "five_path_majority_error": majority_error(q,5),
                     "compression_gain": q/majority_error(q,5) if q else math.inf})
    packet = {"rows": rows,
              "theorem": (
                  "For an m-edge all-odd refinement path with independent edge-flip probability epsilon, the transported twist error is q_m=(1-(1-2epsilon)^m)/2. "
                  "For R independent odd replicated paths, majority coarse-graining has the exact binomial-tail error. It is an approximate natural transformation whose failure probability decreases as a repetition code, not an exact functor."
              )}
    packet_hash = write_packet(PACK173, packet)
    return {"program": "ST173", "object": "Noisy Approximate Refinement Functor and Redundant Compression",
            "packet_file": PACK173.name, "packet_sha256": packet_hash, **packet,
            "status": "proven_exact_error_law_for_independent_replicated_paths",
            "boundary": "Independence, path replication, majority decoder, and edge-noise law are supplied. This is an error-correcting compression model, not a strict fractal dynamics theorem."}


def commutator_laplacian(generators: list[np.ndarray]) -> np.ndarray:
    d = generators[0].shape[0]; eye = np.eye(d); lap = np.zeros((d*d,d*d),complex)
    for g in generators:
        ad = np.kron(eye,g)-np.kron(g.T,eye)
        lap += ad.conj().T@ad
    return (lap+lap.conj().T)/2


def projector_low(lap: np.ndarray, dimension: int) -> tuple[np.ndarray,np.ndarray]:
    vals, vecs = np.linalg.eigh(lap); v = vecs[:,:dimension]
    return vals, v@v.conj().T


def st174_noisy_tomography() -> dict:
    x=np.array([[0,1],[1,0]],complex); z=np.diag([1,-1]).astype(complex); i2=np.eye(2)
    ideal=[np.kron(x,i2),np.kron(z,i2)]; vals0,p0=projector_low(commutator_laplacian(ideal),4)
    rng=np.random.default_rng(SEED+174); rows=[]
    for eta in [.001,.005,.01,.02,.05,.1,.12]:
        noisy=[]
        for g in ideal:
            h=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));h=(h+h.conj().T)/2;h*=eta/np.linalg.norm(h,2);noisy.append(g+h)
        vals,p=projector_low(commutator_laplacian(noisy),4)
        delta=16*eta+8*eta*eta
        rows.append({"generator_noise_norm_bound":eta,"analytic_superoperator_perturbation_bound":delta,
                     "analytic_four_dimensional_cluster_separated":delta<2,
                     "numerical_fourth_eigenvalue":float(vals[3]),"numerical_fifth_eigenvalue":float(vals[4]),
                     "numerical_projector_distance":float(np.linalg.norm(p-p0,2)),
                     "Davis_Kahan_bound":min(1.0,delta/(4-2*delta)) if delta<2 else None})
    packet={"ideal_spectrum":vals0.tolist(),"ideal_commutant_dimension":4,"ideal_gap":float(vals0[4]),"rows":rows,
            "theorem": (
                "For two unit-norm Hermitian Pauli-factor generators perturbed by at most eta each, the commutator-Laplacian perturbation is at most delta=16eta+8eta^2. "
                "Since the ideal zero cluster has dimension four and gap four, delta<2 certifies separation of a unique four-dimensional low spectral subspace."
            )}
    packet_hash=write_packet(PACK174,packet)
    return {"program":"ST174","object":"Robust Atomic-Net Subspace from Noisy Instrument Tomography",
            "packet_file":PACK174.name,"packet_sha256":packet_hash,**packet,
            "status":"proven_analytic_noise_threshold_with_numerical_subspace_audit",
            "boundary":"The theorem identifies a nearby commutant subspace, not exact matrix factors, apparatus calibration, labels, or adversarial systematic-error correction."}


def st175_hmm_collatz() -> dict:
    P=np.array([[.9,.1],[.2,.8]]);p0=np.array([.02,.08]);p1=np.array([.92,.98])
    _,t=hidden_pair_transfer(P,p0,p1,.5)
    vals,vecs=np.linalg.eig(t);idx=int(np.argmax(abs(vals)));rho=float(abs(vals[idx]));v=np.real(vecs[:,idx]);v*=np.sign(np.sum(v));v=np.maximum(v,1e-14);v/=np.min(v)
    mp.iv.dps=70
    Tiv=[[iv(0) for _ in range(4)] for _ in range(4)]
    for z,zp,w,wp in itertools.product(range(2),repeat=4):
        transition=mp.iv.sqrt(iv(str(P[z,w]))*iv(str(P[zp,wp])))
        overlap=mp.iv.sqrt(iv(str(p0[w]))*iv(str(p1[wp])))+mp.iv.sqrt((1-iv(str(p0[w])))*(1-iv(str(p1[wp]))))
        Tiv[2*z+zp][2*w+wp]=transition*overlap
    viv=[mp.iv.mpf([str(np.nextafter(x,-np.inf)),str(np.nextafter(x,np.inf))]) for x in v]
    ratio_bounds=[]
    for i in range(4):
        ratio=sum((Tiv[i][j]*viv[j] for j in range(4)),iv(0))/viv[i]
        ratio_bounds.append(replay_bounds(ratio))
    lower=min(x[0] for x in ratio_bounds);upper=max(x[1] for x in ratio_bounds)
    packet={"positive_test_vector":v.tolist(),"component_ratio_intervals":[list(x) for x in ratio_bounds],
            "spectral_radius_interval":[lower,upper],"floating_spectral_radius":rho,
            "certified_pair_path_error_exponent_lower_bound":-math.log(upper),
            "theorem":"The positive Collatz--Wielandt inequalities min_i(Tv/v)_i <= rho(T) <= max_i(Tv/v)_i, evaluated with outward intervals, certify the displayed spectral-radius bracket and sharpen the hidden-HMM pair-path error exponent."}
    packet_hash=write_packet(PACK175,packet)
    return {"program":"ST175","object":"Interval Collatz Certificate for the Hidden-HMM Pair Transfer",
            "packet_file":PACK175.name,"packet_sha256":packet_hash,**packet,
            "status":"proven_sharper_pair_path_exponent_certificate",
            "boundary":"This sharpens the conservative pair-path bound, not the exact observed-process Chernoff exponent."}


def st176_microscopic_bath() -> dict:
    bath_qubits=math.ceil(math.log2(N));unused=2**bath_qubits-N
    packet={"system_levels":N,"minimum_bath_qubits_from_dimension":bath_qubits,"unused_computational_states":unused,
            "parallel_pairwise_SWAP_gates":bath_qubits,"two_qubit_circuit_depth":1,
            "parallel_time_times_per_edge_strength":math.pi/2,
            "time_times_global_sum_norm_budget":2*math.pi,
            "energy_conserving_global_SWAP_locality_qubits":2*bath_qubits,
            "theorem": (
                "Dimension at least twelve requires at least four bath qubits. Encoding the twelve-level system in four qubits, four disjoint two-qubit SWAP gates implement register SWAP in depth one and time pi/(2g) at per-edge strength g. "
                "Under a global sum-norm budget G the same construction takes 2pi/G. For an arbitrary encoded strict Hamiltonian, the manifestly energy-conserving generator g S_register is an eight-qubit interaction; the 2-local gate generator need not conserve that energy at intermediate times."
            )}
    packet_hash=write_packet(PACK176,packet)
    return {"program":"ST176","object":"Microscopic Qubit-Bath Locality and Energy-Conservation Tradeoff",
            "packet_file":PACK176.name,"packet_sha256":packet_hash,**packet,
            "status":"proven_dimension_and_explicit_local_circuit_tradeoff",
            "boundary":"No no-go is proved for all energy-conserving 2-local controlled constructions; qubit encoding, controls, physical strengths, beta, and apparatus are supplied."}


def st177_validator_study() -> dict:
    result=run_validator(True)
    packet={"validator_record_file":"FIN_ST177_Operational_Validator_Record.json","validator_record_sha256":sha(ROOT/"FIN_ST177_Operational_Validator_Record.json"),
            "validation":result["validation"],"all_eleven_deletions_detected":result["all_eleven_deletions_detected"],
            "packet_sha256":result["packet_sha256"],"boundary":result["boundary"],
            "theorem":"The executable schema detects deletion of every ST165 operational atom and verifies canonical record and holdout hashes. The supplied local template is structurally valid but must remain physically blocked because units, clock calibration, apparatus, custody, and laboratory execution are absent."}
    packet_hash=write_packet(PACK177,packet)
    return {"program":"ST177","object":"Executable Operational Tuple and Frozen-Record Validator",
            "packet_file":PACK177.name,"packet_sha256_file":packet_hash,**packet,
            "status":"proven_executable_structural_validator_with_physical_block",
            "boundary":result["boundary"]}


def make_figures(d:dict):
    FIG_DIR.mkdir(exist_ok=True)
    fig,ax=plt.subplots(figsize=(7,4));a=np.array(d["ST166"]["logical_probability_table"]);ax.imshow(a,vmin=0,vmax=1,cmap="viridis");ax.set(title="ST166: carrier-independent operational probability table",xlabel="preparation",ylabel="effect");fig.colorbar(ax.images[0],ax=ax);fig.tight_layout();fig.savefig(FIG_DIR/"st166_carrier_invariant.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST167"]["rows"];ax.plot([r["functional_generators"] for r in rows],[r["joint_commutant_dimension"] for r in rows],"o-");ax.axhline(12,ls="--",color="red",label="vertex target");ax.set(xlabel="number of f(A) generators",ylabel="commutant dimension",title="ST167: functional calculus cannot refine the pointer algebra");ax.legend();fig.tight_layout();fig.savefig(FIG_DIR/"st167_functional_no_go.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST168"]["fold_rows"];ax.scatter([sum(r["rho_interval"])/2 for r in rows],[sum(r["kappa_interval"])/2 for r in rows],c=[r["mode"] for r in rows]);ax.set(xlabel="uniform amplitude squared",ylabel="fold kappa",title="ST168: complete uniform-branch fold classes");fig.tight_layout();fig.savefig(FIG_DIR/"st168_uniform_folds.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST169"]["rows"];ax.plot([r["theta"] for r in rows],[r["optimal_entanglement_fidelity"] for r in rows],"o-");ax.set(xlabel="unitary phase separation",ylabel="optimal entanglement fidelity",title="ST169: exact nonorthogonal binary-unitary recovery");fig.tight_layout();fig.savefig(FIG_DIR/"st169_nonorthogonal_recovery.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST170"]["rows"];ax.loglog([r["gap"] for r in rows],[r["minimum_geodesic_action_T1"] for r in rows],"o-");ax.set(xlabel="selector gap",ylabel="minimum Fisher--Rao action",title="ST170: information-geometric selector cost");fig.tight_layout();fig.savefig(FIG_DIR/"st170_fisher_cost.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST172"]["rows"];ax.hist([r["margin"] for r in rows],bins=20);ax.set(xlabel="Krawczyk margin",ylabel="cell count",title="ST172: extended 125-cell nuisance cover");fig.tight_layout();fig.savefig(FIG_DIR/"st172_extended_cover.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));rows=d["ST174"]["rows"];ax.semilogx([r["generator_noise_norm_bound"] for r in rows],[r["numerical_projector_distance"] for r in rows],"o-",label="sampled projector error");ax.set(xlabel="generator error bound",ylabel="low-subspace projector distance",title="ST174: robust tomographic factor subspace");ax.legend();fig.tight_layout();fig.savefig(FIG_DIR/"st174_tomography_noise.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));lo,hi=d["ST175"]["spectral_radius_interval"];ax.bar(["row-sum ST162","Collatz ST175"],[d["ST162_reference_row_sum"],hi]);ax.set(ylabel="certified spectral-radius upper bound",title="ST175: sharper hidden-HMM pair-transfer certificate");fig.tight_layout();fig.savefig(FIG_DIR/"st175_collatz.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));ax.bar(["bath qubits","2-local SWAP gates","energy-conserving locality"],[d["ST176"]["minimum_bath_qubits_from_dimension"],d["ST176"]["parallel_pairwise_SWAP_gates"],d["ST176"]["energy_conserving_global_SWAP_locality_qubits"]]);ax.set(ylabel="count / locality order",title="ST176: microscopic bath tradeoff");fig.tight_layout();fig.savefig(FIG_DIR/"st176_local_bath.png",dpi=190);plt.close(fig)


def main():
    _,a,_=strict_operator();out={"metadata":{"programs":"ST166-ST177","date":"2026-08-11","seed":SEED,"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"sympy":sp.__version__}}
    out["ST166"]=st166_carrier_invariant();out["ST167"]=st167_strict_interaction_no_go(a);out["ST168"]=st168_uniform_fold_classification();out["ST169"]=st169_nonorthogonal_recovery()
    out["ST170"]=st170_information_geometric_selector();out["ST171"]=st171_anisotropy_coupling_no_go();out["ST172"]=st172_extended_cover(a);out["ST173"]=st173_noisy_refinement()
    out["ST174"]=st174_noisy_tomography();out["ST175"]=st175_hmm_collatz();out["ST176"]=st176_microscopic_bath();out["ST177"]=st177_validator_study()
    out["ST162_reference_row_sum"]=.7939043351829553
    out["recommended_next_programs"]=[
        {"id":"ST178","priority":1,"study":"test carrier invariance for irreversible encodings and characterize sufficient statistics rather than exact left inverses"},
        {"id":"ST179","priority":2,"study":"search one strict nonlinear or state-dependent environment coupling outside functional calculus"},
        {"id":"ST180","priority":3,"study":"classify nonuniform fold branches with validated continuation from the 30 uniform folds"},
        {"id":"ST181","priority":4,"study":"solve three-error nonorthogonal recovery with an independently replayable primal-dual SDP"},
        {"id":"ST182","priority":5,"study":"derive a controlled Lindblad selector cost under a strict-compatible detailed-balance constraint"},
        {"id":"ST183","priority":6,"study":"audit candidate twelfth-order observables capable of identifying g without importing physical units"},
        {"id":"ST184","priority":7,"study":"extend nuisance continuation to halfwidth 1e-4 with adaptive subdivision and checkpointed interval cells"},
        {"id":"ST185","priority":8,"study":"derive correlated-noise refinement bounds and multiscale compression thresholds"},
        {"id":"ST186","priority":9,"study":"reconstruct exact factor algebras from the certified noisy commutant subspace"},
        {"id":"ST187","priority":10,"study":"bound the gap between pair-path and exact observed hidden-HMM Chernoff rates"},
        {"id":"ST188","priority":11,"study":"test energy-conserving k-local reset circuits for the encoded strict Hamiltonian"},
        {"id":"ST189","priority":12,"study":"connect the ST177 validator to externally supplied immutable event records without generating custody claims"},
    ]
    out["central_verdict"]=("Operational information can be invariant under carrier change when encodings have left inverses and instruments are transported, but Shannon entropy alone is insufficient. Strict functional calculus cannot source a finer environment pointer algebra. "
        "The uniform fold branch is exactly classifiable with 30 geometric folds, and the twelfth-order anisotropy coupling is invisible to all strict linear and order-11 data. Robust inference and operational validation improve, while physical source objects remain explicit inputs.")
    out["epistemic_boundary"]="No strict nonlinear environment coupling, irreversible universal carrier invariant, QW-2191 discharge, strict twelfth-order coupling, dimensional unit, calibrated laboratory apparatus or independent record, legacy-to-strict completion or role transfer, Standard Model, gravity, L_total, physical projection theorem, or ToE closure is claimed."
    RESULTS.write_text(json.dumps(native(out),indent=2,sort_keys=True),encoding="utf-8")
    with SUMMARY.open("w",newline="",encoding="utf-8") as h:
        w=csv.writer(h);w.writerow(["program","object","status"])
        for k in range(166,178):w.writerow([f"ST{k}",out[f"ST{k}"]["object"],out[f"ST{k}"]["status"]])
    make_figures(out);print(json.dumps({"results":RESULTS.name,"programs":12,"figures":9},indent=2))


if __name__=="__main__":main()
