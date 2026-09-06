"""FIN ST8549--ST8578: thirty bounded, adaptive analytical investigations.

Numbers here are floating verification, not interval certificates. Proofs and
assumptions are in REPORT_EN.md. No function manufactures a proof from a hash.
"""
from __future__ import annotations
import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.stats import binom


def radial(values):
    n = 2*len(values)
    return np.array([[0 if i == j else values[min(abs(i-j), n-abs(i-j))-1]
                      for j in range(n)] for i in range(n)], dtype=float)


def kernels():
    d = np.arange(1, 7, dtype=float)
    strict = radial(np.cos(.18575*d+.16250)/(1+d**1.8))
    legacy = radial(4*math.log(2)*np.cos(math.pi*d/4+math.pi/6)/(1+.01*d))
    return strict, legacy


def lap(W):
    return np.diag(W.sum(axis=1))-W


def signed_cover(W):
    pos, neg = np.maximum(W, 0), np.maximum(-W, 0)
    C = np.block([[pos, neg], [neg, pos]])
    I = np.eye(len(W))/math.sqrt(2)
    return C, np.vstack((I, I)), np.vstack((I, -I))


def record(k, question, result, evidence, next_reason, status="Proven (conditional; see proof)"):
    return dict(round=k, program=f"ST{8548+k}", question=question,
                result=result, status=status, evidence=evidence,
                next_reason=next_reason)


def first_batch():
    W, V = kernels()
    A, raw = lap(W), lap(V)
    out = []
    assert np.min(W[W != 0]) > 0
    negative = [(d, float(V[0, d])) for d in range(1, 7) if V[0, d] < 0]
    assert negative
    out.append(record(1, "Can both raw kernels be Markov conductances?",
        "Strict is positive; the declared raw Legacy* profile is signed and fails the rate test.",
        dict(strict_weights=W[0, 1:7].tolist(), legacy_weights=V[0, 1:7].tolist(),
             negative_legacy=negative, strict_row_sum=float(W[0].sum()),
             legacy_raw_min_eigenvalue=float(np.linalg.eigvalsh(raw)[0])),
        "Instead of taking absolute values silently, retain signs in a doubled state space."))
    C, E, O = signed_cover(V)
    B = lap(C)
    signed = np.diag(np.abs(V).sum(axis=1))-V
    errs = [np.linalg.norm(B@E-E@lap(abs(V))), np.linalg.norm(B@O-O@signed)]
    assert max(errs) < 1e-12
    shift = 2*np.maximum(-V, 0).sum(axis=1)
    assert np.linalg.norm(signed-raw-np.diag(shift)) < 1e-12
    out.append(record(2, "Can legacy signs be retained in a positive stochastic cover?",
        "A 24-state positive double cover has unsigned-even and signed-odd Laplacian sectors.",
        dict(intertwining_errors=errs, odd_shift=shift.tolist(),
             odd_min_eigenvalue=float(np.linalg.eigvalsh(signed)[0])),
        "A shifted signed amplitude is not a probability: test whether a positive gauge removes the distinction."))
    triangles = [(i,j,k) for i,j,k in itertools.combinations(range(12), 3)
                 if V[i,j]*V[j,k]*V[k,i] < 0]
    assert triangles
    out.append(record(3, "Can a diagonal sign or positive Doob gauge make raw legacy positive?",
        "A negative triangle product obstructs sign switching to all-positive weights; positive diagonal similarity preserves each sign.",
        dict(frustrated_triangles=len(triangles), witness=triangles[0],
             witness_product=float(np.prod([V[triangles[0][0],triangles[0][1]],
                 V[triangles[0][1],triangles[0][2]],V[triangles[0][2],triangles[0][0]]]))),
        "The sign degree of freedom must remain explicit; determine what a coarse observer loses."))
    C2, _, _ = signed_cover(abs(V))
    collapse = np.hstack((np.eye(12), np.eye(12)))
    coarseerr = np.linalg.norm(collapse@lap(C)-lap(abs(V))@collapse)
    scale = float(np.sum(W*abs(V))/np.sum(V*V))
    mismatch = float(np.linalg.norm(W-scale*abs(V))/np.linalg.norm(W))
    assert coarseerr < 1e-12 and mismatch > .1
    out.append(record(4, "Does the sign-cover automatically complete legacy into strict?",
        "Coarse mass is blind to every sign assignment, but its weights are |legacy|, not strict. The construction is not a completion map.",
        dict(coarse_error=float(coarseerr), different_cover_norm=float(np.linalg.norm(C-C2)),
             best_absolute_legacy_scale=scale, relative_strict_mismatch=mismatch),
        "A mathematical dilation is insufficient: test the claimed dual dynamics at the level of vertex records."))
    ts = [1e-2, 1e-3, 1e-4]
    rows = []
    for t in ts:
        heat = expm(-t*A)
        quantum = abs(expm(-1j*t*A))**2
        rows.append(dict(t=t, heat_over_t=float(heat[0,1]/t),
                         born_over_t2=float(quantum[0,1]/t**2)))
    assert abs(rows[-1]["heat_over_t"]-W[0,1]) < 2e-4
    assert abs(rows[-1]["born_over_t2"]-W[0,1]**2) < 1e-7
    out.append(record(5, "Are heat and unitary vertex records the same process?",
        "Off-diagonal heat is t W_ij+O(t^2); a finite Hamiltonian gives t^2 |H_ij|^2+O(t^3). No positive linear clock identifies them.",
        dict(short_time=rows, W01=float(W[0,1]), W01_squared=float(W[0,1]**2)),
        "Try an explicit repeated-measurement limit instead of Wick continuation."))
    return out


def lindblad_action(W, rho, H=None, dephasing=0.):
    """Rank-one edge jumps sqrt(W_ij)|j><i| and vertex dephasing."""
    rates = W.sum(axis=1)
    value = np.diag(W.T@np.diag(rho))-.5*(rates[:,None]+rates[None,:])*rho
    value = value+dephasing*(np.diag(np.diag(rho))-rho)
    if H is not None:
        value = value-1j*(H@rho-rho@H)
    return value


def measurement_batch():
    W, _ = kernels()
    A = lap(W)
    out = []
    Q2 = -lap(W*W)
    target = expm(.7*Q2)
    errors = []
    for n in [20, 80, 320, 1280]:
        delta = .7/n
        step = abs(expm(-1j*math.sqrt(delta)*A))**2
        errors.append(dict(steps=n, error=float(np.linalg.norm(np.linalg.matrix_power(step,n)-target))))
    assert errors[-1]["error"] < errors[0]["error"]/30
    ratio = W[0,1:7]
    out.append(record(6, "Does rapid projective measurement recover the strict heat operator?",
        "With H scaled as delta^(-1/2), repeated vertex measurements converge to rates |H_ij|^2. For H=A these are W_ij^2, not W_ij.",
        dict(convergence=errors, squared_to_original_edge_ratios=ratio.tolist()),
        "Try to compile the target rates explicitly, then test uniqueness of the compiled Hamiltonian."))
    H = np.sqrt(W)
    H2 = H.copy(); H2[0,1] *= -1; H2[1,0] *= -1
    assert np.array_equal(abs(H)**2, abs(H2)**2)
    spectral_difference = np.linalg.norm(np.linalg.eigvalsh(H)-np.linalg.eigvalsh(H2))
    assert spectral_difference > .05
    out.append(record(7, "Can one compile strict heat from a Hamiltonian, uniquely?",
        "Off-diagonal H_ij=sqrt(W_ij) exp(i theta_ij) compiles W in that limit. Cycle phases and diagonals remain free and change coherent dynamics.",
        dict(rate_error=float(np.linalg.norm(abs(H)**2-W)),
             phase_changed_spectral_difference=float(spectral_difference),
             triangle_fluxes=[0, math.pi]),
        "A compiled measurement limit introduces a Hamiltonian and instrument; test an exact quantum channel extension instead."))
    p = np.arange(1,13,dtype=float)/78
    diagonal = np.diag(p).astype(complex)
    change = lindblad_action(W, diagonal)
    assert np.linalg.norm(change-np.diag(-A@p)) < 1e-14
    v = np.sqrt(p).astype(complex); v[1] *= 1j
    rho = np.outer(v,v.conj())
    assert np.linalg.norm(np.diag(lindblad_action(W,rho))+A@p) < 1e-14
    out.append(record(8, "Is there an exact completely positive heat embedding?",
        "Rank-one jumps sqrt(W_ij)|j><i| give a trace-preserving Lindblad semigroup with exactly autonomous strict heat populations for all density matrices.",
        dict(population_error=float(np.linalg.norm(np.diag(lindblad_action(W,rho))+A@p)),
             trace_derivative=float(np.trace(lindblad_action(W,rho)).real)),
        "Exact embedding exists; test whether vertex records determine its coherence law."))
    base = lindblad_action(W,rho)
    extra = lindblad_action(W,rho,dephasing=3.)
    assert np.linalg.norm(np.diag(extra-base)) < 1e-14
    assert np.linalg.norm(extra-base) > .1
    out.append(record(9, "Do complete heat population records select a quantum completion?",
        "Arbitrary nonnegative vertex dephasing changes coherences while preserving every autonomous population trajectory.",
        dict(population_difference=float(np.linalg.norm(np.diag(extra-base))),
             full_derivative_difference=float(np.linalg.norm(extra-base))),
        "Test the desired original coherent generator H=A in the same extension."))
    full = lindblad_action(W,rho,H=A)
    witness = np.diag(full-base)
    assert np.linalg.norm(witness) > .01
    out.append(record(10, "Can H=A be added without altering exact autonomous heat populations?",
        "Not in the declared rank-one-jump plus vertex-dephasing class: autonomous heat for all density matrices forces H to be vertex-diagonal.",
        dict(coherent_population_witness=witness.real.tolist(),
             witness_norm=float(np.linalg.norm(witness))),
        "Do not identify two semigroups by a density-state Wick rotation; test normalization and linearity explicitly."))
    return out


def population_generator(W, N=2, theta=0.):
    """Independent coordinate jumps mixed with common vertex transpositions."""
    n = len(W)
    states = list(itertools.product(range(n), repeat=N))
    index = {x:i for i,x in enumerate(states)}
    G = np.zeros((len(states),len(states)))
    for row,x in enumerate(states):
        for k in range(N):
            for j in range(n):
                if j != x[k]:
                    y = list(x); y[k] = j
                    G[row,index[tuple(y)]] += (1-theta)*W[x[k],j]
        for i in range(n):
            for j in range(i+1,n):
                y = tuple(j if v == i else i if v == j else v for v in x)
                if y != x:
                    G[row,index[y]] += theta*W[i,j]
        G[row,row] = -G[row].sum()
    return G, states


def microscopic_batch():
    W, _ = kernels()
    out = []
    P = np.diag([1., math.exp(-1.)])
    r0, r1 = np.diag([1.,0.]), np.diag([0.,1.])
    def filt(rho):
        X = P@rho@P
        return X/np.trace(X)
    gap = np.linalg.norm(filt((r0+r1)/2)-(filt(r0)+filt(r1))/2)
    assert gap > .1
    out.append(record(11, "Does normalized imaginary-time filtering define a universal physical channel?",
        "It is nonlinear in the density matrix except for scalar filters. Unnormalized P rho P is a trace-decreasing CP operation, not deterministic heat on states.",
        dict(two_level_affinity_defect=float(gap)),
        "Both stochastic and quantum lifts need extra premises. Test whether population consistency removes the stochastic freedom."))
    G0, states = population_generator(W,theta=0.)
    G1, _ = population_generator(W,theta=1.)
    E = np.array([[int(x[0] == i) for i in range(12)] for x in states], dtype=float)
    marginal_error = np.linalg.norm(G1@E-E@(-lap(W)))
    assert marginal_error < 1e-12
    swap = np.array([states.index((y,x)) for x,y in states])
    assert np.linalg.norm(G1[np.ix_(swap,swap)]-G1) < 1e-12
    out.append(record(12, "Does exchangeability plus exact projective consistency force independent walkers?",
        "No. Common vertex transpositions give a consistent exchangeable family with exactly the same single-particle generator Q and synchronized multi-particle jumps.",
        dict(single_marginal_error=float(marginal_error),
             different_two_particle_generator=float(np.linalg.norm(G1-G0))),
        "The fully synchronous law is reducible. Try irreducibility, detailed balance and the full equilibrium law as possible selectors."))
    theta = .35
    G = (1-theta)*G0+theta*G1
    symmetry_error = np.linalg.norm(G-G.T)
    eig = np.linalg.eigvalsh(-G)
    assert symmetry_error < 1e-12 and eig[1] > 0
    out.append(record(13, "Do irreducibility, reversibility and a product-uniform equilibrium select the kinetics?",
        "No. Every mixture 0<=theta<1 is irreducible and reversible with the identical product-uniform equilibrium and exact singleton heat.",
        dict(theta=theta, symmetry_error=float(symmetry_error), gap=float(eig[1]),
             uniform_stationary_error=float(np.linalg.norm(np.ones(144)@G))),
        "Static Shannon entropy is therefore insufficient. Locate a dynamical observable separating these realizations."))
    initial = states.index((0,0))
    f = np.array([x[0] == 1 for x in states], dtype=float)
    g = np.array([x[1] == 1 for x in states], dtype=float)
    def cross(Q):
        return (Q@(f*g)-f*(Q@g)-g*(Q@f))[initial]
    independent, dependent = cross(G0), cross(G)
    assert abs(dependent-theta*W[0,1]) < 1e-14 and independent == 0
    out.append(record(14, "Which infinitesimal observation distinguishes the hidden common dynamics?",
        "Pair coincidence / cross quadratic variation distinguishes the generators even when every single-particle trajectory and equilibrium law agrees.",
        dict(independent_pair_rate=float(independent), common_pair_rate=float(dependent),
             predicted=float(theta*W[0,1])),
        "Ask exactly which extra dynamical premise makes singleton rates determine the full generator."))
    smallW = np.array([[0.,2.,3.],[2.,0.,5.],[3.,5.,0.]])
    small, xs = population_generator(smallW, N=3, theta=0.)
    changed = np.array([[sum(a != b for a,b in zip(x,y)) for y in xs] for x in xs])
    assert np.max(abs(small[(changed>1)])) == 0
    errs = []
    for coordinate in range(3):
        proj = np.array([[x[coordinate] == j for j in range(3)] for x in xs],dtype=float)
        errs.append(float(np.linalg.norm(small@proj-proj@(-lap(smallW)))))
    out.append(record(15, "What minimal local premise yields a unique population generator?",
        "Autonomous singleton generator Q at every configuration plus no simultaneous coordinate changes forces G_N=sum_k Q^(k). Exchangeability is unnecessary.",
        dict(three_coordinate_replay_errors=errs, configurations=len(xs)),
        "Replace the qualitative no-simultaneous-jump premise by a measurable, nonnegative witness."))
    return out


def coincidence_budget(G, states):
    distance = np.array([[sum(a != b for a,b in zip(x,y))
                          for y in states] for x in states])
    return np.sum(G*distance*(distance-1)/2, axis=1)


def locality_batch():
    W, _ = kernels()
    out = []
    smallW = np.array([[0., .7],[.7,0.]])
    Gi, states = population_generator(smallW,N=3,theta=0.)
    G, _ = population_generator(smallW,N=3,theta=.2)
    B = coincidence_budget(G,states)
    assert np.allclose(B,.42)
    out.append(record(16, "Can absence of common jumps be expressed without cancellation?",
        "The pair budget B(x)=sum_y g_xy binom(K(x,y),2) is nonnegative and is zero exactly when every jump changes at most one coordinate.",
        dict(pair_budget=B.tolist(), independent_budget=coincidence_budget(Gi,states).tolist()),
        "A zero witness gives an exact theorem; derive a stable version for small nonzero coincidence rates."))
    matrix_bound = .5*np.max(np.sum(abs(G-Gi),axis=1))
    t = .8
    tv = .5*np.max(np.sum(abs(expm(t*G)-expm(t*Gi)),axis=1))
    assert matrix_bound <= 2*max(B)+1e-14 and tv <= 2*t*max(B)
    out.append(record(17, "Does approximate pair locality bound the full finite-time model error?",
        "With identical autonomous singleton rates, sup_x B(x)<=epsilon implies sup_x TV(P_t(x,.),P_ind,t(x,.))<=min(1,2 epsilon t).",
        dict(generator_tv_norm=float(matrix_bound), epsilon=float(max(B)),
             time=t, measured_tv=float(tv), upper_bound=min(1.,2*t*float(max(B)))),
        "Challenge the bound and its assumptions using many exchangeable particles with small jumps."))
    Ns = [2,4,16,64]
    theta = .35
    rows = []
    for N in Ns:
        independent_rate = (1-theta)*N*.7
        pair_rate = theta*N*.7/2
        mean_rate = independent_rate+2*pair_rate
        variance_rate = independent_rate+4*pair_rate
        rows.append(dict(N=N, mean_per_particle=mean_rate/N,
                         variance_per_particle=variance_rate/N,
                         variance_of_fraction=variance_rate/N**2))
    assert all(abs(r["variance_per_particle"]/.7-1-theta)<1e-14 for r in rows)
    out.append(record(18, "Do exchangeability and a deterministic large-population limit select shot noise?",
        "No. Uniformly chosen pairs with rate theta W/(N-1), mixed with single jumps, have the same mean heat and O(1/N) empirical variance but an altered noise prefactor.",
        dict(pair_model=rows, variance_factor=1+theta),
        "This family is not projectively consistent across N. Test whether moment data can certify the missing event structure directly."))
    rates = {1:.3,2:.2,3:.1}
    mu = sum(k*v for k,v in rates.items())
    variance = sum(k*k*v for k,v in rates.items())
    fourth = sum(k**4*v for k,v in rates.items())
    assert variance >= mu and fourth >= variance
    a,b = W[0,1]/78, 12*W[0,1]/78
    lam = (a-b)/math.log(a/b)
    out.append(record(19, "Which raw-event inequalities survive arbitrary nonnegative jump rates?",
        "Local integer-count cumulants obey kappa2>=|kappa1| and kappa4>=kappa2; equality in the latter excludes all marks with |k|>=2.",
        dict(compound_mean=mu,compound_variance=variance,compound_fourth=fourth,
             strict_quadratic_variance_over_drift=float(2*lam/abs(a-b))),
        "These are instantaneous raw-count statements. Attempt to falsify their use with ordinary detector binning."))
    intensity, window = 2., .4
    q = -math.expm1(-intensity*window)
    variance_over_mean = 1-q
    assert variance_over_mean < 1
    out.append(record(20, "Can an ordinary detector make valid Poisson events appear to violate the local count bound?",
        "Yes. A one-click-per-bin threshold maps Poisson arrivals to Bernoulli clicks, with variance/mean=exp(-lambda Delta)<1. Finite-bin data cannot be inserted into the instantaneous theorem.",
        dict(intensity=intensity,bin_width=window,click_probability=q,
             click_variance_over_mean=variance_over_mean),
        "Separate removable loss from irreversible thresholding and test whether a robust coincidence protocol exists."))
    return out


def observed_double_rate(w, theta, efficiency, delta):
    """Two-state/two-particle example; counts retain binomial jump marks."""
    c1 = 2*w*efficiency*(1-theta*efficiency)
    c2 = theta*w*efficiency**2
    q = -math.expm1(-delta*(c1+c2))-math.exp(-delta*(c1+c2))*delta*c1
    return c1,c2,q


def detector_batch():
    out = []
    eta = .6
    checks = []
    for K in range(1,9):
        value = sum(math.comb(m,2)*math.comb(K,m)*eta**m*(1-eta)**(K-m)
                    for m in range(K+1))
        expected = eta**2*math.comb(K,2)
        assert abs(value-expected)<1e-13
        checks.append(dict(mark=K, detected_pair_mean=value, prediction=expected))
    out.append(record(21, "Which detector imperfection preserves the zero-coincidence theorem?",
        "Independent resolution-preserving Bernoulli thinning gives B_observed=eta^2 B. For eta>0 it preserves the exact zero set; a quantitative bound needs a lower efficiency bound.",
        dict(thinning=checks),
        "Finite bins mix multiple events. Bound that contamination instead of assuming it vanishes."))
    w,theta,delta = .7,.2,.01
    Lambda = 2*w
    c1,c2,q = observed_double_rate(w,theta,eta,delta)
    null = (Lambda*delta)**2/2
    lower = math.exp(-Lambda*delta)*delta*theta*w*eta**2
    assert q > lower and q > null
    out.append(record(22, "Can finite time resolution still test absence of simultaneous changes?",
        "With a finite escape bound Lambda, a single-change null has Pr(at least two detected changes)<= (Lambda Delta)^2/2. A resolved joint event gives a positive order-Delta lower bound.",
        dict(delta=delta,escape_bound=Lambda,null_upper_bound=null,
             alternative_probability=q,alternative_lower_bound=lower),
        "The rate orders separate analytically; compute a finite-sample test without synthetic-data fitting."))
    M,alpha = 200000,.001
    numerical_minimum = int(binom.isf(alpha,M,null))+1
    # A slightly conservative threshold admits simple exact rational bounds.
    threshold = 40
    rational_size = Fraction(11,16)**20
    rational_miss = Fraction(3**49,2**108)
    size = float(binom.sf(threshold-1,M,null))
    power = float(binom.sf(threshold-1,M,q))
    assert rational_size < Fraction(1,1000) and rational_miss < Fraction(1,10**9)
    assert size <= float(rational_size) and power >= 1-float(rational_miss)
    out.append(record(23, "Does the proposed rate-order test have a finite resource specification?",
        "For independently reset trials a predeclared binomial upper-tail test controls the null and has high calculated power in the supplied two-particle model. No physical trial was performed.",
        dict(trials=M,threshold=threshold,numerical_minimum_threshold=numerical_minimum,
             nominal_alpha=alpha,calculated_size=size,
             calculated_power=power,null_probability_bound=null,alternative_probability=q),
        "Attempt to destroy identifiability with an uncalibrated but nondegenerate detector." ,
        status="Proven conditional test rule; numerical tails supplemented by exact rational bounds"))
    out[-1]['evidence'].update(exact_size_upper_bound=str(rational_size),
        exact_miss_upper_bound=str(rational_miss),
        size_upper_bound_float=float(rational_size),miss_upper_bound_float=float(rational_miss))
    theta2 = .3
    eta_null = (2-theta2)/2
    detected_null = 2*w*eta_null
    collapsed_alt = (2-theta2)*w
    assert abs(detected_null-collapsed_alt)<1e-15
    out.append(record(24, "Can unknown event merging exactly hide common dynamics?",
        "Yes. One-click-per-event reporting of the common-flip mixture and independently thinned single flips can produce the same entire nonzero Poisson record law despite different microscopic generators.",
        dict(theta=theta2,null_efficiency=eta_null,
             common_observed_rate=collapsed_alt,independent_observed_rate=detected_null),
        "Even a mark-preserving detector may leave a clock/efficiency ambiguity. Test the full record likelihood, not just its mean."))
    base = observed_double_rate(w,.2,.6,delta)
    other = observed_double_rate(2*w,.4,.3,delta)
    assert np.linalg.norm(np.array(base)-np.array(other))<1e-14
    out.append(record(25, "Does full resolved-count likelihood determine both clock and efficiency?",
        "No. The two-particle count law has only c1,c2; a simultaneous change of w,theta,eta preserves both and hence every counting distribution.",
        dict(original=dict(w=w,theta=.2,eta=.6,c1=base[0],c2=base[1]),
             indistinguishable=dict(w=2*w,theta=.4,eta=.3,c1=other[0],c2=other[1])),
        "Operational nonuniqueness is explicit. Audit the earlier nonlinear/prism and inverse-action alternatives before the final synthesis."))
    return out


def parity_rates(rule, theta=.3, field=0.):
    states = list(itertools.product((-1,1),repeat=3))
    energy = np.array([-theta*np.prod(x)-field*sum(x) for x in states])
    Q = np.zeros((8,8))
    for i,x in enumerate(states):
        for k in range(3):
            y = list(x); y[k] *= -1
            j = states.index(tuple(y))
            de = energy[j]-energy[i]
            Q[i,j] = min(1.,math.exp(-de)) if rule == "metropolis" else 1/(1+math.exp(de))
        Q[i,i] = -Q[i].sum()
    return Q,energy


def final_batch():
    from scipy.optimize import root
    W,_ = kernels()
    out = []
    hb,E = parity_rates("heatbath",field=.17)
    met,_ = parity_rates("metropolis",field=.17)
    mask = hb>0
    ratios = np.unique(np.round(met[mask]/hb[mask],12))
    pi = np.exp(-E); pi /= pi.sum()
    assert len(ratios)==2 and np.linalg.norm(pi@met)<1e-14
    hb0,_ = parity_rates("heatbath")
    met0,_ = parity_rates("metropolis")
    olderr = np.linalg.norm(met0-(1+math.exp(-.6))*hb0)
    assert olderr < 1e-14
    out.append(record(26, "Can the failed pure-parity kinetic counterexample be repaired honestly?",
        "Adding an explicitly declared symmetric field gives two flip-energy magnitudes, so Metropolis and heat-bath cease to differ only by a clock. The original pure-parity proportionality is retained.",
        dict(field=.17,rate_ratios=ratios.tolist(),stationary_error=float(np.linalg.norm(pi@met)),
             old_proportionality_error=float(olderr)),
        "Equilibrium geometry fails to select kinetics. Test the analogous inverse-propagator-to-action claim."))
    M = lap(W)+.5*np.eye(12)
    Green = np.linalg.inv(M)
    J = np.linspace(-.5,1.,12)
    v = Green@J
    coupling = .2
    rows = []
    for eps in [.1,.05,.025]:
        sol = root(lambda f:M@f+4*coupling*f**3-eps*J,eps*v,
                   jac=lambda f:M+np.diag(12*coupling*f*f),tol=1e-11)
        assert np.linalg.norm(M@sol.x+4*coupling*sol.x**3-eps*J)<1e-10
        expansion = eps*v-4*coupling*eps**3*(Green@(v**3))
        rows.append(dict(source_scale=eps,cubic_expansion_error=float(np.linalg.norm(sol.x-expansion))))
    assert rows[-1]["cubic_expansion_error"] < rows[0]["cubic_expansion_error"]/500
    out.append(record(27, "Does a supplied propagator determine the nonlinear variational theory?",
        "No. Every S_lambda=phi^T M phi/2-J^T phi+lambda sum phi_i^4, lambda>=0, has the identical zero-source Hessian and linear response M^-1 but different cubic response.",
        dict(coupling=coupling,cubic_checks=rows),
        "Test whether hidden-state compression can select the missing dynamical response from the same static Green data."))
    a = 2.
    examples=[]
    for c,d in [(1.,1.),(2.,3.)]:
        M2=np.array([[a+c*c/d,c],[c,d]])
        value = float(np.linalg.inv(M2)[0,0])
        slope = 1+c*c/(d*d)
        assert abs(value-1/a)<1e-14 and np.linalg.eigvalsh(M2)[0]>0
        examples.append(dict(c=c,d=d,static_Green=value,effective_slope=slope,
                             memory_amplitude=c*c,memory_decay=d))
    out.append(record(28, "Does exact static Green/Schur compression select a hidden memory law?",
        "No. Positive two-variable actions have identical coarse static Green function but different rational response and exponential memory; iteration does not remove the unsourced hidden parameters.",
        dict(positive_completions=examples),
        "The surviving ambiguities occur in signs, quantum instruments, collective jumps, interactions and memory. Formulate only the no-go actually supported by explicit models."))
    G0,xs=population_generator(W,theta=0.)
    G,_=population_generator(W,theta=.35)
    initial=xs.index((0,0)); destination=xs.index((1,1))
    t=.02
    alt=float(expm(t*G)[initial,destination])
    null=float(expm(t*G0)[initial,destination])
    assert alt > 3*null
    out.append(record(29, "Can the audited finite FIN data logically determine a unique collective physical realization?",
        "Not in the stated finite Markov completion class: same W, all singleton paths, product equilibrium, detailed balance, projective consistency and FIN graph symmetries coexist with different pair records.",
        dict(time=t,independent_pair_probability=null,common_pair_probability=alt),
        "This is not a universal impossibility of future physics. Close with a premise-removal and robustness audit of the smallest demonstrated positive bridge."))
    rng=np.random.default_rng(8549)
    stress=[]
    for _ in range(12):
        R=rng.uniform(.01,1.,(3,3)); R=(R+R.T)/2; np.fill_diagonal(R,0)
        theta=float(rng.uniform(0,1))
        G0,xs=population_generator(R,N=2,theta=0.)
        G,xs=population_generator(R,N=2,theta=theta)
        E=np.array([[x[0]==j for j in range(3)] for x in xs],float)
        err=np.linalg.norm(G@E-E@(-lap(R)))
        B=coincidence_budget(G,xs)
        tv=.5*np.max(np.sum(abs(expm(.3*G)-expm(.3*G0)),axis=1))
        bound=min(1.,.6*float(max(B)))
        assert err<1e-12 and tv<=bound+1e-13
        stress.append(dict(marginal_error=float(err),tv=float(tv),bound=bound))
    out.append(record(30, "Which rigorously narrowed bridge survives premise removal and perturbations?",
        "In a declared finite CTMC, autonomous singleton Q plus zero pair-coincidence budget uniquely determines the generator; a small budget gives a stable TV bound. Initial law, clock, detector and the source of zero budget remain separate obligations.",
        dict(random_graph_stress=stress,removed_no_joint_witness="rounds 12--14",
             removed_singleton_rate_witness="2 times the independent generator",
             remaining_source="No FIN theorem presently forces the pair budget to vanish."),
        "Next: seek a strict-derived composition law forcing the nonnegative pair budget to vanish, or prove its independence in a precisely declared enlarged source class."))
    return out


def run(through=30):
    out = first_batch()
    if through >= 10:
        out += measurement_batch()
    if through >= 15:
        out += microscopic_batch()
    if through >= 20:
        out += locality_batch()
    if through >= 25:
        out += detector_batch()
    if through >= 30:
        out += final_batch()
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--through", type=int, choices=(5,10,15,20,25,30),default=30)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("results.json"))
    args = parser.parse_args()
    results = run(args.through)
    args.output.write_text(json.dumps(results, indent=2, allow_nan=False)+"\n")
    for r in results:
        print(f"{r['round']:02d} {r['program']}: {r['result']}")
