"""ST8548: Poisson current action versus quadratic log-mean representation.

All microscopic conclusions assume independent continuous-time Markov walkers
with supplied rates W. Neither that interpretation nor a physical clock is
derived from the finite FIN kernel.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.special import gammaln


def strict_weights():
    W = np.zeros((12, 12))
    for i in range(12):
        for j in range(i+1, 12):
            d = min(j-i, 12-j+i)
            W[i, j] = W[j, i] = math.cos(.18575*d+.1625)/(1+d**(9/5))
    return W


def generator(W):
    return W-np.diag(W.sum(axis=1))  # row-generator convention


def current_covariance(W, p, quadratic=False):
    """Covariance rate of empirical occupations per particle/exposure."""
    C = np.zeros_like(W)
    for i in range(len(p)):
        for j in range(i+1,len(p)):
            value = 2*log_mean(p[i],p[j]) if quadratic else p[i]+p[j]
            C[i,j] = C[j,i] = W[i,j]*value
    return -generator(C)


def log_mean(a, b):
    if min(a, b) <= 0:
        raise ValueError("Rates must be strictly positive")
    if a == b:
        return a
    hi, lo = max(a,b), min(a,b)
    z = (hi-lo)/lo
    return lo*z/math.log1p(z) if z < .1 else (hi-lo)/(math.log(hi)-math.log(lo))


def hamiltonian(s, a, b):
    return a*math.expm1(s)+b*math.expm1(-s)


def optimum_tilt(j, a, b):
    return math.asinh(j/(2*math.sqrt(a*b)))-.5*math.log(a/b)


def current_cost(j, a, b):
    if min(a,b) <= 0:
        raise ValueError("Interior rate formula only")
    g = 2*math.sqrt(a*b)
    return j*math.asinh(j/g)-.5*j*math.log(a/b)-math.hypot(j,g)+a+b


def poisson_entropy(q, rate):
    if q < 0 or rate <= 0:
        raise ValueError("Invalid flux or rate")
    return rate if q == 0 else q*math.log(q/rate)-q+rate


def optimal_traffic(j, a, b):
    root = math.hypot(j, 2*math.sqrt(a*b))
    return .5*(root+j), .5*(root-j)


def cosh_dual(force, a, b):
    # subtracting cosh(0) in a cancellation-safe form
    return 4*math.sqrt(a*b)*math.sinh(force/2)**2


def symmetric_cost(j, a, b):
    g = 2*math.sqrt(a*b)
    return j*math.asinh(j/g)-math.hypot(j,g)+g


def quadratic_cost(j, a, b):
    # Chosen normalization matches both the mean and the time-reversal asymmetry.
    return (j-(a-b))**2/(4*log_mean(a,b))


def contraction_fluxes(j, forward, reverse):
    s = optimum_tilt(j, sum(forward), sum(reverse))
    return np.asarray(forward)*math.exp(s)-np.asarray(reverse)*math.exp(-s)


def net_count_moments(W, p, edge, time, order=4):
    """Finite-time moments from derivatives of a tilted full Markov generator.

    This does not replace finite-time closed-chain counts by independent
    Poisson variables. The backward tilted semigroup is differentiated by
    a block-triangular linear system, evaluated by a matrix exponential.
    """
    n = len(p)
    i,j = edge
    derivatives = [generator(W)]
    for k in range(1, order+1):
        D = np.zeros_like(W)
        D[i,j] = W[i,j]
        D[j,i] = (-1)**k*W[j,i]
        derivatives.append(D)
    T = np.zeros(((order+1)*n,(order+1)*n))
    for row in range(order+1):
        for col in range(row+1):
            T[row*n:(row+1)*n,col*n:(col+1)*n] = (
                math.comb(row,col)*derivatives[row-col])
    initial = np.zeros((order+1)*n)
    initial[:n] = 1
    blocks = (expm(time*T)@initial).reshape(order+1,n)
    raw = blocks@p
    mu1, mu2, mu3, mu4 = raw[1:5]
    variance = mu2-mu1**2
    fourth = mu4-4*mu3*mu1-3*mu2**2+12*mu2*mu1**2-6*mu1**4
    return {"mean": float(mu1), "variance": float(variance),
            "fourth_cumulant": float(fourth)}


def configuration_rule(rule, theta=.3):
    """Recompute the old three-bit parity example to audit the scale claim."""
    import itertools
    states = list(itertools.product((-1,1), repeat=3))
    energy = np.array([-theta*np.prod(s) for s in states])
    Q = np.zeros((8,8))
    for i,s in enumerate(states):
        for k in range(3):
            x = list(s); x[k] *= -1
            j = states.index(tuple(x))
            d = energy[j]-energy[i]
            if rule == "heatbath":
                Q[i,j] = 1/(1+math.exp(d))
            elif rule == "metropolis":
                Q[i,j] = min(1,math.exp(-d))
            elif rule == "symmetric":
                Q[i,j] = math.exp(-d/2)
            else:
                raise ValueError("Unknown rule")
        Q[i,i] = -Q[i].sum()
    return Q


def run():
    W = strict_weights()
    p = np.arange(1,13,dtype=float)/78
    edge = (0,11)
    a,b = W[edge]*p[edge[0]],W[edge]*p[edge[1]]
    j0 = a-b
    m = log_mean(a,b)
    z = .5*math.log(a/b)
    fin_times = []
    for t in [.1,.01,.001,.0001,.00001]:
        moment = net_count_moments(W,p,edge,t)
        fin_times.append({"time":t,**moment,
                         "variance_rate":moment["variance"]/t,
                         "fourth_cumulant_rate":moment["fourth_cumulant"]/t})
    # Correlated fine state: different child laws at edge endpoints.
    forward = a*np.array([.8,.2])
    reverse = b*np.array([.25,.75])
    target = .04
    split = contraction_fluxes(target,forward,reverse)
    contracted = sum(current_cost(j,aa,bb) for j,aa,bb in zip(split,forward,reverse))
    coarse = current_cost(target,a,b)
    fine_m = sum(log_mean(aa,bb) for aa,bb in zip(forward,reverse))
    quad_contract = (target-j0)**2/(4*fine_m)
    quad_coarse = quadratic_cost(target,a,b)
    noise_defect = current_covariance(W,p)-current_covariance(W,p,quadratic=True)
    P = p[:,None]*np.array([[.5,.5]])
    v = np.cos(2*np.pi*np.arange(12)/12)
    P *= 1+.6*v[:,None]*np.array([[1,-1]])
    Wf = np.kron(W,np.eye(2))+np.kron(np.eye(12),np.array([[0.,.2],[.2,0.]]))
    C = np.kron(np.eye(12),np.ones((1,2)))
    micro_projection_error = float(np.max(abs(C@current_covariance(Wf,P.ravel())@C.T-current_covariance(W,P.sum(axis=1)))))
    quadratic_projection_defect = current_covariance(W,P.sum(axis=1),True)-C@current_covariance(Wf,P.ravel(),True)@C.T
    Qh = configuration_rule("heatbath")
    old_counter = {
        "theta":.3, "metropolis_factor":1+math.exp(-.6),
        "symmetric_factor":2*math.cosh(.3),
        "metropolis_residual":float(np.max(abs(configuration_rule("metropolis")-(1+math.exp(-.6))*Qh))),
        "symmetric_residual":float(np.max(abs(configuration_rule("symmetric")-2*math.cosh(.3)*Qh)))}
    # Stationary independent-walker count statistics: finite Stirling check.
    pi = np.full(12,1/12)
    D = float(np.dot(p,np.log(p/pi)))
    sanov = []
    for N in [78,780,7800,78000]:
        counts = np.rint(N*p).astype(int)
        log_prob = gammaln(N+1)-np.sum(gammaln(counts+1))+np.dot(counts,np.log(pi))
        sanov.append({"particles":N,"minus_log_probability_per_particle":float(-log_prob/N),
                      "relative_entropy":D})
    return {
        "programme":"ST8548",
        "scope":"Conditional independent Markov walkers; exact local action and finite-time tilted-generator checks.",
        "strict_edge":{
            "vertices":list(edge),"W":float(W[edge]),"p_i":float(p[edge[0]]),"p_j":float(p[edge[1]]),
            "a":float(a),"b":float(b),"drift":float(j0),
            "poisson_local_variance":float(a+b),"logmean_quadratic_variance":float(2*m),
            "variance_ratio":float((a+b)/(2*m)),
            "shot_variance_to_absolute_drift":float((a+b)/abs(j0)),
            "quadratic_variance_to_absolute_drift":float(2*m/abs(j0)),
            "quadratic_violates_integer_jump_bound":bool(2*m<abs(j0)),
            "z_coth_z":float(z/math.tanh(z)),
            "poisson_fourth_cumulant_rate":float(a+b),"quadratic_fourth_cumulant_rate":0,
            "local_cost_at_drift":current_cost(j0,a,b)},
        "finite_time_counts":fin_times,
        "correlated_refinement":{
            "forward_rates":forward.tolist(),"backward_rates":reverse.tolist(),"target_current":target,
            "optimal_split":split.tolist(),"current_constraint_error":float(abs(sum(split)-target)),
            "poisson_contraction_error":float(abs(contracted-coarse)),
            "poisson_coarse_cost":float(coarse),
            "quadratic_contracted_cost":float(quad_contract),"quadratic_coarse_cost":float(quad_coarse),
            "coarse_minus_fine_logmean":float(m-fine_m)},
        "graph_covariance":{
            "shot_minus_onsager_eigenvalues":np.linalg.eigvalsh(noise_defect).tolist(),
            "exact_covariance_projection_error":micro_projection_error,
            "quadratic_projection_defect_norm":float(np.linalg.norm(quadratic_projection_defect,2))},
        "old_three_bit_rate_claim_correction":old_counter,
        "stationary_multinomial_check":sanov,
        "new_physics_claim":False
    }


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=Path("results.json"))
    args=parser.parse_args()
    results=run()
    args.output.write_text(json.dumps(results,indent=2,sort_keys=True)+"\n")
    print(json.dumps(results,indent=2))
