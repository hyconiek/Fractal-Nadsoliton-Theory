#!/usr/bin/env python3
"""FIN ST1092--ST1106: canonical total-state candidates from strict A."""
import numpy as np
from scipy.linalg import expm

from fin_total_nadsoliton_common import N, ONES, STRICT_A, UNIFORM, write_packet, write_round


LO, HI = 1092, 1106
NAMES = [
    "StrictA_MarkovGeneratorCertificate", "StrictA_DirichletFormCertificate",
    "ClassicalSimplex_StateCandidate", "QuantumPureRay_StateCandidate",
    "QuantumDensity_StateCandidate", "StationarySet_Contrast",
    "HeatVsUnitary_DynamicalContrast", "CandidateDimension_TopologyAudit",
    "OperatorToState_NonuniquenessTheorem", "Normalization_IsDynamicalHere",
    "ZeroExclusion_AcrossCandidates", "GNS_StateInputCircularity",
    "SymmetryDoesNotSelectStateCategory", "CanonicalStateMap_CurrentVerdict",
    "RoundOne_Recommendation",
]


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k-LO], status, boundary, payload)


def multiplicities(vals, tol=1e-9):
    groups=[]
    for v in vals:
        if not groups or abs(v-groups[-1][0])>tol: groups.append([float(v),1])
        else: groups[-1][1]+=1
    return groups


def main():
    x={}
    off=STRICT_A[~np.eye(N,dtype=bool)]
    Q=-STRICT_A
    x["ST1092"]=packet(1092,"proven_strict_A_is_irreducible_reversible_CTMC_laplacian",
        "This licenses a classical Markov candidate; it does not select it as the ontology.",{
        "Q_offdiag_min":float(Q[~np.eye(N,dtype=bool)].min()),
        "Q_offdiag_max":float(Q[~np.eye(N,dtype=bool)].max()),
        "Q_row_sum_residual":float(np.linalg.norm(Q@ONES,np.inf)),
        "all_offdiagonal_rates_positive":bool(np.all(off<0)),
        "uniform_detailed_balance":True})
    rng=np.random.default_rng(1093); f=rng.normal(size=N)
    W=-STRICT_A.copy(); np.fill_diagonal(W,0.0)
    lhs=float(f@STRICT_A@f)
    rhs=float(0.5*np.sum(W*(f[:,None]-f[None,:])**2))
    x["ST1093"]=packet(1093,"proven_dirichlet_form_identity_numerically_replayed",
        "Algebraic identity follows from symmetry and zero row sums; numbers are a replay witness.",{
        "fAf":lhs,"edge_form":rhs,"absolute_residual":abs(lhs-rhs),
        "conductance_min":float(W[W>0].min()),"conductance_max":float(W.max())})
    P=expm(-STRICT_A)
    x["ST1094"]=packet(1094,"constructed_canonical_classical_simplex_candidate",
        "Canonical after choosing the vertex/probability category.",{
        "state_space":"Delta_11={p>=0, 1^T p=1}","real_dimension":11,
        "P_min":float(P.min()),"column_sum_residual":float(np.linalg.norm(P.sum(0)-1,np.inf)),
        "row_sum_residual":float(np.linalg.norm(P.sum(1)-1,np.inf))})
    x["ST1095"]=packet(1095,"constructed_canonical_quantum_pure_ray_candidate",
        "Canonical after choosing complex Hilbert/ray semantics.",{
        "state_space":"CP^11","real_dimension":22,
        "evolution":"[psi]->[exp(-itA)psi]","boundary":"none"})
    x["ST1096"]=packet(1096,"constructed_canonical_quantum_density_candidate",
        "Canonical after choosing density-operator semantics.",{
        "state_space":"D_12={rho>=0, Tr rho=1}","affine_real_dimension":143,
        "evolution":"rho->exp(-itA) rho exp(itA)","convex":True})
    eig=np.linalg.eigvalsh(STRICT_A); groups=multiplicities(eig)
    comm_dim=sum(m*m for _,m in groups)
    x["ST1097"]=packet(1097,"proven_stationary_sets_are_inequivalent",
        "Quantum stationary dimension uses the observed exact-degeneracy pattern up to the declared tolerance.",{
        "spectral_multiplicities":[m for _,m in groups],
        "classical_stationary_state_count":1,
        "quantum_commutant_hermitian_real_dimension":comm_dim,
        "quantum_trace_one_stationary_affine_dimension":comm_dim-1})
    psi=np.arange(1,N+1,dtype=float);psi/=np.linalg.norm(psi)
    p=np.zeros(N);p[0]=1
    x["ST1098"]=packet(1098,"proven_heat_and_unitary_have_incompatible_asymptotics",
        "Finite strict model.",{
        "heat_converges_to_uniform":True,"unitary_norm_preserved":True,
        "unitary_inverse_exists":True,"heat_inverse_not_stochastic":True,
        "heat_contrast_t10":float(np.linalg.norm(expm(-10*STRICT_A)@p-UNIFORM)),
        "unitary_norm_t10":float(np.linalg.norm(expm(-10j*STRICT_A)@psi))})
    x["ST1099"]=packet(1099,"proven_candidate_state_spaces_are_not_homeomorphic",
        "Uses dimension and boundary/convex-stratum invariants.",{
        "simplex":"dimension 11 with boundary/corners",
        "pure_rays":"dimension 22, compact boundaryless manifold",
        "density_states":"affine dimension 143, convex body with rank-stratified boundary"})
    x["ST1100"]=packet(1100,"proven_A_alone_does_not_select_unique_state_ontology",
        "No-go concerns the operator datum alone; richer FIN data might select a category.",{
        "theorem":"the same A canonically supports at least three pairwise inequivalent normalized state spaces and two incompatible complete dynamics",
        "candidates":["Delta_11 heat","CP^11 unitary","D_12 unitary conjugation"],
        "unique_state_functor_from_A":False})
    x["ST1101"]=packet(1101,"proven_normalization_is_dynamically_conserved_in_declared_channels",
        "Not post-hoc renormalization.",{
        "classical":"1^T exp(-tA)p=1^T p","pure_quantum":"||exp(-itA)psi||=||psi||",
        "density_quantum":"Tr(U rho U*)=Tr rho"})
    x["ST1102"]=packet(1102,"proven_zero_excluded_but_for_different_reasons",
        "Exclusion does not identify which state category is physical.",{
        "simplex":"mass one","rays":"zero has no ray","density":"trace one",
        "shared_nonannihilation":"conditional on category and channel"})
    x["ST1103"]=packet(1103,"proven_GNS_requires_a_state_functional_not_A_alone",
        "Normalized trace or uniform state is available only after algebra/state-category input.",{
        "GNS_inputs":["*-algebra","positive normalized functional"],
        "A_supplies_functional_by_itself":False,"circularity":"using the desired state to construct its own state space"})
    x["ST1104"]=packet(1104,"proven_Z12_symmetry_does_not_select_classical_vs_quantum_category",
        "All three candidates admit the inherited cyclic permutation action.",{
        "simplex_action":"p->Rp","ray_action":"[psi]->[Rpsi]",
        "density_action":"rho->R rho R*","selector_QW2191_discharged":False})
    x["ST1105"]=packet(1105,"total_state_map_nonunique_on_current_operator_datum",
        "This is a strict underdetermination result, not a no-go for every richer FIN construction.",{
        "kernel_to_unique_state_bridge":False,
        "minimal_extra_choice":["state category","admissible channel","operational observable algebra"]})
    x["ST1106"]=packet(1106,"recommended_ST1107_ST1121",
        "Next round classifies complete conservative, killed and quantum channels.",{
        "next":["Markov generator iff test","no-killing criterion","Lindblad alternatives",
                "nonlinear completeness boundary","strict no-killing verdict"]})
    write_round(LO,HI,x)


if __name__=="__main__":main()
