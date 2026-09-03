#!/usr/bin/env python3
"""FIN ST1017--ST1031: invariant manifolds and annihilation boundaries."""
import math

import numpy as np

from fin_total_nadsoliton_common import N, normalized, write_packet, write_round


LO, HI = 1017, 1031
NAMES = [
    "AnnihilationBoundary_Definition", "UnitSphere_PositiveDistanceBarrier",
    "Simplex_MassBarrier", "ConservedCharge_ZeroNoGo",
    "LocallyLipschitz_FiniteTimeHitNoGo", "NonLipschitz_FiniteTimeCollapse",
    "Lipschitz_AsymptoticExtinction", "CompleteFlow_StateSpaceInvariance",
    "ProjectivePhase_Quotient", "CompactInvariantSet_Persistence",
    "CoerciveLyapunov_Sublevel", "Topology_ComponentProtection",
    "PremiseRemoval_Countermodels", "CurrentFIN_PremiseAudit",
    "RoundTwo_Recommendation",
]


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k - LO], status, boundary, payload)


def main():
    x = {}
    x["ST1017"] = packet(1017, "constructed_annihilation_boundary",
        "A new mathematical object, not a strict-core export.", {
        "definition": "partial_0 M={limit points representing loss of normalization/carrier}",
        "normalized_models": "partial_0 M is absent or at positive distance",
    })
    psi = normalized(np.arange(1.0, N+1.0))
    x["ST1018"] = packet(1018, "proven_unit_sphere_distance_to_zero_is_one",
        "Finite or infinite Hilbert dimension; normalized-state premise required.", {
        "distance": float(np.linalg.norm(psi)), "formula": "inf_{||psi||=1}||psi-0||=1",
    })
    p = np.arange(1.0, N+1.0); p /= p.sum()
    x["ST1019"] = packet(1019, "proven_simplex_l1_distance_to_zero_is_one",
        "Probability normalization premise required.", {
        "distance_l1": float(np.linalg.norm(p,1)), "mass": float(p.sum()),
    })
    x["ST1020"] = packet(1020, "proven_positive_conserved_charge_forbids_zero",
        "Charge must be strict-sourced and conserved; current FIN only has channel-dependent norms/mass.", {
        "theorem": "Q(Phi_t x)=Q(x)=q0>0 and Q(0)=0 implies Phi_t x !=0",
        "logical_strength": "sufficient, not necessary",
    })
    x["ST1021"] = packet(1021, "proven_finite_time_equilibrium_hit_no_go_under_uniqueness",
        "Autonomous locally Lipschitz ODE and zero equilibrium; backward uniqueness on the relevant interval.", {
        "theorem": "if x(T)=0 and f(0)=0, uniqueness implies x(t)=0 on the shared interval",
        "nonzero_initial_can_hit_zero": False,
    })
    x0 = 1.0; T = 2.0*math.sqrt(x0)
    sample_t = 1.0
    nonlip = max(math.sqrt(x0)-sample_t/2.0,0.0)**2
    x["ST1022"] = packet(1022, "proven_nonlipschitz_finite_time_collapse_counterexample",
        "Scalar countermodel, not a FIN law.", {
        "law": "x_dot=-sqrt(x), x>=0", "x0": x0, "hitting_time": T,
        "x_at_t1": nonlip, "lesson": "state-space closure or regularity is indispensable",
    })
    x["ST1023"] = packet(1023, "proven_lipschitz_dynamics_can_extinguish_asymptotically",
        "Shows finite-time no-hit is weaker than eternal positive lower bound.", {
        "law": "x_dot=-x", "solution": "x(t)=x0 exp(-t)", "finite_hit": False,
        "limit_zero": True,
    })
    x["ST1024"] = packet(1024, "proven_complete_invariant_flow_has_no_exit_event",
        "This conclusion is conditional on invariance and completeness.", {
        "premise": "Phi:R x M -> M", "conclusion": "no trajectory leaves M at finite time",
        "warning": "completeness is not implied by self-reference",
    })
    phase = np.exp(1j*0.73)
    overlap = abs(np.vdot(psi, phase*psi))
    x["ST1025"] = packet(1025, "proven_global_phase_is_not_state_annihilation",
        "Applies to projective/ray semantics.", {
        "ray_overlap": float(overlap), "same_ray": bool(abs(overlap-1.0)<1e-12),
        "zero_ray_exists": False,
    })
    x["ST1026"] = packet(1026, "proven_compact_invariant_set_precludes_escape",
        "Compactness does not by itself prove the set is invariant or select its dynamics.", {
        "theorem": "continuous complete orbit in compact invariant K remains in K",
        "does_not_prove": ["localization", "nontrivial recurrence period", "physical realization"],
    })
    x["ST1027"] = packet(1027, "constructed_coercive_barrier_criterion",
        "Candidate theorem schema; FIN source for V remains open.", {
        "criterion": "V(x)->infinity near partial_0 M and V(Phi_t x)<=V(x0)",
        "conclusion": "orbit cannot approach partial_0 M",
    })
    x["ST1028"] = packet(1028, "proven_component_label_protects_only_if_zero_is_outside_component",
        "Topology requires a specified state space and continuous dynamics.", {
        "theorem": "continuous paths cannot cross connected components",
        "failure_mode": "gap closure/singularity changes the admissible space",
    })
    x["ST1029"] = packet(1029, "proven_each_minimal_premise_is_removable",
        "Countermodels establish necessity for the proposed total-state theorem, not metaphysical necessity.", {
        "removals": {
            "exclude_zero": "include absorbing zero and allow decay",
            "invariance": "x_dot=-x on M=(0,infinity) approaches boundary",
            "completeness": "x_dot=-sqrt(x) reaches boundary",
            "nonempty_carrier": "empty model contains no being to persist",
        },
    })
    x["ST1030"] = packet(1030, "current_FIN_supplies_channel_invariants_not_total_ontology",
        "Strict A is real symmetric Laplacian-like; total-state manifold is not yet an exported axiom/theorem.", {
        "supplied": ["unitary norm conservation", "heat mass conservation", "conditional wave energy"],
        "missing": ["canonical total-state type", "proof this type is the nadsoliton", "global completeness for full nonlinear law", "no-killing theorem"],
    })
    x["ST1031"] = packet(1031, "recommended_ST1032_ST1046",
        "Next round tests whether bootstrap/self-reference can derive rather than assume persistence.", {
        "next": ["fixed-point existence versus stability", "internal environment and closed resource circulation",
                 "global versus subsystem entropy", "recurrence and irreversible effective dynamics"],
    })
    write_round(LO, HI, x)


if __name__ == "__main__":
    main()
