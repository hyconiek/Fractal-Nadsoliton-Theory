#!/usr/bin/env python3
"""FIN ST1032--ST1046: bootstrap, self-reference, and internal circulation."""
import math

import numpy as np
from scipy.linalg import expm

from fin_total_nadsoliton_common import write_packet, write_round


LO, HI = 1032, 1046
NAMES = [
    "Bootstrap_FixedPoint_NotStability", "IdentityMap_MaximalNonselection",
    "UnstableSelfConsistentFixedPoint", "Contraction_StableBootstrap",
    "SelfSourcing_CollapseCounterexample", "ClosedSkewCirculation",
    "InternalPump_GlobalConservation", "SubsystemDissipation_GlobalPersistence",
    "GlobalPure_SubsystemMixed", "CoarsePinching_InformationLoss",
    "FiniteClosed_AlmostRecurrence", "PositiveFeedback_NotPersistence",
    "BootstrapPersistence_MinimalTheorem", "CurrentFIN_BootstrapGap",
    "RoundThree_Recommendation",
]


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k - LO], status, boundary, payload)


def entropy(eig):
    eig = np.asarray(eig, dtype=float)
    eig = eig[eig > 0]
    return float(-np.sum(eig * np.log(eig)))


def main():
    x = {}
    x["ST1032"] = packet(1032, "proven_fixed_point_existence_does_not_imply_dynamical_stability",
        "General mathematical obstruction; no FIN-specific fixed-point law is assumed.", {
        "counterexample": "F(x)=2x has fixed point 0 but multiplier 2",
        "fixed_point": 0.0, "spectral_radius_derivative": 2.0,
    })
    x["ST1033"] = packet(1033, "proven_self_consistency_can_fail_to_select",
        "Identity-map toy model.", {
        "map": "F(x)=x", "fixed_point_set": "all x", "selector_exported": False,
        "lesson": "bootstrap equality alone supplies neither uniqueness nor persistence",
    })
    x["ST1034"] = packet(1034, "proven_unique_bootstrap_can_be_unstable",
        "Discrete iteration example.", {
        "map": "F(x)=2x", "unique_fixed_point": 0.0,
        "orbit_x0_0.1_after_5": 3.2, "stable": False,
    })
    orbit = [0.0]
    for _ in range(12): orbit.append(0.5*orbit[-1]+0.5)
    x["ST1035"] = packet(1035, "proven_contraction_bootstrap_is_unique_and_attracting",
        "Conditional architecture; contraction is not strict-sourced.", {
        "map": "F(x)=x/2+1/2", "fixed_point": 1.0,
        "iteration_12": orbit[-1], "error_12": 1.0-orbit[-1],
        "theorem": "Banach contraction gives existence, uniqueness, and convergence",
    })
    x["ST1036"] = packet(1036, "proven_self_coupling_can_still_decay_to_zero",
        "Counterexample to verbal self-sourcing arguments.", {
        "law": "x_dot=-x (state determines its own decay rate)",
        "x0": 1.0, "x_t10": math.exp(-10.0), "limit_zero": True,
    })
    J = np.array([[0.0,1.0,-1.0],[-1.0,0.0,1.0],[1.0,-1.0,0.0]])
    z0 = np.array([1.0,0.2,-0.4])
    zt = expm(2.3*J) @ z0
    x["ST1037"] = packet(1037, "proven_closed_skew_circulation_preserves_global_quadratic_resource",
        "Three-sector toy model, not derived from strict A.", {
        "J_skew_residual": float(np.linalg.norm(J+J.T)),
        "J_row_sum_residual": float(np.linalg.norm(J@np.ones(3))),
        "norm2_initial": float(z0@z0), "norm2_after": float(zt@zt),
        "total_initial": float(z0.sum()), "total_after": float(zt.sum()),
    })
    transfer = np.array([-0.3, 0.1, 0.2])
    x["ST1038"] = packet(1038, "constructed_internal_pump_as_zero_sum_transfer",
        "A subsystem pump is admissible only with its source sector and backreaction included.", {
        "sector_rate": transfer.tolist(), "global_rate": float(transfer.sum()),
        "external_creation": False,
    })
    x["ST1039"] = packet(1039, "proven_effective_subsystem_decay_is_compatible_with_total_persistence",
        "Requires an internal partition and coupling; it does not derive a FIN environment.", {
        "model": "amplitude/information leaves subsystem S and enters internal sector E",
        "global_conserved": True, "subsystem_monotone_possible": True,
        "lesson": "local death does not entail death of the whole nadsoliton",
    })
    x["ST1040"] = packet(1040, "proven_global_purity_can_coexist_with_subsystem_entropy",
        "Two-qubit Bell-state witness, used only as a finite information theorem.", {
        "global_spectrum": [1.0,0.0,0.0,0.0], "reduced_spectrum": [0.5,0.5],
        "global_entropy": entropy([1,0,0,0]), "subsystem_entropy": entropy([0.5,0.5]),
    })
    rho = np.array([[0.5,0.5],[0.5,0.5]])
    pinched = np.diag(np.diag(rho))
    x["ST1041"] = packet(1041, "proven_pinching_erases_accessible_coherence_without_erasing_state",
        "Choice of observable algebra/basis is an additional operational structure.", {
        "purity_before": float(np.trace(rho@rho)), "purity_after": float(np.trace(pinched@pinched)),
        "trace_before": float(np.trace(rho)), "trace_after": float(np.trace(pinched)),
    })
    x["ST1042"] = packet(1042, "proven_finite_closed_unitary_has_almost_recurrence_not_immortality_theorem",
        "Recurrence requires a finite/discrete spectral setting and says nothing about ontology outside the model.", {
        "mechanism": "compact torus of finitely many eigenphases",
        "protects": "arbitrarily close returns", "does_not_protect": "a named localized identity at every time",
    })
    x["ST1043"] = packet(1043, "proven_positive_feedback_is_neither_necessary_nor_sufficient",
        "Scalar counterexamples.", {
        "not_sufficient": "x_dot=x blows up and need not define a complete total state",
        "not_necessary": "x_dot=Jx with J skew preserves norm without positive gain",
        "correct_object": "invariant complete dynamics, not the sign word 'positive'",
    })
    x["ST1044"] = packet(1044, "constructed_bootstrap_persistence_theorem_schema",
        "Conditional theorem, not a strict FIN closure.", {
        "premises": ["F maps a complete invariant set K into itself", "F is a contraction or has certified stable index",
                     "the generated evolution is globally well posed", "K excludes the annihilation boundary"],
        "conclusion": "a unique/stable self-consistent state persists inside K",
    })
    x["ST1045"] = packet(1045, "blocked_no_strict_bootstrap_stability_source",
        "Current repository has self-reference language but no single certified full-state contraction/index/completeness theorem.", {
        "strict_fixed_point_map": False, "strict_stability_certificate": False,
        "strict_internal_resource_partition": False,
    })
    x["ST1046"] = packet(1046, "recommended_ST1047_ST1061",
        "Next round asks what total annihilation could mean to an observer who is inside the total state.", {
        "next": ["record and observer typing", "vacuum versus nonexistence", "coarse blindness and pattern death",
                 "dual unitary/heat channels", "internal falsifiability obstruction"],
    })
    write_round(LO, HI, x)


if __name__ == "__main__":
    main()
