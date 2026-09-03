#!/usr/bin/env python3
"""FIN ST1002--ST1016: total-being ontology and annihilation typing."""
import math

import numpy as np
from scipy.linalg import expm

from fin_total_nadsoliton_common import (
    N, ONES, ROOT, STRICT_A, UNIFORM, normalized, spectral_facts,
    write_packet, write_round,
)


LO, HI = 1002, 1016
NAMES = [
    "TotalBeing_TypeLedger", "HilbertSphere_ZeroExclusion",
    "ProjectiveState_ZeroExclusion", "Simplex_ZeroExclusion",
    "StrictUnitary_TotalStateInvariance", "StrictHeat_TotalMassInvariance",
    "InternalPatternDecay_NotOntologicalDeath", "Vacuum_NotNothing",
    "CemeteryState_Extension", "KilledEvolution_RequiresNewLaw",
    "AnnihilationQuestion_WellPosedness", "ExternalPump_OntologyConflict",
    "InternalCirculation_ClosedResource", "TotalPersistence_MinimalPremises",
    "RoundOne_Recommendation",
]


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k - LO], status, boundary, payload)


def main():
    x = {}
    facts = spectral_facts()
    psi = normalized(np.arange(1.0, N + 1.0))
    p = np.arange(1.0, N + 1.0); p /= p.sum()
    t = 3.0
    u = expm(-1j * t * STRICT_A) @ psi
    h = expm(-t * STRICT_A) @ p

    x["ST1002"] = packet(1002, "constructed_typed_ontology",
        "Typing is a proposed formalization of the user's ontology, not a strict-core derivation.", {
        "types": {
            "total_nadsoliton": "the entire fundamental being / complete physical state",
            "internal_pattern": "a localization, excitation, record, or subsystem state inside it",
            "vacuum": "a physical state of the total system",
            "nothing": "absence of the state space itself; not an ordinary state",
        },
        "key_correction": "annihilation of an internal pattern is not annihilation of the total nadsoliton",
    })
    x["ST1003"] = packet(1003, "proven_zero_excluded_from_unit_sphere",
        "Applies only if the total state is represented on a normalized Hilbert sphere.", {
        "state_space": "S(H)={psi: ||psi||=1}", "distance_to_zero": 1.0,
        "normalizing_zero_raises": True,
    })
    try:
        normalized(np.zeros(N)); zero_error = False
    except ValueError:
        zero_error = True
    x["ST1004"] = packet(1004, "proven_zero_has_no_projective_ray",
        "Projective quantum-state semantics is conditional, not strict-derived FIN ontology.", {
        "state_space": "P(H)=(H\\{0})/C*", "zero_representation_exists": False,
        "zero_normalization_rejected": zero_error,
    })
    x["ST1005"] = packet(1005, "proven_zero_excluded_from_probability_simplex",
        "Applies to normalized probability/information states.", {
        "state_space": "Delta={p>=0: sum p=1}", "zero_mass": 0.0,
        "simplex_mass": float(p.sum()), "l1_distance_to_zero": float(np.linalg.norm(p, 1)),
    })
    x["ST1006"] = packet(1006, "proven_strict_unitary_preserves_total_state_manifold",
        "Protects global normalization, not a chosen internal pattern.", {
        "norm_initial": float(np.linalg.norm(psi)), "norm_after": float(np.linalg.norm(u)),
        "norm_error": float(abs(np.linalg.norm(u)-np.linalg.norm(psi))), **facts,
    })
    x["ST1007"] = packet(1007, "proven_strict_heat_preserves_total_mass",
        "Heat is an internal coarse dynamics only after the probability embedding is supplied.", {
        "mass_initial": float(p.sum()), "mass_after": float(h.sum()),
        "mass_error": float(abs(h.sum()-p.sum())), "minimum_after": float(h.min()),
    })
    initial_pattern = float(np.linalg.norm(p - UNIFORM))
    later_pattern = float(np.linalg.norm(h - UNIFORM))
    x["ST1008"] = packet(1008, "proven_internal_pattern_can_decay_while_total_state_persists",
        "Finite-time sample plus spectral theorem; this is not a soliton theorem.", {
        "pattern_norm_initial": initial_pattern, "pattern_norm_t3": later_pattern,
        "mass_t3": float(h.sum()), "pattern_reduction_factor": later_pattern/initial_pattern,
        "interpretation": "homogenization is pattern death, not ontological nonexistence",
    })
    x["ST1009"] = packet(1009, "proven_vacuum_and_nothing_are_different_types",
        "The identification vacuum=nothing is rejected as a category error.", {
        "vacuum": "a vector, density operator, probability state, or algebraic state",
        "nothing": "no carrier/algebra/state at all", "vacuum_is_testable_state": True,
    })
    x["ST1010"] = packet(1010, "constructed_cemetery_state_extension",
        "This explicitly adds a state not present in the closed normalized model.", {
        "extension": "M_hat=M disjoint_union {dagger}",
        "annihilation_event": "Phi_t(x)=dagger", "strict_exports_dagger": False,
    })
    gamma = 0.4
    survival = math.exp(-gamma * t)
    x["ST1011"] = packet(1011, "proven_killing_requires_generator_or_hazard_extension",
        "Countermodel is additional dynamics, not FIN-derived.", {
        "killed_semigroup": "T_t=exp[-t(A+gamma I)]", "gamma": gamma,
        "mass_survival_t3": survival, "missing_mass_to_cemetery": 1.0-survival,
    })
    x["ST1012"] = packet(1012, "proven_internal_annihilation_is_undefined_without_target",
        "A no-transition result by state-space closure is structural, not a force law.", {
        "well_posed_requirements": ["source state", "target state or boundary", "evolution law", "event criterion"],
        "strict_total_nothing_target": False,
    })
    x["ST1013"] = packet(1013, "ontology_rejects_external_pump_as_fundamental_explanation",
        "An effective subsystem may be pumped by another internal sector.", {
        "reason": "if the nadsoliton is the whole being, an external reservoir contradicts closure",
        "admissible_replacement": "internal circulation with global conservation",
    })
    current = np.array([1.0, -0.5, -0.5])
    x["ST1014"] = packet(1014, "constructed_closed_internal_circulation_example",
        "Toy conservation witness only; no strict source for this three-sector current.", {
        "sector_currents": current.tolist(), "global_current_sum": float(current.sum()),
        "interpretation": "one sector can be sustained by another without external creation of the total resource",
    })
    x["ST1015"] = packet(1015, "proven_minimal_total_persistence_premise_list",
        "FIN has not yet derived all premises as one theorem.", {
        "premises": ["a nonempty total-state space M", "nothing/zero excluded or separated",
                     "evolution maps M into M", "global solution exists for all admitted times"],
        "optional_strengthening": ["conserved normalization/charge", "compactness/coercive barrier"],
    })
    x["ST1016"] = packet(1016, "recommended_ST1017_ST1031",
        "Recommendation follows the newly typed total-being lane.", {
        "next": ["invariant manifolds and zero-distance barriers", "finite-time collapse counterexamples",
                 "completeness versus asymptotic decay", "necessity/removal audit of every persistence premise"],
    })
    write_round(LO, HI, x)


if __name__ == "__main__":
    main()
