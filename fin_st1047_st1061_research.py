#!/usr/bin/env python3
"""FIN ST1047--ST1061: the internal observer and operational meaning."""
import numpy as np
from scipy.linalg import expm

from fin_total_nadsoliton_common import N, STRICT_A, UNIFORM, normalized, write_packet, write_round


LO, HI = 1047, 1061
NAMES = [
    "InternalObserver_TypedInclusion", "Record_IsInternalPattern",
    "TotalAnnihilation_NoPostEventRecord", "VacuumOperationallyDistinct",
    "CoarseObserver_PatternBlindness", "TraceOut_LocalDeath",
    "DualDynamics_SameGeneratorDifferentQuestions", "UnitaryHeat_NotPhysicalEquivalence",
    "WaveClock_Conditional", "TotalState_ObserverIncluded",
    "NearBoundary_PrecursorOnly", "SelfMeasurement_Limitation",
    "OperationalAnnihilation_Definitions", "CurrentFIN_OperationalBoundary",
    "RoundFour_Recommendation",
]


def packet(k, status, boundary, payload):
    return write_packet(k, NAMES[k - LO], status, boundary, payload)


def main():
    x = {}
    x["ST1047"] = packet(1047, "constructed_internal_observer_inclusion",
        "Operational structure remains an added typed layer, as already required by FIN measurement audits.", {
        "typing": "observer O, apparatus M, environment E, and record R are subsystems/patterns of total state Omega",
        "external_observer": False,
    })
    x["ST1048"] = packet(1048, "proven_record_persistence_is_weaker_than_total_persistence",
        "A record is an encoded correlation, not the whole being.", {
        "record": "correlation between memory and measured sector", "record_can_erase": True,
        "global_state_can_persist": True,
    })
    x["ST1049"] = packet(1049, "proven_no_internal_post_annihilation_record",
        "Logical/operational no-go, not evidence that annihilation cannot occur.", {
        "premise": "total annihilation removes observer, apparatus, clock, and record",
        "conclusion": "no internal record can be produced after the event",
        "empirical_status": "direct post-event verification from inside is impossible",
    })
    x["ST1050"] = packet(1050, "proven_vacuum_is_operationally_distinguishable_from_nothing",
        "A physical vacuum requires an observable algebra and state; nothing has neither.", {
        "vacuum_has_correlators": True, "vacuum_can_be_prepared_or_approached": True,
        "nothing_has_measurement_statistics": False,
    })
    p = np.zeros(N); p[0] = 1.0
    pt = expm(-4.0*STRICT_A) @ p
    coarse = float(pt.sum())
    fine_contrast = float(np.linalg.norm(pt-UNIFORM))
    x["ST1051"] = packet(1051, "proven_coarse_observer_can_miss_surviving_fine_pattern",
        "Specific coarse statistic is total mass; instrument choice is not strict-generated.", {
        "coarse_mass": coarse, "fine_contrast_t4": fine_contrast,
        "uniform_mass": float(UNIFORM.sum()), "coarse_records_equal": abs(coarse-1.0)<1e-12,
    })
    x["ST1052"] = packet(1052, "proven_local_purity_loss_does_not_imply_global_information_loss",
        "Subsystem factorization is conditional.", {
        "example": "Bell state: global pure, each reduced state maximally mixed",
        "local_purity": 0.5, "global_purity": 1.0,
    })
    psi = normalized(np.arange(1.0,N+1.0))
    unit = expm(-1j*STRICT_A)@psi
    prob = np.abs(psi)**2; heat = expm(-STRICT_A)@prob
    x["ST1053"] = packet(1053, "proven_dual_channels_answer_distinct_persistence_questions",
        "Both are functional calculi of A, but their state types differ.", {
        "unitary_norm": float(np.linalg.norm(unit)), "heat_mass": float(heat.sum()),
        "unitary_preserves": "amplitude norm/coherent reversibility",
        "heat_preserves": "mass while damping nonzero spectral modes",
    })
    x["ST1054"] = packet(1054, "proven_same_generator_does_not_make_unitary_and_heat_physically_equivalent",
        "No Wick rotation, bath, preparation, clock, or measurement bridge is strict-derived.", {
        "mathematical_relation": "same spectral projectors, multipliers exp(-it lambda) vs exp(-t lambda)",
        "physical_equivalence": False,
    })
    x["ST1055"] = packet(1055, "conditional_frequency_clock_does_not_prove_total_time",
        "A periodic internal mode can serve as a relational clock only with a state, phase readout, and monotonic/counting convention.", {
        "clock_candidate": "relative phase exp(-it(lambda_j-lambda_k))",
        "requires": ["two populated modes", "phase-coherent instrument", "calibration scale"],
        "strict_dimensional_seconds": False,
    })
    x["ST1056"] = packet(1056, "proven_total_state_description_must_include_observer_backreaction",
        "Conceptual consistency condition; no complete FIN measurement dynamics is exported.", {
        "closed_description": "Omega includes system+observer+apparatus+environment+record",
        "ignored_backreaction_is_effective": True,
    })
    x["ST1057"] = packet(1057, "proven_only_pre_event_precursors_are_internally_recordable",
        "Precursors are model-dependent and cannot certify literal nonexistence afterward.", {
        "recordable": ["norm leakage in a subsystem", "gap closure", "loss of coherence", "approach to a model boundary"],
        "not_recordable_after_total_event": ["survival probability", "post-event state", "external timestamp"],
    })
    x["ST1058"] = packet(1058, "proven_internal_observer_cannot_operationally_compare_with_absence_of_all_observers",
        "An epistemic no-go; it does not settle ontology.", {
        "comparison_requires": "a common record space containing both outcomes",
        "total_nothing_supplies_common_record_space": False,
    })
    x["ST1059"] = packet(1059, "constructed_operational_annihilation_hierarchy",
        "Only the first three are internally testable in principle.", {
        "testable": ["pattern disappearance", "record/coherence loss", "subsystem vacuum transition"],
        "not_directly_internal_testable": ["annihilation of the complete carrier/state space"],
    })
    x["ST1060"] = packet(1060, "blocked_no_complete_FIN_operational_total_state",
        "State, clock, preparation, instrument, environment, apparatus, and record remain unclosed.", {
        "complete_operational_model": False, "laboratory_record": False,
        "total_annihilation_test": False,
    })
    x["ST1061"] = packet(1061, "recommended_ST1062_ST1076",
        "Next round attempts to destroy every proposed protection package.", {
        "next": ["absorbing-zero extension", "stochastic killing", "finite-time collapse", "topological gap closure",
                 "normalization-gauge critique", "observational equivalence and premise-removal matrix"],
    })
    write_round(LO, HI, x)


if __name__ == "__main__":
    main()
