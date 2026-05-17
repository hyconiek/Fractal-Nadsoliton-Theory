#!/usr/bin/env python3
"""P1956 S906 strict gauge-gauge physical projector certificate.

This executor exports the exact, local on-shell transverse polarization
projector for the two gauge bosons in the graviton -> gauge_gauge Cutkosky
channel.  It is deliberately narrower than a BRST theorem: it proves the
polarization-sum algebra in one flat reference cut kinematics, while keeping
BRST nilpotency, ghost cancellation, and cohomology as open requirements.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1956_s906_strict_gauge_gauge_physical_projector_certificate.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest(obj: object) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


ETA = sp.diag(-1, 1, 1, 1)
ETA_CONTRA = sp.diag(-1, 1, 1, 1)


def vec(entries: tuple[int, int, int, int]) -> sp.Matrix:
    return sp.Matrix([sp.Integer(item) for item in entries])


def dot(u: sp.Matrix, v: sp.Matrix) -> sp.Expr:
    return sp.simplify((u.T * ETA * v)[0])


def lower(v: sp.Matrix) -> sp.Matrix:
    return ETA * v


def outer(u: sp.Matrix, v: sp.Matrix) -> sp.Matrix:
    return u * v.T


def axial_projector(k: sp.Matrix, n: sp.Matrix) -> sp.Matrix:
    """Return P^{mu nu} for signature (-,+,+,+).

    For k^2=0 and k.n != 0:
      P^{mu nu}=eta^{mu nu}
        -(k^mu n^nu+n^mu k^nu)/(k.n)
        +n^2 k^mu k^nu/(k.n)^2.
    In the exported reference frame n is lightlike, so the last term vanishes,
    but it is kept to make the exact convention explicit.
    """

    kd = dot(k, n)
    n2 = dot(n, n)
    return sp.simplify(
        ETA_CONTRA
        - (outer(k, n) + outer(n, k)) / kd
        + n2 * outer(k, k) / kd**2
    )


def polarization_sum(pols: list[sp.Matrix]) -> sp.Matrix:
    acc = sp.zeros(4, 4)
    for eps in pols:
        acc += outer(eps, eps)
    return sp.simplify(acc)


def mat_zero(m: sp.Matrix) -> bool:
    return all(sp.simplify(item) == 0 for item in list(m))


def matrix_to_nested_strings(m: sp.Matrix) -> list[list[str]]:
    return [[str(sp.simplify(m[i, j])) for j in range(m.cols)] for i in range(m.rows)]


def vector_to_strings(v: sp.Matrix) -> list[str]:
    return [str(sp.simplify(v[i])) for i in range(v.rows)]


def projector_checks(label: str, k: sp.Matrix, n: sp.Matrix, pols: list[sp.Matrix]) -> dict[str, Any]:
    projector = axial_projector(k, n)
    pol_sum = polarization_sum(pols)
    completeness_residual = sp.simplify(pol_sum - projector)
    projector_mixed = sp.simplify(projector * ETA)
    idempotence_residual = sp.simplify(projector_mixed * projector_mixed - projector_mixed)
    k_transverse_residual = sp.simplify(projector * lower(k))
    n_transverse_residual = sp.simplify(projector * lower(n))

    pol_norms = [dot(eps, eps) for eps in pols]
    pol_cross = [dot(pols[i], pols[j]) for i in range(len(pols)) for j in range(i + 1, len(pols))]
    pol_k_transverse = [dot(k, eps) for eps in pols]
    pol_n_transverse = [dot(n, eps) for eps in pols]

    return {
        "label": label,
        "k": vector_to_strings(k),
        "n_reference": vector_to_strings(n),
        "polarizations": [vector_to_strings(eps) for eps in pols],
        "kinematic_scalars": {
            "k_squared": str(dot(k, k)),
            "n_squared": str(dot(n, n)),
            "k_dot_n": str(dot(k, n)),
        },
        "projector_contravariant": matrix_to_nested_strings(projector),
        "projector_mixed": matrix_to_nested_strings(projector_mixed),
        "polarization_sum": matrix_to_nested_strings(pol_sum),
        "machine_checks": {
            "k_is_lightlike": dot(k, k) == 0,
            "reference_is_lightlike": dot(n, n) == 0,
            "reference_not_collinear": dot(k, n) != 0,
            "polarization_norms": [str(item) for item in pol_norms],
            "polarization_cross_terms": [str(item) for item in pol_cross],
            "polarizations_transverse_to_k": [str(item) for item in pol_k_transverse],
            "polarizations_transverse_to_n": [str(item) for item in pol_n_transverse],
            "completeness_residual_zero": mat_zero(completeness_residual),
            "projector_k_transverse_zero": mat_zero(k_transverse_residual),
            "projector_n_transverse_zero": mat_zero(n_transverse_residual),
            "idempotence_residual_zero": mat_zero(idempotence_residual),
            "projector_trace_mixed": str(sp.trace(projector_mixed)),
            "projector_rank_mixed": int(projector_mixed.rank()),
        },
        "residuals": {
            "polarization_sum_minus_axial_projector": matrix_to_nested_strings(completeness_residual),
            "projector_on_k_lower": vector_to_strings(k_transverse_residual),
            "projector_on_n_lower": vector_to_strings(n_transverse_residual),
            "mixed_idempotence_residual": matrix_to_nested_strings(idempotence_residual),
        },
    }


def main() -> None:
    p1767 = load("p1767_s717_strict_bianchi_ward_to_brst_cutkosky_gate_sequencing_checkpoint.json")
    p1801 = load("p1801_s751_strict_brst_nilpotency_witness_intake_and_gate_link_checkpoint.json")
    p1802 = load("p1802_s752_strict_cutkosky_unitarity_witness_intake_and_tg3_gate_link_checkpoint.json")
    p1954 = load("p1954_s904_strict_dressed_amplitude_nonavailability_theorem.json")
    p1955 = load("p1955_s905_strict_minimal_hAA_vertex_export.json")

    # Center-of-mass massless two-body cut frame, E normalized to one.
    k1 = vec((1, 0, 0, 1))
    k2 = vec((1, 0, 0, -1))
    n1 = vec((1, 0, 0, -1))
    n2 = vec((1, 0, 0, 1))
    eps_x = vec((0, 1, 0, 0))
    eps_y = vec((0, 0, 1, 0))

    check_1 = projector_checks("gauge_boson_1", k1, n1, [eps_x, eps_y])
    check_2 = projector_checks("gauge_boson_2", k2, n2, [eps_x, eps_y])

    p1_mixed = sp.Matrix([[sp.Integer(item) for item in row] for row in check_1["projector_mixed"]])
    p2_mixed = sp.Matrix([[sp.Integer(item) for item in row] for row in check_2["projector_mixed"]])
    two_particle_projector = sp.kronecker_product(p1_mixed, p2_mixed)
    two_particle_idempotence = sp.simplify(
        two_particle_projector * two_particle_projector - two_particle_projector
    )
    two_particle_rank = int(two_particle_projector.rank())
    two_particle_trace = sp.trace(two_particle_projector)

    local_checks = [
        check_1["machine_checks"]["k_is_lightlike"],
        check_2["machine_checks"]["k_is_lightlike"],
        check_1["machine_checks"]["reference_is_lightlike"],
        check_2["machine_checks"]["reference_is_lightlike"],
        check_1["machine_checks"]["reference_not_collinear"],
        check_2["machine_checks"]["reference_not_collinear"],
        check_1["machine_checks"]["completeness_residual_zero"],
        check_2["machine_checks"]["completeness_residual_zero"],
        check_1["machine_checks"]["projector_k_transverse_zero"],
        check_2["machine_checks"]["projector_k_transverse_zero"],
        check_1["machine_checks"]["projector_n_transverse_zero"],
        check_2["machine_checks"]["projector_n_transverse_zero"],
        check_1["machine_checks"]["idempotence_residual_zero"],
        check_2["machine_checks"]["idempotence_residual_zero"],
        two_particle_rank == 4,
        two_particle_trace == 4,
        mat_zero(two_particle_idempotence),
    ]

    out = {
        "packet_id": "P1956",
        "stage_id": "S906",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "local_verdict": "GAUGE_GAUGE_TRANSVERSE_POLARIZATION_PROJECTOR_EXPORTED__BRST_COHOMOLOGY_STILL_OPEN",
        "route": "strict_only",
        "legacy_bridge_used": False,
        "higher_reasoning_used": True,
        "pre_execution_grep_summary": {
            "english_terms": [
                "BRST",
                "physical-state projector",
                "polarization projector",
                "transverse polarization",
                "ghost sector",
                "ghost cancellation",
                "gauge_gauge",
                "Cutkosky",
            ],
            "polish_terms": [
                "projektor",
                "polaryzacja",
                "poprzeczny",
                "sektor duchow",
                "kasowanie duchow",
                "unitarnosc",
                "kohomologia",
            ],
            "key_existing_sources_found": [
                "P1767: BRST/Cutkosky gate sequencing requires G_BW PASS, ghost_sector_nonproxy_export, and BV_BRST_operator_map",
                "P1801: BRST nilpotency intake template, not a channel projector",
                "P1802: Cutkosky unitarity intake template, not a computed projector",
                "P1854: B1 BRST cochain and first Cutkosky channel proxy, still pre-theorem",
                "P1954: dressed amplitude nonavailability theorem lists missing M2 BRST projector",
                "P1955: minimal tree-level hAA vertex exported, leaving BRST projection open",
            ],
            "negative_search_result": (
                "No existing BRSTPhysicalProjector_gauge_gauge_strict_B1_v1 or "
                "channel-level ghost-subtracted BRST cohomology projector was found."
            ),
        },
        "depends_on": {
            "p1767_present": "gate_sequencing_contract" in p1767,
            "p1801_present": "tg2_pass_requirements" in p1801,
            "p1802_present": "tg3_pass_requirements" in p1802,
            "p1954_present": "formal_nonavailability_theorem" in p1954,
            "p1955_present": "V_hAA_tensor_strict_B1_v1" in p1955,
        },
        "input_hashes": {
            "p1767_sha256": digest(p1767),
            "p1801_sha256": digest(p1801),
            "p1802_sha256": digest(p1802),
            "p1954_sha256": digest(p1954),
            "p1955_sha256": digest(p1955),
        },
        "kinematic_scope": {
            "metric_signature": "eta_mu_nu = diag(-1,1,1,1)",
            "frame": "flat two-body massless cut frame, E=1",
            "channel": "graviton -> gauge_gauge external cut states",
            "kernel_lane": "K_strict operational lane only; no legacy kernel role transfer",
            "scheme_scope": "external on-shell polarization sum compatible with the P1955 minimal tree-level hAA seed; no loop dressing",
        },
        "projector_formula": {
            "name": "TransverseAxialPolarizationProjector_v1",
            "formula": (
                "P^{mu nu}(k,n)=eta^{mu nu}-(k^mu n^nu+n^mu k^nu)/(k.n)"
                "+(n.n)k^mu k^nu/(k.n)^2"
            ),
            "reference_choice": "n_i is the opposite lightlike cut momentum, so n_i^2=0 and k_i.n_i=-2",
            "polarization_sum_identity": "sum_{lambda=x,y} eps_lambda^mu eps_lambda^nu = P^{mu nu}(k,n)",
        },
        "single_particle_projector_checks": [check_1, check_2],
        "two_particle_gauge_gauge_projector": {
            "construction": "P_gauge_gauge = P_gauge_boson_1 mixed kron P_gauge_boson_2 mixed",
            "rank": two_particle_rank,
            "trace": str(two_particle_trace),
            "physical_state_count": 4,
            "idempotence_residual_zero": mat_zero(two_particle_idempotence),
        },
        "machine_check_summary": {
            "all_local_projector_checks_zero": all(local_checks),
            "single_particle_projector_ranks": [
                check_1["machine_checks"]["projector_rank_mixed"],
                check_2["machine_checks"]["projector_rank_mixed"],
            ],
            "two_particle_rank_expected_4": two_particle_rank == 4,
            "two_particle_trace_expected_4": two_particle_trace == 4,
            "two_particle_idempotence_zero": mat_zero(two_particle_idempotence),
        },
        "p1954_requirement_update": {
            "M1_4D_hAA_vertex": "DISCHARGED_FOR_MINIMAL_TREE_LEVEL_FIELD_STRENGTH_VERTEX_BY_P1955",
            "M2_BRST_projector": "DISCHARGED_ONLY_FOR_LOCAL_ONSHELL_TRANSVERSE_POLARIZATION_SUM",
            "remaining_M2_limitations": [
                "no BV/BRST charge Q exported",
                "no Q^2=0 nilpotency proof for the interacting strict sector",
                "no ghost-sector nonproxy action and cancellation trace",
                "no proof that the full cut sum is restricted to BRST cohomology",
                "no global Lorentz-covariant reference-vector family over all cut momenta",
            ],
            "M3_dressed_propagator_residues": "OPEN",
            "M4_DiscM_CutSum_same_scheme": "OPEN",
            "M5_integral_reductions": "OPEN",
        },
        "gatekeeper_status": {
            "TG2_BRST_GLOBAL_NILPOTENCY": "NOT_PROMOTED",
            "TG3_CUTKOSKY_GLOBAL_UNITARITY": "NOT_PROMOTED",
            "reason": "P1956 supplies only external transverse polarization algebra; P1767/P1801/P1802 still require BRST charge, ghost sector, nilpotency, cohomology, and optical theorem matching.",
        },
        "false_pass_guard": (
            "This certificate is a local on-shell physical polarization projector. "
            "It is not a BRST cohomology theorem, not a ghost-cancellation theorem, "
            "not a dressed Cutkosky equality, and not a global UR_link theorem."
        ),
        "next_honest_step": (
            "Build P1957 with high reasoning: either construct the strict BV/BRST charge, "
            "ghost-sector nonproxy action, Q^2=0 nilpotency check, and gauge_gauge ghost-cancellation trace, "
            "or export a formal obstruction proving those data are not yet available."
        ),
        "lay_explanation": (
            "Wyodrebnilismy dokladnie dwie fizyczne polaryzacje dla kazdego z dwoch bozonow cechowania "
            "i pokazalismy rachunkowo, ze daja cztery dodatnie stany koncowe. To potrzebny klocek "
            "unitarnosci, ale jeszcze nie pelny test BRST: nadal trzeba pokazac, ze niefizyczne skladowe "
            "i duchy naprawde sie kasuja w calej teorii."
        ),
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
