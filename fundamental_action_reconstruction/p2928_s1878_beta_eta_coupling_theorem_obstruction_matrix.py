#!/usr/bin/env python3
"""P2928/S1878: beta/eta coupling-theorem obstruction matrix.

P2927 identified three absent obligations for a strict damping beta/eta source
packet: prime-log values, the delta/eta source law, and the beta/eta coupling
theorem.  P2928 attacks the coupling-theorem obligation alone.  It constructs
the formal finite coupling carrier from prime exponents and delta, verifies that
this carrier is multiplicative on the audited products, and then shows that the
carrier is only conditional: without strict L_p values and a strict delta source
it cannot become a sourced beta/eta theorem.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2927 = GEN / "p2927_s1877_strict_damping_beta_eta_source_packet_verifier.json"
OUT = GEN / "p2928_s1878_beta_eta_coupling_theorem_obstruction_matrix.json"
MD = GEN / "p2928_s1878_beta_eta_coupling_theorem_obstruction_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

NODES = list(range(1, 12))
PRIMES = [2, 3, 5, 7, 11]


def factor_vector(n: int) -> list[int]:
    remaining = n
    vector = []
    for prime in PRIMES:
        exponent = 0
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        vector.append(exponent)
    if remaining != 1:
        raise ValueError(f"node {n} has unaudited factor {remaining}")
    return vector


def formal_coupling_rows() -> list[dict[str, Any]]:
    rows = []
    for d in NODES:
        vector = factor_vector(d)
        rows.append({
            "node": d,
            "prime_exponent_vector": vector,
            "formal_log_tail": "delta*(" + (" + ".join(f"{exp}*L_{p}" for p, exp in zip(PRIMES, vector) if exp) or "0") + ")",
            "formal_tail_factor": "exp(delta*sum_p v_p(d)*L_p)",
            "conditional_beta_eta_reading": "beta=1 from y_1=0 and eta=1+delta only if L_p and delta are strictly sourced",
        })
    return rows


def product_coupling_rows() -> list[dict[str, Any]]:
    rows = []
    for d in NODES:
        for e in NODES:
            if d * e <= NODES[-1]:
                vd = factor_vector(d)
                ve = factor_vector(e)
                vde = factor_vector(d * e)
                defect = [a + b - c for a, b, c in zip(vd, ve, vde)]
                rows.append({
                    "d": d,
                    "e": e,
                    "de": d * e,
                    "coupling_exponent_defect_vector": defect,
                    "passes_formal_multiplicative_coupling": defect == [0, 0, 0, 0, 0],
                })
    return rows


def theorem_schema() -> dict[str, Any]:
    return {
        "name": "Strict_Beta_Eta_Coupling_Theorem",
        "premises": [
            "Strict_Prime_Log_Value_Source_Law exports nonzero L_p values",
            "Strict_Damping_Delta_Source_Law exports delta=4/5 and eta=9/5",
            "P2601/P2923 carrier readiness supplies y_1=0 and prime-exponent product additivity",
        ],
        "conclusion": "the strict damping/compression tail is beta-normalized and eta-sourced rather than imported from fitted or legacy parameters",
        "nonconclusion_without_premises": "formal coupling readiness alone does not export strict beta/eta, L_total, bridge, role transfer, or ToE",
    }


def candidate_coupling_theorems() -> list[dict[str, Any]]:
    return [
        {
            "candidate": "formal_exponent_coupling_carrier",
            "formal_coupling_passes": True,
            "strict_Lp_values_present": False,
            "strict_delta_source_present": False,
            "accepted_as_coupling_theorem": False,
            "failure": "the carrier is algebraically correct but conditional on missing sourced values and slope",
        },
        {
            "candidate": "external_strict_gate_eta_1_8_beta_1",
            "formal_coupling_passes": True,
            "strict_Lp_values_present": False,
            "strict_delta_source_present": False,
            "accepted_as_coupling_theorem": False,
            "failure": "imports the later operational gate tuple rather than deriving the source law inside FAR",
        },
        {
            "candidate": "legacy_beta_tors_bridge_reuse",
            "formal_coupling_passes": False,
            "strict_Lp_values_present": False,
            "strict_delta_source_present": False,
            "accepted_as_coupling_theorem": False,
            "failure": "legacy linear torsion damping cannot be silently substituted for strict nonlinear beta/eta coupling",
        },
        {
            "candidate": "P2927_packet_verifier_name",
            "formal_coupling_passes": False,
            "strict_Lp_values_present": False,
            "strict_delta_source_present": False,
            "accepted_as_coupling_theorem": False,
            "failure": "the verifier names acceptance obligations but does not instantiate the coupling theorem",
        },
        {
            "candidate": "future_Strict_Beta_Eta_Coupling_Theorem",
            "formal_coupling_passes": False,
            "strict_Lp_values_present": False,
            "strict_delta_source_present": False,
            "accepted_as_coupling_theorem": False,
            "failure": "admissible future theorem shape, but no current proof artifact is supplied",
        },
    ]


def build_payload(p2927: dict[str, Any]) -> dict[str, Any]:
    coupling_rows = formal_coupling_rows()
    product_rows = product_coupling_rows()
    candidates = candidate_coupling_theorems()
    accepted = [candidate for candidate in candidates if candidate["accepted_as_coupling_theorem"]]
    return {
        "status": "P2928_BETA_ETA_COUPLING_THEOREM_OBSTRUCTION_MATRIX_CONDITIONAL_ONLY",
        "input_hashes": {"P2927": hashlib.sha256(P2927.read_bytes()).hexdigest() if P2927.exists() else None},
        "constructed_theoretical_objects": {
            "missing_theorem_schema": theorem_schema(),
            "formal_coupling_rows": coupling_rows,
            "product_coupling_rows": product_rows,
            "candidate_coupling_theorems": candidates,
        },
        "coupling_matrix": {
            "node_count": len(coupling_rows),
            "product_pair_count_de_le_11": len(product_rows),
            "formal_product_coupling_failures": sum(1 for row in product_rows if not row["passes_formal_multiplicative_coupling"]),
            "candidate_coupling_theorem_count": len(candidates),
            "accepted_coupling_theorem_count": len(accepted),
            "p2927_packet_absent_obligations_inherited": p2927.get("acceptance_matrix", {}).get("strict_damping_beta_eta_source_packet_exported") is False,
            "formal_coupling_carrier_exported": True,
            "strict_beta_eta_coupling_theorem_exported": False,
            "strict_damping_beta_eta_source_packet_exported": False,
        },
        "decision": {
            "positive_witnesses": {
                "formal_coupling_carrier_constructed": True,
                "formal_product_coupling_verified": all(row["passes_formal_multiplicative_coupling"] for row in product_rows),
                "missing_coupling_theorem_schema_exported": True,
            },
            "negative_export_flags": {
                "strict_prime_log_value_source_exported": False,
                "strict_delta_source_law_exported": False,
                "strict_beta_eta_coupling_theorem_exported": False,
                "strict_damping_beta_eta_source_packet_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_hamiltonian_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2928 constructs the formal beta/eta coupling carrier and verifies zero multiplicative coupling defects on all audited products.  This is conditional readiness only: without strict L_p values and a strict delta=4/5/eta=9/5 source law, no strict coupling theorem or damping source packet is exported.",
            "next_honest_step": "A next admissible move must provide a concrete strict source artifact for L_p values or delta/eta, or a proof of the Strict_Beta_Eta_Coupling_Theorem using such sources.  Otherwise preserve the P2927/P2928 no-new-live-frontier certificate and pivot to a fresh state-map object outside this closed damping replay.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    matrix = payload["coupling_matrix"]
    lines = [
        "# P2928/S1878 beta/eta coupling theorem obstruction matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Coupling matrix",
        f"- node count: `{matrix['node_count']}`",
        f"- product pairs de<=11: `{matrix['product_pair_count_de_le_11']}`",
        f"- formal product coupling failures: `{matrix['formal_product_coupling_failures']}`",
        f"- candidate coupling theorems: `{matrix['candidate_coupling_theorem_count']}`",
        f"- accepted coupling theorems: `{matrix['accepted_coupling_theorem_count']}`",
        f"- strict beta/eta coupling theorem exported: `{matrix['strict_beta_eta_coupling_theorem_exported']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2927))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2928/S1878 beta/eta coupling theorem obstruction matrix", "## P2928/S1878 beta/eta coupling theorem obstruction matrix\n\n`P2928/S1878` attacks the P2927 coupling-theorem obligation directly.  It constructs the formal carrier `exp(delta*sum_p v_p(d)*L_p)` and verifies zero formal multiplicative coupling defects on all `29` audited products `d*e<=11`.  This is conditional readiness only: current artifacts still lack strict `L_p` values, a strict `delta=4/5`/`eta=9/5` source law, and therefore a strict beta/eta coupling theorem.  No strict damping source packet, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2928/S1878 beta/eta coupling theorem `L_total` guard", "## P2928/S1878 beta/eta coupling theorem `L_total` guard\n\n`P2928/S1878` supplies only a conditional formal coupling carrier for the damping term.  Without sourced `L_p` values and sourced `delta=4/5`/`eta=9/5`, the carrier cannot be promoted to a role-bearing damping contribution in nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current beta/eta coupling theorem obstruction guardrail (P2928/S1878, 2026-06-19)", "## Current beta/eta coupling theorem obstruction guardrail (P2928/S1878, 2026-06-19)\n\n- P2928 attacks the P2927 coupling-theorem obligation by constructing the formal carrier `exp(delta*sum_p v_p(d)*L_p)`.\n- The formal carrier has zero product-coupling defects on all audited products `d*e<=11`, but this is conditional readiness only because strict `L_p` values and strict `delta=4/5`/`eta=9/5` remain unsourced.\n- Do not promote the formal carrier to strict damping `beta/eta`, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n- A next admissible move must provide a concrete strict source artifact for one missing input or pivot to a fresh state-map object outside this closed damping replay.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
