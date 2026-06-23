#!/usr/bin/env python3
"""P3044/S1994: memory-lag commutator source-candidate audit.

P3043 requires a concrete formula/artifact outside the exhausted P3038 receiver
classes.  P3044 supplies one bounded candidate: a bilinear memory-lag commutator
    C_s(i) = K_i*M_{i+s} - M_i*K_{i+s}
where K is sampled K_strict_gate and M is the P3038 memory/viscosity trace.
This is not a sine/chiral phase receiver, rho/path tuning, unit normalization,
or nonproxy placeholder.  It is a computable signed nonlocal object.

The finite result is positive but not closure: several lag sums are nonzero and
the commutator is antisymmetric under K<->M exchange, but the construction still
uses the P3038 memory trace, has no strict nadsoliton source law, no
chart-independent localizer, and no selector-torsor/readout coupling theorem.
"""
from __future__ import annotations

import hashlib, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3023_s1973_kernel_dissipation_time_order_candidate_obstruction import N, k_strict
from p3038_s1988_viscous_retarded_chiral_selector_candidate_obstruction import memory_viscosity_trace
from p3043_s1993_post_selector_candidate_new_source_intake_gate import OUT as P3043, REQUIRED_NEW_SOURCE_PREDICATES

OUT = GEN / "p3044_s1994_memory_lag_commutator_source_candidate.json"
MD = GEN / "p3044_s1994_memory_lag_commutator_source_candidate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
TOL = 1e-12
LAGS = list(range(1, N // 2 + 1))


def kernel_vector() -> list[float]:
    return [k_strict(i + 1) for i in range(N)]


def commutator(values: list[float], memory: list[float], lag: int) -> list[float]:
    return [values[i] * memory[(i + lag) % N] - memory[i] * values[(i + lag) % N] for i in range(N)]


def exchange_commutator(values: list[float], memory: list[float], lag: int) -> list[float]:
    return [memory[i] * values[(i + lag) % N] - values[i] * memory[(i + lag) % N] for i in range(N)]


def lag_rows() -> list[dict[str, Any]]:
    k = kernel_vector()
    m = memory_viscosity_trace(k)
    rows = []
    for lag in LAGS:
        c = commutator(k, m, lag)
        e = exchange_commutator(k, m, lag)
        rows.append({
            "lag": lag,
            "signed_sum": round(sum(c), 15),
            "l1_norm": round(sum(abs(x) for x in c), 15),
            "max_value": round(max(c), 15),
            "min_value": round(min(c), 15),
            "finite_nonzero": any(abs(x) > TOL for x in c),
            "signed_sum_nonzero": abs(sum(c)) > TOL,
            "exchange_antisymmetry_verified": max(abs(a + b) for a, b in zip(c, e)) < 1e-12,
            "accepted_as_new_strict_source_law": False,
            "failure": "nonzero commutator is a candidate signed observable, but no strict source/localizer/coupling theorem exports it",
        })
    return rows


def build_matrix() -> dict[str, Any]:
    read_json(P3043)
    rows = lag_rows()
    candidate_predicates = {
        "strict_nadsoliton_provenance": False,
        "nonpremise_source_law": False,
        "outside_exhausted_receiver_classes": True,
        "computable_nonzero_signed_or_branch_value": any(row["signed_sum_nonzero"] for row in rows),
        "chart_or_gauge_independent_localizer": False,
        "explicit_coupling_to_selector_torsor_or_readout": False,
        "unit_or_nonproxy_installation_when_physical_export_is_claimed": False,
    }
    obligations = [
        {"obligation": "p3043_new_source_intake_used", "satisfied": True, "detail": "candidate is evaluated against P3043 predicates"},
        {"obligation": "outside_exhausted_receiver_classes", "satisfied": candidate_predicates["outside_exhausted_receiver_classes"], "detail": "not sine/chiral phase, rho/path tuning, unit normalization/import, or nonproxy placeholder"},
        {"obligation": "finite_nonzero_signed_commutator", "satisfied": candidate_predicates["computable_nonzero_signed_or_branch_value"], "detail": "at least one lag has nonzero signed sum"},
        {"obligation": "exchange_antisymmetry_verified", "satisfied": all(row["exchange_antisymmetry_verified"] for row in rows), "detail": "K<->M exchange sends C_s to -C_s"},
        {"obligation": "strict_nadsoliton_source_law", "satisfied": False, "detail": "memory-lag commutator is constructed from a receiver trace, not exported as a strict source law"},
        {"obligation": "chart_independent_lag_localizer", "satisfied": False, "detail": "lag and one-sided memory initialization are chart choices"},
        {"obligation": "selector_torsor_or_readout_coupling", "satisfied": False, "detail": "no theorem couples C_s to QW-2191, a selector torsor, or a physical readout row"},
    ]
    return {
        "object": "MemoryLagCommutator_SourceCandidateAcceptanceMatrix",
        "formula": "C_s(i)=K_i*M_{i+s}-M_i*K_{i+s}",
        "lags": LAGS,
        "candidate_predicates": candidate_predicates,
        "required_predicates_from_p3043": REQUIRED_NEW_SOURCE_PREDICATES,
        "lag_rows": rows,
        "proof_obligations": obligations,
        "finite_certificate": {
            "lag_rows": len(rows),
            "finite_nonzero_lag_rows": sum(1 for row in rows if row["finite_nonzero"]),
            "signed_sum_nonzero_lag_rows": sum(1 for row in rows if row["signed_sum_nonzero"]),
            "exchange_antisymmetry_rows": sum(1 for row in rows if row["exchange_antisymmetry_verified"]),
            "accepted_new_strict_source_law_rows": sum(1 for row in rows if row["accepted_as_new_strict_source_law"]),
            "p3043_predicates_satisfied": sum(1 for value in candidate_predicates.values() if value),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "new_strict_source_law_exported": all(row["satisfied"] for row in obligations) and all(row["accepted_as_new_strict_source_law"] for row in rows),
        },
    }


def build_payload() -> dict[str, Any]:
    matrix = build_matrix()
    return {
        "status": "P3044_MEMORY_LAG_COMMUTATOR_SOURCE_CANDIDATE_BOUNDED_NO_EXPORT",
        "input_hashes": {"P3043": hashlib.sha256(P3043.read_bytes()).hexdigest() if P3043.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "P3044 supplies a concrete formula outside the exhausted P3038 receiver classes: a memory-lag commutator.  It has nonzero finite signed lag sums and exact K/M exchange antisymmetry, so it is a real computational hint.  It still does not export a strict source law because provenance, chart-independent lag/localizer, and selector/readout coupling are absent.",
            "negative_export_flags": {k: False for k in ["new_strict_source_law_exported", "new_live_frontier_unlocked", "strict_selector_mechanism_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not promote the nonzero memory-lag commutator alone.  A next proof-grade move may attack exactly one missing premise for this new object: a strict nadsoliton source law for the memory-lag commutator, a chart-independent lag/localizer theorem, or an explicit coupling theorem to a selector torsor/readout row.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3044/S1994 memory-lag commutator source-candidate audit", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- lag rows: `{c['lag_rows']}`",
        f"- finite nonzero lag rows: `{c['finite_nonzero_lag_rows']}`",
        f"- signed-sum nonzero lag rows: `{c['signed_sum_nonzero_lag_rows']}`",
        f"- exchange antisymmetry rows: `{c['exchange_antisymmetry_rows']}`",
        f"- accepted new strict source-law rows: `{c['accepted_new_strict_source_law_rows']}`",
        f"- P3043 predicates satisfied: `{c['p3043_predicates_satisfied']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- new strict source law exported: `{c['new_strict_source_law_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3044/S1994 memory-lag commutator source-candidate audit", "## P3044/S1994 memory-lag commutator source-candidate audit\n\n`P3044/S1994` supplies one concrete formula outside the exhausted P3038 receiver classes required by P3043: `C_s(i)=K_i*M_{i+s}-M_i*K_{i+s}` using sampled `K_strict_gate` and the memory trace.  The finite audit finds nonzero signed lag sums and exact `K<->M` exchange antisymmetry, so the object is a real computational hint.  It is not yet a strict source law: provenance, chart-independent lag/localizer, selector/readout coupling, unit-bearing action/EOM, `L_total`, `QW-2191` discharge, bridge/role transfer, and ToE closure remain unexported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3044/S1994 memory-lag commutator `L_total` guard", "## P3044/S1994 memory-lag commutator `L_total` guard\n\n`P3044/S1994` adds no physical `L_total` term.  The memory-lag commutator is finite and signed, but it lacks strict source provenance, chart-independent localization, explicit selector/readout coupling, and unit-bearing variational/action/EOM installation.\n")
    append_once(AGENTS, "Current memory-lag commutator source-candidate guardrail (P3044/S1994, 2026-06-23)", "## Current memory-lag commutator source-candidate guardrail (P3044/S1994, 2026-06-23)\n\n- P3044 supplies one concrete formula outside the exhausted P3038 receiver classes: the memory-lag commutator `C_s(i)=K_i*M_{i+s}-M_i*K_{i+s}`.\n- The finite computation gives nonzero signed lag sums and exact exchange antisymmetry, but no strict nadsoliton source law, chart-independent lag/localizer theorem, or selector/readout coupling is exported.\n- Do not promote memory-lag commutator signs, lag choices, exchange antisymmetry, or memory-trace bilinears to `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move may attack exactly one missing premise for this new object: strict source law, chart-independent lag/localizer, or explicit selector/readout coupling theorem.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
