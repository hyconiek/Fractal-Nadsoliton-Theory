#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from itertools import combinations, product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2529_s1479_strict_damping_numeric_subkey_rank_lattice_certificate.json"
MD = GEN / "p2529_s1479_strict_damping_numeric_subkey_rank_lattice_certificate.md"

SOURCE_FILES = {
    "P2528_PRIME_SLOPE_ANCHOR_EQUIVALENCE": GEN / "p2528_s1478_strict_damping_prime_slope_anchor_equivalence_certificate.json",
}

NODE_DOMAIN = list(range(1, 12))
PRIMES = [2, 3, 5, 7, 11]
STRICT_DELTA = 4.0 / 5.0
TOL = 1e-10
NUMERIC_SUBKEYS = ["M_multiplicative_character_law", "P_prime_log_proportionality", "A_single_prime_slope_anchor"]
SOURCE_KEYS = NUMERIC_SUBKEYS + ["O_m2_operator_signature_source"]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:40]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2529|S1479|numeric subkey rank lattice|strict damping numeric-source lattice|multiplicative prime anchor lattice|prime anchor source lattice",
        "intended_research_nonduplication": "numeric subkey rank lattice|prime anchor source lattice|multiplicative.*prime.*anchor|slope anchor source lattice|source-key lattice",
        "precursor_packets": "P2528|S1478|prime slope anchor equivalence|P2527|prime log proportionality slope line|P2526|finite monoid prime character nullity|P2516|dual-key source acceptance",
        "source_obligation_language": "multiplicative_character_law_source|prime_log_proportionality_source|slope_value_or_prime_anchor_source|m2_operator_signature_source|beta_eta_numeric_source",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def strict_vector() -> list[float]:
    return [STRICT_DELTA * math.log(d) for d in NODE_DOMAIN]


def multiplicative_pairs() -> list[tuple[int, int, int]]:
    return [(d, e, d * e) for d, e in product(NODE_DOMAIN, repeat=2) if d * e in NODE_DOMAIN]


def row_zero() -> list[float]:
    return [0.0 for _ in NODE_DOMAIN]


def multiplicative_constraints() -> list[tuple[list[float], float, str]]:
    rows = []
    for d, e, de in multiplicative_pairs():
        row = row_zero()
        row[de - 1] += 1.0
        row[d - 1] -= 1.0
        row[e - 1] -= 1.0
        rows.append((row, 0.0, f"y_{de}-y_{d}-y_{e}=0"))
    return rows


def prime_log_proportionality_constraints() -> list[tuple[list[float], float, str]]:
    rows = []
    log2 = math.log(2)
    for p in PRIMES[1:]:
        row = row_zero()
        row[p - 1] = log2
        row[2 - 1] = -math.log(p)
        rows.append((row, 0.0, f"log(2)*y_{p}-log({p})*y_2=0"))
    return rows


def single_prime_anchor_constraint(anchor_prime: int = 2) -> list[tuple[list[float], float, str]]:
    row = row_zero()
    row[anchor_prime - 1] = 1.0
    rhs = STRICT_DELTA * math.log(anchor_prime)
    return [(row, rhs, f"y_{anchor_prime}=(4/5)*log({anchor_prime})")]


def constraints_for_keys(keys: set[str]) -> list[tuple[list[float], float, str]]:
    rows: list[tuple[list[float], float, str]] = []
    if "M_multiplicative_character_law" in keys:
        rows.extend(multiplicative_constraints())
    if "P_prime_log_proportionality" in keys:
        rows.extend(prime_log_proportionality_constraints())
    if "A_single_prime_slope_anchor" in keys:
        rows.extend(single_prime_anchor_constraint(2))
    return rows


def matrix_rank(matrix: list[list[float]], tol: float = TOL) -> int:
    if not matrix:
        return 0
    work = [row[:] for row in matrix]
    row_count = len(work)
    col_count = len(work[0])
    rank = 0
    for col in range(col_count):
        pivot = max(range(rank, row_count), key=lambda r: abs(work[r][col]), default=rank)
        if rank >= row_count or abs(work[pivot][col]) <= tol:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][col]
        work[rank] = [value / pivot_value for value in work[rank]]
        for r in range(row_count):
            if r == rank or abs(work[r][col]) <= tol:
                continue
            factor = work[r][col]
            work[r] = [a - factor * b for a, b in zip(work[r], work[rank])]
        rank += 1
        if rank == row_count:
            break
    return rank


def max_strict_residual(rows: list[tuple[list[float], float, str]]) -> float:
    y = strict_vector()
    if not rows:
        return 0.0
    return max(abs(sum(coef * value for coef, value in zip(row, y)) - rhs) for row, rhs, _ in rows)


def numeric_subkey_rank_lattice() -> list[dict[str, Any]]:
    lattice = []
    for size in range(len(NUMERIC_SUBKEYS) + 1):
        for keys_tuple in combinations(NUMERIC_SUBKEYS, size):
            keys = set(keys_tuple)
            rows = constraints_for_keys(keys)
            coeffs = [row for row, _, _ in rows]
            rank = matrix_rank(coeffs)
            nullity = len(NODE_DOMAIN) - rank
            lattice.append({
                "active_numeric_subkeys": list(keys_tuple),
                "constraint_count": len(rows),
                "rank": rank,
                "nullity": nullity,
                "strict_target_satisfies_constraints": max_strict_residual(rows) < TOL,
                "max_strict_residual": max_strict_residual(rows),
                "unique_strict_numeric_target_on_node_space": rank == len(NODE_DOMAIN) and nullity == 0 and max_strict_residual(rows) < TOL,
                "proper_numeric_subkey_subset": set(keys_tuple) != set(NUMERIC_SUBKEYS),
            })
    return lattice


def source_key_lattice() -> list[dict[str, Any]]:
    rows = []
    for bits in product([False, True], repeat=len(SOURCE_KEYS)):
        active = [key for key, enabled in zip(SOURCE_KEYS, bits) if enabled]
        active_set = set(active)
        beta_eta = all(key in active_set for key in NUMERIC_SUBKEYS)
        strict_source = beta_eta and "O_m2_operator_signature_source" in active_set
        rows.append({
            "active_source_keys": active,
            "beta_eta_numeric_source_accepts": beta_eta,
            "strict_damping_beta_eta_source_accepts": strict_source,
            "proper_subset_rejected_for_strict_source": not strict_source and len(active) < len(SOURCE_KEYS),
        })
    return rows


def append_doc_sections() -> None:
    eq_section = """
`P2529/S1479` consolidates the P2525-P2528 numeric subkeys into a finite rank lattice on raw node variables `y_1,...,y_11`.  The three conditional numeric subkeys have rank/nullity chain: multiplicativity gives rank/nullity `6/5`, adding prime log-proportionality gives `10/1`, and adding one prime slope anchor gives `11/0`; every proper numeric subset leaves nullity positive.  The source-key lattice then records that `beta_eta_numeric_source` requires all three numeric subkeys, while strict damping source still additionally requires the independent `m=2` operator-signature source.
""".strip()
    lag_section = """
`P2529/S1479` records the current strict damping source normal form after the multiplicative/prime-character audits: `beta_eta_numeric_source = multiplicative_character_law_source AND prime_log_proportionality_source AND slope_value_or_prime_anchor_source`, and `strict_damping_beta_eta_source` additionally requires `m2_operator_signature_source`.  Since these are only conditional subkeys here, no nonlinear compression-flow source or role-bearing `L_total` term is licensed.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "P2529/S1479", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "P2529/S1479", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2528 = theorem(sources["P2528_PRIME_SLOPE_ANCHOR_EQUIVALENCE"], "strict_damping_prime_slope_anchor_equivalence_certificate")
    numeric_lattice = numeric_subkey_rank_lattice()
    source_lattice = source_key_lattice()
    full_numeric = next(row for row in numeric_lattice if set(row["active_numeric_subkeys"]) == set(NUMERIC_SUBKEYS))
    proper_numeric = [row for row in numeric_lattice if row["proper_numeric_subkey_subset"]]
    accepting_source_rows = [row for row in source_lattice if row["strict_damping_beta_eta_source_accepts"]]
    theorem_export = {
        "frontier_atom_under_attack": "strict_damping_numeric_source_subkey_lattice_after_prime_anchor_equivalence",
        "p2528_prime_anchor_equivalence_inherited": bool(p2528.get("conditional_prime_slope_anchor_equivalence_exported", False)),
        "node_variable_count": len(NODE_DOMAIN),
        "numeric_subkeys": NUMERIC_SUBKEYS,
        "source_keys": SOURCE_KEYS,
        "numeric_subkey_lattice_row_count": len(numeric_lattice),
        "source_key_lattice_row_count": len(source_lattice),
        "full_numeric_subkey_rank": full_numeric["rank"],
        "full_numeric_subkey_nullity": full_numeric["nullity"],
        "full_numeric_subkeys_conditionally_pin_unique_strict_target": full_numeric["unique_strict_numeric_target_on_node_space"],
        "every_proper_numeric_subkey_subset_leaves_positive_nullity": all(row["nullity"] > 0 for row in proper_numeric),
        "beta_eta_numeric_source_accepting_rows": sum(1 for row in source_lattice if row["beta_eta_numeric_source_accepts"]),
        "strict_damping_source_accepting_rows": len(accepting_source_rows),
        "strict_damping_source_accepts_only_all_four_keys": len(accepting_source_rows) == 1 and set(accepting_source_rows[0]["active_source_keys"]) == set(SOURCE_KEYS),
        "numeric_subkey_rank_lattice_exported": True,
        "multiplicative_character_law_source_exported": False,
        "prime_log_proportionality_source_exported": False,
        "slope_value_or_prime_anchor_source_exported": False,
        "beta_eta_numeric_source_exported": False,
        "m2_operator_signature_source_exported": False,
        "strict_damping_beta_eta_source_exported": False,
        "damping_compression_bridge_component_ready": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_claimed": False,
        "strict_damping_numeric_subkey_rank_lattice_certificate": {
            "numeric_subkey_rank_lattice": numeric_lattice,
            "source_key_lattice": source_lattice,
            "source_obligation_normal_form": "beta_eta_numeric_source = multiplicative_character_law_source AND prime_log_proportionality_source AND slope_value_or_prime_anchor_source; strict_damping_beta_eta_source additionally requires m2_operator_signature_source",
        },
    }
    gatekeepers = {
        "p2528_inherited": theorem_export["p2528_prime_anchor_equivalence_inherited"],
        "full_numeric_chain_unique_only_with_all_numeric_subkeys": theorem_export["full_numeric_subkeys_conditionally_pin_unique_strict_target"] and theorem_export["every_proper_numeric_subkey_subset_leaves_positive_nullity"],
        "strict_source_requires_all_four_keys": theorem_export["strict_damping_source_accepts_only_all_four_keys"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "multiplicative_character_law_source_exported",
            "prime_log_proportionality_source_exported",
            "slope_value_or_prime_anchor_source_exported",
            "beta_eta_numeric_source_exported",
            "m2_operator_signature_source_exported",
            "strict_damping_beta_eta_source_exported",
            "damping_compression_bridge_component_ready",
            "full_bridge_theorem_exported",
            "role_transfer_theorem_exported",
            "selector_closure_exported",
            "qw2191_discharged_by_this_certificate",
            "role_bearing_ltotal_exported",
            "toe_closure_claimed",
        ]),
    }
    return {
        "packet_id": "P2529",
        "stage_id": "S1479",
        "status": "STRICT_DAMPING_NUMERIC_SUBKEY_RANK_LATTICE_CERTIFICATE_CONDITIONAL_LATTICE_ONLY_NO_NUMERIC_SOURCE_NO_OPERATOR_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_numeric_subkey_rank_lattice_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_numeric_subkey_rank_lattice_certificate"]["theorem_export"]
    lines = [
        "# P2529/S1479 strict damping numeric subkey rank lattice certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2528 prime-anchor equivalence inherited: `{t['p2528_prime_anchor_equivalence_inherited']}`.",
        f"- Numeric subkey lattice rows: `{t['numeric_subkey_lattice_row_count']}`.",
        f"- Full numeric subkey rank/nullity: `{t['full_numeric_subkey_rank']}/{t['full_numeric_subkey_nullity']}`.",
        f"- Full numeric subkeys conditionally pin unique strict target: `{t['full_numeric_subkeys_conditionally_pin_unique_strict_target']}`.",
        f"- Every proper numeric subset leaves positive nullity: `{t['every_proper_numeric_subkey_subset_leaves_positive_nullity']}`.",
        f"- Source-key lattice rows: `{t['source_key_lattice_row_count']}`.",
        f"- Strict damping source accepts only all four keys: `{t['strict_damping_source_accepts_only_all_four_keys']}`.",
        f"- Beta/eta numeric source exported: `{t['beta_eta_numeric_source_exported']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only a conditional rank/source-key lattice. It does not source multiplicativity, prime-log proportionality, a slope/prime-anchor value, the m=2 operator signature, bridge completion, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_numeric_subkey_rank_lattice_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_numeric_subkey_rank_lattice_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
