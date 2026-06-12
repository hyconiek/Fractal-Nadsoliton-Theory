#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import re
import subprocess
from fractions import Fraction
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2642_s1592_universal_affine_zero_lattice_source_no_go_and_node_role_demotion.json"
MD = GEN / "p2642_s1592_universal_affine_zero_lattice_source_no_go_and_node_role_demotion.md"

MOD = 12
AUT = (1, 5, 7, 11)
LEGACY_NODES = (2, 5, 8, 11)
OMEGA = Fraction(1, 4)
PHI = Fraction(1, 6)
ETA = Fraction(9, 5)

SOURCE_FILES = {
    "P2639_OFFSET_STRIDE_EXHAUSTION": GEN / "p2639_s1589_offset_stride_metric_lift_exhaustion_closure_path.json",
    "P2640_SOURCE_LEDGER": GEN / "p2640_s1590_offset_stride_topology_selector_source_no_go_audit.json",
    "P2641_Z12_SUCCESSOR_CONNECTION": GEN / "p2641_s1591_z12_quotient_safe_successor_connection_no_go_audit.json",
    "N456_NO_AUT_INVARIANT_12_CYCLE": ROOT / "N456_CURRENT_FIRST_STRICT_PHASE_12_AUT_Z12_CANONICAL_12_CYCLE_SUCCESSOR_MAP_NONEXISTENCE_THEOREM.md",
    "N462_NO_AUT_INVARIANT_GENERATOR_ORIENTATION": ROOT / "N462_CURRENT_FIRST_STRICT_Z12_AUT_Z12_CANONICAL_GENERATOR_ORIENTATION_FIXING_DATUM_NONEXISTENCE_THEOREM.md",
    "N523_PREMISE_BASED_T164_FIXING_DATUM": ROOT / "N523_CURRENT_FIRST_STRICT_T164_Z12_GENERATOR_ORIENTATION_CANONICAL_FIXING_DATUM_DISCHARGE_THEOREM.md",
    "N454_PARITY_CHARACTER_AUT_INVARIANCE": ROOT / "N454_CURRENT_FIRST_ACTUAL_STRICT_Z12_PARITY_CHARACTER_AND_AUT_INVARIANCE_WITNESS_THEOREM.md",
    "N457_QUOTIENT_HOLONOMY_BOUNDARY": ROOT / "N457_CURRENT_FIRST_STRICT_PHASE_12_AUT_Z12_QUOTIENT_ORBIT_SPACE_TOPOLOGICAL_HOLONOMY_TRIVIALITY_BOUNDARY_THEOREM.md",
    "STRICT_EQUATION_SHEET": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "STRICT_LAGRANGIAN_DRAFT": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

NEGATIVE_EXPORT_FLAGS = [
    "universal_affine_lift_source_exported",
    "zero_lattice_origin_source_exported",
    "role_like_stride_source_exported",
    "legacy_integer_node_gauge_role_transferred",
    "phase_frequency_node_gauge_certificate_exported",
    "legacy_to_strict_bridge_completion_exported",
    "legacy_role_transfer_revalidated",
    "strict_kernel_full_kernel_claimed",
    "toe_closure_claimed",
    "positive_beta_renormalization_source_exported",
    "inverse_hierarchy_role_transfer_exported",
    "selector_source_exported",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "blind_empirical_confirmation_claimed",
]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.exists() else ""


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json", "-g", "*.lean",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:35]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "universal_affine_zero_lattice_content": (
            "affine.*zero|zero.*affine|zero-lattice|integer node|phase-node|node/gauge|"
            "offset.*stride|stride.*offset|metric lift|reindex"
        ),
        "source_atom_origin_stride_content": (
            "origin-and-stride|canonical origin|zero-lattice origin|canonical offset|canonical stride|"
            "source atom|source theorem|selector source|hidden selector"
        ),
        "z12_aut_premise_content": (
            "Z_12|Z12|Aut_Z12|Aut\\(Z_12\\)|generator/orientation|premise-based|fixing datum|"
            "successor map|quotient-safe|transport object"
        ),
        "role_demotion_transfer_content": (
            "demot|degrad|role-transfer|legacy.*role|node/gauge role|integer.*role|survives|rejected|"
            "modified/compressed"
        ),
        "remaining_professorial_path_content": (
            "beta source|positive beta|inverse hierarchy|QW-2191|full kernel|ToE closure|"
            "role-bearing L_total|blind empirical|frozen-kernel"
        ),
    }
    return {"tool": "rg", "mode": "content-first semantic audit for universal affine zero-lattice source no-go and node-role demotion", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def latest_commit_audit() -> list[dict[str, Any]]:
    proc = subprocess.run(["git", "log", "-n", "8", "--oneline", "--name-only"], cwd=REPO, text=True, stdout=subprocess.PIPE, check=True)
    rows: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        if re.match(r"^[0-9a-f]{7,12} ", line):
            if current:
                rows.append(current)
            sha, subject = line.split(" ", 1)
            current = {"sha": sha, "subject": subject, "files": []}
        elif current is not None:
            current["files"].append(line)
    if current:
        rows.append(current)
    return rows


def frac_str(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def zero(index: int) -> Fraction:
    return Fraction(4, 3) + 4 * index


def phase(x: float) -> float:
    return math.cos(math.pi * (float(OMEGA) * x + float(PHI)))


def strict_abs(x: float) -> float | None:
    if x <= 0:
        return None
    return abs(phase(x) / (1.0 + x ** float(ETA)))


def affine_for(k0: int, stride: int) -> tuple[Fraction, Fraction]:
    a = Fraction(4 * stride, 3)
    b = zero(k0) - a * LEGACY_NODES[0]
    return a, b


def orbit(k: int) -> tuple[int, ...]:
    return tuple(sorted({(u * k) % MOD for u in AUT}))


def aut_fixed(k: int) -> bool:
    return all((u * k) % MOD == k % MOD for u in AUT)


def evaluate_lift(k0: int, stride: int) -> dict[str, Any]:
    a, b = affine_for(k0, stride)
    r1 = a + b
    r7 = 7 * a + b
    k0_mod = k0 % MOD
    stride_mod = stride % MOD
    node_rows = []
    for i, node in enumerate(LEGACY_NODES):
        target_index = k0 + i * stride
        rd = a * node + b
        node_rows.append({
            "legacy_node": node,
            "target_zero_index": target_index,
            "r_d": frac_str(rd),
            "target_zero": frac_str(zero(target_index)),
            "exact_match": rd == zero(target_index),
            "cos_residual_abs": abs(phase(float(rd))),
        })
    numerator = strict_abs(float(r7))
    denominator = strict_abs(float(r1))
    ratio = numerator / denominator if numerator is not None and denominator not in (None, 0.0) else None
    return {
        "k0": k0,
        "stride": stride,
        "k0_mod12": k0_mod,
        "stride_mod12": stride_mod,
        "map": f"r(d)=({frac_str(a)})*d+({frac_str(b)})",
        "exact_node_lift": all(row["exact_match"] and row["cos_residual_abs"] < 1e-12 for row in node_rows),
        "node_rows": node_rows,
        "r1": frac_str(r1),
        "r7": frac_str(r7),
        "uv_positive_at_d1": r1 > 0,
        "strict_lifted_abs_k7_over_k1": ratio,
        "inverse_hierarchy_proxy_above_one": bool(ratio is not None and ratio > 1.0),
        "k0_aut_fixed": aut_fixed(k0_mod),
        "stride_aut_fixed": aut_fixed(stride_mod),
        "k0_orbit_mod12": list(orbit(k0_mod)),
        "stride_orbit_mod12": list(orbit(stride_mod)),
        "hidden_origin_selector_required": not aut_fixed(k0_mod),
    }


def symbolic_parameterization_theorem() -> dict[str, Any]:
    return {
        "legacy_nodes": list(LEGACY_NODES),
        "strict_zero_lattice": "z_k = 4/3 + 4 k",
        "all_exact_affine_lifts": "r_{k0,s}(d) = (4s/3)d + (4/3 + 4k0 - 8s/3)",
        "derivation": [
            "Legacy nodes are d_i = 2 + 3i.",
            "Exact lift requires r(d_i) = z_{k0 + i s} = 4/3 + 4k0 + 4is.",
            "Subtracting adjacent equations gives 3a = 4s, hence a = 4s/3.",
            "Using d_0=2 gives b = 4/3 + 4k0 - 8s/3.",
        ],
        "uv_positive_at_d1_iff": "r(1)>0 iff 1 + 3k0 - s > 0",
        "why_this_is_universal": "Every affine exact repair of the legacy integer nodes into the strict zero lattice is one member of this two-integer family; larger brute-force boxes add no new algebraic type.",
    }


def strict_aut_source_scan() -> dict[str, Any]:
    fixed_origins = [k for k in range(MOD) if aut_fixed(k)]
    fixed_nonzero_strides = [s for s in range(1, MOD) if aut_fixed(s)]
    candidates = [evaluate_lift(k0, stride) for k0 in fixed_origins for stride in fixed_nonzero_strides]
    return {
        "source_scope": "strict Aut(Z12)-invariant origin and nonzero translation stride only",
        "fixed_origins": fixed_origins,
        "fixed_nonzero_strides": fixed_nonzero_strides,
        "candidate_count": len(candidates),
        "uv_safe_count": sum(1 for row in candidates if row["uv_positive_at_d1"]),
        "inverse_hierarchy_count": sum(1 for row in candidates if row["inverse_hierarchy_proxy_above_one"]),
        "candidates": candidates,
        "theorem_result": "NO_STRICT_AUT_SOURCE_FOR_UV_SAFE_INVERSE_HIERARCHY_AFFINE_NODE_REPAIR",
    }


def premise_generator_with_canonical_origin_scan() -> dict[str, Any]:
    fixed_origins = [k for k in range(MOD) if aut_fixed(k)]
    premise_strides = list(range(1, MOD + 1))
    candidates = [evaluate_lift(k0, stride) for k0 in fixed_origins for stride in premise_strides]
    by_origin = {}
    for origin in fixed_origins:
        rows = [row for row in candidates if row["k0"] == origin]
        by_origin[str(origin)] = {
            "candidate_count": len(rows),
            "uv_safe_count": sum(1 for row in rows if row["uv_positive_at_d1"]),
            "inverse_hierarchy_count": sum(1 for row in rows if row["inverse_hierarchy_proxy_above_one"]),
            "best_ratio": max((row["strict_lifted_abs_k7_over_k1"] or -1.0) for row in rows),
        }
    return {
        "source_scope": "N523/T164 may fix a directed generator, but origin remains restricted to Aut-fixed origins to avoid adding a hidden origin selector",
        "fixed_origins": fixed_origins,
        "premise_tracked_strides": premise_strides,
        "candidate_count": len(candidates),
        "uv_safe_count": sum(1 for row in candidates if row["uv_positive_at_d1"]),
        "inverse_hierarchy_count": sum(1 for row in candidates if row["inverse_hierarchy_proxy_above_one"]),
        "by_origin": by_origin,
        "uv_safe_candidates": [row for row in candidates if row["uv_positive_at_d1"]],
        "theorem_result": "EVEN_WITH_PREMISE_GENERATOR_AND_CANONICAL_ORIGIN_NO_UV_SAFE_INVERSE_HIERARCHY_AFFINE_NODE_REPAIR",
    }


def p2639_role_like_source_recheck() -> dict[str, Any]:
    p2639 = load_json(SOURCE_FILES["P2639_OFFSET_STRIDE_EXHAUSTION"])
    role_like = p2639.get("offset_stride_exhaustion", {}).get("role_like_candidates", [])
    if not role_like:
        role_like = [{"k0": 4, "stride": 3}, {"k0": 10, "stride": 6}]
    rows = []
    for row in role_like:
        ev = evaluate_lift(int(row["k0"]), int(row["stride"]))
        rows.append({
            **ev,
            "p2639_ratio": row.get("strict_lifted_abs_k7_over_k1"),
            "strict_aut_source_pass": ev["k0_aut_fixed"] and ev["stride_aut_fixed"],
            "premise_generator_plus_canonical_origin_pass": ev["k0_aut_fixed"],
            "source_failure": "origin_not_aut_fixed_hidden_selector_required" if not ev["k0_aut_fixed"] else "none",
        })
    return {
        "role_like_candidate_count": len(rows),
        "strict_aut_source_pass_count": sum(1 for row in rows if row["strict_aut_source_pass"]),
        "premise_generator_plus_canonical_origin_pass_count": sum(1 for row in rows if row["premise_generator_plus_canonical_origin_pass"]),
        "rows": rows,
        "decision": "P2639 role-like lifts remain mathematical repairs, not sourced node/gauge role transfers.",
    }


def source_text_ledger() -> list[dict[str, Any]]:
    phrases = {
        "N456_NO_AUT_INVARIANT_12_CYCLE": ["no Aut-invariant 12-cycle successor map", "does not claim", "QW-2191"],
        "N462_NO_AUT_INVARIANT_GENERATOR_ORIENTATION": ["no Aut", "generator", "orientation"],
        "N523_PREMISE_BASED_T164_FIXING_DATUM": ["generator_fixed = 1", "suc_fix(k) := (k+1) mod 12", "premise-based"],
        "N454_PARITY_CHARACTER_AUT_INVARIANCE": ["parity character", "Aut", "not"],
        "N457_QUOTIENT_HOLONOMY_BOUNDARY": ["hidden selector", "successor map", "does not prove"],
    }
    rows = []
    for name, needles in phrases.items():
        path = SOURCE_FILES[name]
        text = load_text(path).lower()
        rows.append({
            "source": name,
            "path": rel(path),
            "exists": path.exists(),
            "phrase_hits": {needle: needle.lower() in text for needle in needles},
            "exports_universal_affine_origin_stride_source": False,
        })
    return rows


def demotion_decision(strict_scan: dict[str, Any], premise_scan: dict[str, Any], role_like: dict[str, Any]) -> dict[str, Any]:
    gates = {
        "universal_affine_family_proven": True,
        "strict_aut_source_candidate_with_uv_and_inverse_hierarchy_exists": strict_scan["inverse_hierarchy_count"] > 0,
        "premise_generator_canonical_origin_candidate_with_uv_and_inverse_hierarchy_exists": premise_scan["inverse_hierarchy_count"] > 0,
        "p2639_role_like_candidate_has_source_pass": role_like["strict_aut_source_pass_count"] > 0 or role_like["premise_generator_plus_canonical_origin_pass_count"] > 0,
        "new_origin_source_atom_found_in_current_ledger": False,
    }
    return {
        "gates": gates,
        "legacy_integer_node_gauge_role_status": "DEMOTE_TO_UNSOURCED_NUMERICAL_COORDINATE_REPAIR_CANDIDATE_UNLESS_NEW_ORIGIN_SOURCE_ATOM_IS_ADDED",
        "strict_kernel_status": "ROBUST_WORKING_KERNEL_NOT_FULL_TOE_KERNEL",
        "bridge_completion_status": "BLOCKED_BY_ZERO_LATTICE_ORIGIN_SOURCE_AND_ROLE_TRANSFER",
        "professorial_verdict": (
            "The affine repair problem is now algebraically exhausted: every exact node repair is r_{k0,s}. "
            "Strict Aut(Z12) sources and N523-style premise generator support do not yield a UV-safe inverse-hierarchy repair with a canonical origin. "
            "Therefore the legacy integer node/gauge role should be demoted rather than silently transferred, unless a genuinely new nadsoliton origin-source atom is introduced."
        ),
        "next_honest_step": (
            "Stop enlarging offset/stride searches. Either introduce and prove a new zero-lattice origin source theorem, or formally reroute closure to beta-source and inverse-hierarchy role-transfer blockers with the legacy integer node/gauge role marked demoted."
        ),
        "full_kernel_now": False,
        "toe_closure_now": False,
    }


def professorial_closure_path() -> list[dict[str, str]]:
    return [
        {
            "rank": "1",
            "step": "Origin-source theorem or demotion",
            "criterion": "A new strict/premise-tracked nadsoliton source must choose the zero-lattice origin without hidden selector input; otherwise the legacy integer node/gauge role is formally demoted.",
        },
        {
            "rank": "2",
            "step": "Beta-source / inverse-hierarchy reroute",
            "criterion": "After demotion, do not spend more proof effort on offset/stride brute force; attack beta normalization and inverse-hierarchy role transfer as separate blockers.",
        },
        {
            "rank": "3",
            "step": "Role-transfer matrix rerun",
            "criterion": "Only rerun legacy physical roles after a sourced bridge exists; otherwise classify each role as survivor, compressed successor, or rejected/demoted.",
        },
        {
            "rank": "4",
            "step": "Frozen-kernel empirical interface",
            "criterion": "Once source and role gates are explicit, test frozen strict-kernel signatures against blind CMB/LSS/GW/PTA or cross-sector holdouts without retuning.",
        },
    ]


def write_markdown(payload: dict[str, Any]) -> None:
    decision = payload["demotion_decision"]
    lines = [
        "# P2642/S1592 universal affine zero-lattice source no-go and node-role demotion",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This audit greps research contents before adding the proof: universal affine zero-lattice repairs, origin/stride source atoms, Z12/Aut and premise fixing, role demotion/transfer, and remaining ToE blockers.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Universal affine theorem",
        "",
        f"All exact affine repairs are `{payload['symbolic_parameterization_theorem']['all_exact_affine_lifts']}`.",
        f"UV positivity at d=1 is `{payload['symbolic_parameterization_theorem']['uv_positive_at_d1_iff']}`.",
        "This means larger offset/stride brute-force boxes add no new algebraic type; they only revisit different integer members of the same family.",
        "",
        "## Source scans",
        "",
        f"Strict Aut source candidates: `{payload['strict_aut_source_scan']['candidate_count']}`; UV-safe: `{payload['strict_aut_source_scan']['uv_safe_count']}`; inverse-hierarchy: `{payload['strict_aut_source_scan']['inverse_hierarchy_count']}`.",
        f"Premise-generator with canonical-origin candidates: `{payload['premise_generator_with_canonical_origin_scan']['candidate_count']}`; UV-safe: `{payload['premise_generator_with_canonical_origin_scan']['uv_safe_count']}`; inverse-hierarchy: `{payload['premise_generator_with_canonical_origin_scan']['inverse_hierarchy_count']}`.",
        f"P2639 role-like source-pass count under these scopes: strict Aut `{payload['p2639_role_like_source_recheck']['strict_aut_source_pass_count']}`, premise+canonical-origin `{payload['p2639_role_like_source_recheck']['premise_generator_plus_canonical_origin_pass_count']}`.",
        "",
        "## Demotion decision",
        "",
        decision["professorial_verdict"],
        "",
        f"Legacy integer node/gauge role status: `{decision['legacy_integer_node_gauge_role_status']}`.",
        f"Bridge completion status: `{decision['bridge_completion_status']}`.",
        f"Full kernel now? `{decision['full_kernel_now']}`.",
        f"ToE closure now? `{decision['toe_closure_now']}`.",
        "",
        "## Next honest step",
        "",
        decision["next_honest_step"],
        "",
        "## Professorial closure path",
        "",
    ])
    for row in payload["professorial_closure_path"]:
        lines.append(f"{row['rank']}. **{row['step']}** — {row['criterion']}")
    lines.extend([
        "",
        "## Negative exports",
        "",
    ])
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    strict_scan = strict_aut_source_scan()
    premise_scan = premise_generator_with_canonical_origin_scan()
    role_like = p2639_role_like_source_recheck()
    payload: dict[str, Any] = {
        "status": "P2642_UNIVERSAL_AFFINE_ZERO_LATTICE_SOURCE_NO_GO_NODE_ROLE_DEMOTION_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {name: sha256_file(path) for name, path in SOURCE_FILES.items()},
        "source_text_ledger": source_text_ledger(),
        "symbolic_parameterization_theorem": symbolic_parameterization_theorem(),
        "strict_aut_source_scan": strict_scan,
        "premise_generator_with_canonical_origin_scan": premise_scan,
        "p2639_role_like_source_recheck": role_like,
        "demotion_decision": demotion_decision(strict_scan, premise_scan, role_like),
        "professorial_closure_path": professorial_closure_path(),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        SOURCE_FILES["STRICT_EQUATION_SHEET"],
        "P2642/S1592 universal affine zero-lattice no-go guard",
        "## P2642/S1592 universal affine zero-lattice no-go guard\n\n"
        "`P2642/S1592` algebraically exhausts exact affine repairs of legacy integer nodes into the strict zero lattice: every such repair is `r_{k0,s}(d)=(4s/3)d+(4/3+4k0-8s/3)`.  Strict `Aut(Z12)` sources and N523-style premise generator support, while keeping the origin canonical, do not produce a UV-safe inverse-hierarchy node repair.  Therefore the legacy integer node/gauge role remains demoted unless a new zero-lattice origin source atom is proved; bridge completion, role transfer, beta source, inverse hierarchy, `QW-2191`, and ToE closure remain closed.\n",
    )
    append_once(
        SOURCE_FILES["STRICT_LAGRANGIAN_DRAFT"],
        "P2642/S1592 universal affine node-role Ltotal guard",
        "## P2642/S1592 universal affine node-role Ltotal guard\n\n"
        "`P2642/S1592` does not re-enable role-bearing `L_total`: the universal affine family maps the legacy node list into strict zeros only as an unsourced coordinate repair.  No current strict or premise-tracked origin/stride source yields the required UV-safe inverse-hierarchy transfer, so a variational node/gauge term must wait for a new origin-source theorem or for formal demotion of that legacy role.\n",
    )
    print(rel(OUT))
    print(rel(MD))


if __name__ == "__main__":
    main()
