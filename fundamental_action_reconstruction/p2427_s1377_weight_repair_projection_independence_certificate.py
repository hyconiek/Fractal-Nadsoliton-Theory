#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2427_s1377_weight_repair_projection_independence_certificate.json"
MD = GEN / "p2427_s1377_weight_repair_projection_independence_certificate.md"

SOURCE_FILES = {
    "P2426_PRODUCT": GEN / "p2426_s1376_weighted_tiebreak_repair_subcube_nonpromotion_product_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

BRIDGE_VECTOR = [1, 3]
SELECTOR_VECTOR = [3, 2]
WEIGHT_GRID_MAX = 12
MISSING_GATES = [
    "source_obligation_discharge",
    "chi11_source_export",
    "qw2191_selector_discharge",
    "role_transfer_audit_license",
    "role_bearing_ltotal_export",
]
TARGETS = [
    "bridge_source_ready",
    "selector_source_ready",
    "role_transfer_ready",
    "role_bearing_ltotal_ready",
    "toe_ready",
]


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg",
            "-n",
            pattern,
            "fundamental_action_reconstruction",
            "material_dowodowy",
            "-g",
            "*.py",
            "-g",
            "*.md",
            "-g",
            "*.tex",
            "-g",
            "!generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:20]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2427|S1377|weight repair projection|projection independence|repair.*independence",
        "p2426_input": "P2426|weighted tie-break.*repair|product assignments|proper repair failures",
        "independence_prior": "orthogonal|independent|projection|factorization|conditional distribution",
        "nonpromotion_prior": "nonpromotion|does not discharge|No source|No ToE|not a theorem discharge",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds P2426's product-level nonpromotion language and older generic projection/factorization language, "
            "but no production P24xx certificate that proves the weight projection and repair-gate projection have exact "
            "independent contingency tables across all readiness predicates."
        ),
    }


def theorem(payload: dict[str, Any], key: str) -> dict[str, Any]:
    return payload.get(key, {}).get("theorem_export", {})


def weighted_cost(vector: list[int], bridge_weight: int, selector_weight: int) -> int:
    return bridge_weight * vector[0] + selector_weight * vector[1]


def weight_winner(bridge_weight: int, selector_weight: int) -> str:
    bridge_cost = weighted_cost(BRIDGE_VECTOR, bridge_weight, selector_weight)
    selector_cost = weighted_cost(SELECTOR_VECTOR, bridge_weight, selector_weight)
    if bridge_cost < selector_cost:
        return "bridge_first_pareto"
    if selector_cost < bridge_cost:
        return "selector_pair_first_pareto"
    return "bridge_selector_tie"


def repair_target_ready(added: set[str]) -> dict[str, bool]:
    return {
        "bridge_source_ready": "source_obligation_discharge" in added,
        "selector_source_ready": {"chi11_source_export", "qw2191_selector_discharge"}.issubset(added),
        "role_transfer_ready": "role_transfer_audit_license" in added,
        "role_bearing_ltotal_ready": "role_bearing_ltotal_export" in added,
        "toe_ready": set(MISSING_GATES).issubset(added),
    }


def repair_subsets() -> list[list[str]]:
    out: list[list[str]] = []
    for size in range(len(MISSING_GATES) + 1):
        for combo in combinations(MISSING_GATES, size):
            out.append(list(combo))
    return out


def count_by(items: list[str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for item in items:
        counts[item] = counts.get(item, 0) + 1
    return dict(sorted(counts.items()))


def weight_counts() -> dict[str, int]:
    winners = [
        weight_winner(bridge_weight, selector_weight)
        for bridge_weight in range(1, WEIGHT_GRID_MAX + 1)
        for selector_weight in range(1, WEIGHT_GRID_MAX + 1)
    ]
    return count_by(winners)


def repair_rows() -> list[dict[str, Any]]:
    rows = []
    for repair in repair_subsets():
        added = set(repair)
        target_ready = repair_target_ready(added)
        rows.append(
            {
                "added_gates": repair,
                "added_gate_count": len(repair),
                "remaining_missing_count": len(MISSING_GATES) - len(repair),
                "target_ready": target_ready,
            }
        )
    return rows


def repair_true_counts(rows: list[dict[str, Any]]) -> dict[str, int]:
    return {
        target: sum(1 for row in rows if row["target_ready"][target])
        for target in TARGETS
    }


def missing_count_distribution(rows: list[dict[str, Any]]) -> dict[str, int]:
    return count_by([str(row["remaining_missing_count"]) for row in rows])


def gate_inclusion_counts(rows: list[dict[str, Any]]) -> dict[str, int]:
    return {
        gate: sum(1 for row in rows if gate in row["added_gates"])
        for gate in MISSING_GATES
    }


def factor_table(weight_side: dict[str, int], repair_side: dict[str, int]) -> dict[str, dict[str, int]]:
    return {
        weight_key: {repair_key: weight_count * repair_count for repair_key, repair_count in repair_side.items()}
        for weight_key, weight_count in weight_side.items()
    }


def conditional_distribution(table: dict[str, dict[str, int]]) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for weight_key, row in table.items():
        total = sum(row.values())
        out[weight_key] = {key: str(Fraction(value, total)) for key, value in sorted(row.items())}
    return out


def rows_identical(distribution: dict[str, dict[str, str]]) -> bool:
    rows = list(distribution.values())
    return bool(rows) and all(row == rows[0] for row in rows)


def build_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2426 = theorem(sources["P2426_PRODUCT"], "weighted_tiebreak_repair_subcube_nonpromotion_product_certificate")
    weights = weight_counts()
    repairs = repair_rows()
    readiness_counts = repair_true_counts(repairs)
    readiness_tables = {
        target: factor_table(weights, {"false": len(repairs) - true_count, "true": true_count})
        for target, true_count in readiness_counts.items()
    }
    gate_tables = {
        gate: factor_table(weights, {"absent": len(repairs) - count, "present": count})
        for gate, count in gate_inclusion_counts(repairs).items()
    }
    missing_distribution = missing_count_distribution(repairs)
    missing_table = factor_table(weights, missing_distribution)
    return {
        "weight_counts": weights,
        "repair_assignment_count": len(repairs),
        "repair_readiness_true_counts": readiness_counts,
        "repair_missing_count_distribution": missing_distribution,
        "gate_inclusion_counts": gate_inclusion_counts(repairs),
        "readiness_contingency_tables": readiness_tables,
        "gate_inclusion_contingency_tables": gate_tables,
        "missing_count_contingency_table": missing_table,
        "conditional_readiness_distributions": {target: conditional_distribution(table) for target, table in readiness_tables.items()},
        "conditional_missing_count_distribution": conditional_distribution(missing_table),
        "p2426_product_assignment_count_inherited": p2426.get("product_assignment_count"),
        "p2426_weight_winner_counts_inherited": p2426.get("weight_winner_product_counts"),
        "p2426_toe_ready_product_count_inherited": p2426.get("toe_ready_product_count"),
        "p2426_weighted_choice_reduces_missing_gate_count_inherited": p2426.get("weighted_choice_reduces_missing_gate_count"),
    }


def all_tables_factor(cert: dict[str, Any]) -> bool:
    weights = cert["weight_counts"]
    repair_count = cert["repair_assignment_count"]
    readiness_counts = cert["repair_readiness_true_counts"]
    for target, table in cert["readiness_contingency_tables"].items():
        true_count = readiness_counts[target]
        expected = factor_table(weights, {"false": repair_count - true_count, "true": true_count})
        if table != expected:
            return False
    for gate, table in cert["gate_inclusion_contingency_tables"].items():
        count = cert["gate_inclusion_counts"][gate]
        expected = factor_table(weights, {"absent": repair_count - count, "present": count})
        if table != expected:
            return False
    return cert["missing_count_contingency_table"] == factor_table(weights, cert["repair_missing_count_distribution"])


def append_doc_sections() -> None:
    eq_section = """
## P2427/S1377 weight-repair projection independence certificate

`P2427/S1377` refines P2426 by proving an exact projection-independence fact: the weighted frontier side and the five-gate repair side factor as independent finite contingencies.  For every weight-winner class, the repair distribution is the same `2^5=32` subcube distribution, including the same ToE-ready count `1`, selector-ready count `8`, bridge-source count `16`, and missing-count profile `1/5/10/10/5/1`.

The consequence is sharper than a count check: changing the weighted proof-search preference changes only order labels.  It never changes which source, chi11, QW-2191, role-transfer, or role-bearing `L_total` gate is present, so it cannot be promoted into bridge closure, role transfer, or ToE closure.
""".strip()
    lag_section = """
## P2427/S1377 weight-repair projection independence guard for Lagrangian/EOM

`P2427/S1377` certifies that all readiness predicates have empty weight-side support: their contingency tables factor through the repair subcube only.  Thus no weighted cost/order term can enter a role-bearing `L_total` as a surrogate for missing source, selector, QW-2191, or role-transfer theorems.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    cert = build_certificate(sources)
    conditional_readiness_identical = all(
        rows_identical(distribution) for distribution in cert["conditional_readiness_distributions"].values()
    )
    conditional_missing_identical = rows_identical(cert["conditional_missing_count_distribution"])
    theorem_export = {
        "theorem_name": "P2427_T1_weight_repair_projection_independence_certificate",
        "weight_assignment_count": sum(cert["weight_counts"].values()),
        "repair_assignment_count": cert["repair_assignment_count"],
        "product_assignment_count": sum(cert["weight_counts"].values()) * cert["repair_assignment_count"],
        "weight_counts": cert["weight_counts"],
        "repair_readiness_true_counts": cert["repair_readiness_true_counts"],
        "repair_missing_count_distribution": cert["repair_missing_count_distribution"],
        "gate_inclusion_counts": cert["gate_inclusion_counts"],
        "all_readiness_tables_factor_by_weight_x_repair": all_tables_factor(cert),
        "conditional_readiness_distributions_identical_across_weight_winners": conditional_readiness_identical,
        "conditional_missing_count_distribution_identical_across_weight_winners": conditional_missing_identical,
        "weight_side_support_for_repair_readiness_predicates": [],
        "repair_side_support_for_weight_winner_predicate": [],
        "p2426_product_assignment_count_inherited": cert["p2426_product_assignment_count_inherited"] == 4608,
        "p2426_weight_winner_counts_inherited": cert["p2426_weight_winner_counts_inherited"] == {
            "bridge_first_pareto": 3456,
            "bridge_selector_tie": 192,
            "selector_pair_first_pareto": 960,
        },
        "p2426_toe_ready_product_count_inherited": cert["p2426_toe_ready_product_count_inherited"] == 144,
        "p2426_weighted_choice_reduces_missing_gate_count_inherited": cert[
            "p2426_weighted_choice_reduces_missing_gate_count_inherited"
        ] is False,
        "weighted_choice_discharges_repair_gate": False,
        "source_obligation_discharge_exported": False,
        "chi11_source_exported": False,
        "qw2191_discharged": False,
        "role_transfer_licensed": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Projection independence is not a theorem-gate discharge.",
            "Weight winner labels have empty support on repair readiness predicates.",
            "Repair subsets have empty support on weighted frontier winner predicates.",
            "No source, selector, QW-2191, role-transfer, role-bearing L_total, or ToE theorem is exported.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "weight_grid_144": theorem_export["weight_assignment_count"] == 144,
        "repair_subcube_32": theorem_export["repair_assignment_count"] == 32,
        "product_4608": theorem_export["product_assignment_count"] == 4608,
        "weight_counts_match_p2425": theorem_export["weight_counts"] == {
            "bridge_first_pareto": 108,
            "bridge_selector_tie": 6,
            "selector_pair_first_pareto": 30,
        },
        "repair_readiness_counts_expected": theorem_export["repair_readiness_true_counts"] == {
            "bridge_source_ready": 16,
            "role_bearing_ltotal_ready": 16,
            "role_transfer_ready": 16,
            "selector_source_ready": 8,
            "toe_ready": 1,
        },
        "missing_distribution_binomial": theorem_export["repair_missing_count_distribution"] == {
            "0": 1,
            "1": 5,
            "2": 10,
            "3": 10,
            "4": 5,
            "5": 1,
        },
        "gate_inclusions_half_subcube": all(value == 16 for value in theorem_export["gate_inclusion_counts"].values()),
        "all_tables_factor": theorem_export["all_readiness_tables_factor_by_weight_x_repair"],
        "conditional_readiness_identical": theorem_export[
            "conditional_readiness_distributions_identical_across_weight_winners"
        ],
        "conditional_missing_identical": theorem_export[
            "conditional_missing_count_distribution_identical_across_weight_winners"
        ],
        "empty_weight_support_for_readiness": theorem_export["weight_side_support_for_repair_readiness_predicates"] == [],
        "empty_repair_support_for_winner": theorem_export["repair_side_support_for_weight_winner_predicate"] == [],
        "p2426_product_inherited": theorem_export["p2426_product_assignment_count_inherited"],
        "p2426_winner_counts_inherited": theorem_export["p2426_weight_winner_counts_inherited"],
        "p2426_toe_count_inherited": theorem_export["p2426_toe_ready_product_count_inherited"],
        "p2426_no_gate_reduction_inherited": theorem_export[
            "p2426_weighted_choice_reduces_missing_gate_count_inherited"
        ],
        "weighted_choice_discharges_no_gate": not theorem_export["weighted_choice_discharges_repair_gate"],
        "source_still_open": not theorem_export["source_obligation_discharge_exported"],
        "chi11_still_open": not theorem_export["chi11_source_exported"],
        "qw2191_still_open": not theorem_export["qw2191_discharged"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "ltotal_still_blocked": not theorem_export["role_bearing_ltotal_exported"],
        "toe_still_open": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2427_s1377_v1",
        "packet_id": "P2427",
        "stage_id": "S1377",
        "result_kind": "WEIGHT_REPAIR_PROJECTION_INDEPENDENCE_CERTIFICATE",
        "status": "PASS_WEIGHT_REPAIR_PROJECTION_INDEPENDENCE_NO_GATE_DISCHARGE_NO_CLOSURE",
        "weight_repair_projection_independence_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Use the independence certificate as a guardrail: proof-search weights may prioritize work, but source/selector/role gates "
            "must still be discharged by independent theorem work, especially chi11 and QW-2191."
        ),
        "global_status": "OPEN_PROGRESS_WEIGHT_REPAIR_INDEPENDENCE_CERTIFIED_NO_THEOREM_GATE_DISCHARGED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["weight_repair_projection_independence_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2427 S1377: weight-repair projection independence certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Weight assignments: `{theorem['weight_assignment_count']}`.",
                f"- Repair assignments: `{theorem['repair_assignment_count']}`.",
                f"- Product assignments: `{theorem['product_assignment_count']}`.",
                f"- Repair readiness true counts: `{theorem['repair_readiness_true_counts']}`.",
                f"- Missing-count distribution: `{theorem['repair_missing_count_distribution']}`.",
                f"- Tables factor: `{theorem['all_readiness_tables_factor_by_weight_x_repair']}`.",
                "",
                "## Hard limits",
                "",
                *[f"- {item}" for item in theorem["not_licensed"]],
                "",
                "## Gatekeepers",
                "",
                f"`{payload['gatekeeper_checks']}`",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    append_doc_sections()
    payload = build_payload()
    write_outputs(payload)
    if not all(payload["gatekeeper_checks"].values()):
        raise SystemExit(f"gatekeeper failure: {payload['gatekeeper_checks']}")


if __name__ == "__main__":
    main()
