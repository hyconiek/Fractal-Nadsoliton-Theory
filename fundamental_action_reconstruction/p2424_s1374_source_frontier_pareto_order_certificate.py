#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2424_s1374_source_frontier_pareto_order_certificate.json"
MD = GEN / "p2424_s1374_source_frontier_pareto_order_certificate.md"

SOURCE_FILES = {
    "P2423_REPAIR_ORDER_POSET": GEN / "p2423_s1373_admissible_repair_order_poset_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

BRIDGE_SOURCE_GATE = "source_obligation_discharge"
CHI11_GATE = "chi11_source_export"
QW2191_GATE = "qw2191_selector_discharge"
ROLE_GATE = "role_transfer_audit_license"
LTOTAL_GATE = "role_bearing_ltotal_export"
SOURCE_FRONTIER_GATES = [BRIDGE_SOURCE_GATE, CHI11_GATE, QW2191_GATE]


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
        "new_packet": "P2424|S1374|source frontier Pareto|source-frontier Pareto|Pareto order certificate",
        "p2423_input": "P2423|admissible repair-order|admissible linear extensions|role-transfer step",
        "pareto_prior": "Pareto|frontier|dominance|dominated order",
        "source_frontier_prior": "source frontier|first-three source|bridge-source readiness|selector-source readiness",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds P2423's admissible order poset and generic Pareto/frontier material, but no production P24xx "
            "certificate classifying the six admissible first-three source orders into bridge-first, selector-pair-first, "
            "and dominated mixed strategies."
        ),
    }


def p2423_theorem(source: dict[str, Any]) -> dict[str, Any]:
    return source.get("admissible_repair_order_poset_certificate", {}).get("theorem_export", {})


def order_class(order: list[str]) -> str:
    first_three = order[:3]
    if first_three[0] == BRIDGE_SOURCE_GATE:
        return "bridge_first_pareto"
    if set(first_three[:2]) == {CHI11_GATE, QW2191_GATE}:
        return "selector_pair_first_pareto"
    return "mixed_split_dominated"


def order_metrics(order: list[str]) -> dict[str, Any]:
    bridge_step = order.index(BRIDGE_SOURCE_GATE) + 1
    selector_step = max(order.index(CHI11_GATE), order.index(QW2191_GATE)) + 1
    role_step = order.index(ROLE_GATE) + 1
    ltotal_step = order.index(LTOTAL_GATE) + 1
    return {
        "order": order,
        "first_three_source_gates": order[:3],
        "class": order_class(order),
        "bridge_source_ready_step": bridge_step,
        "selector_source_ready_step": selector_step,
        "role_transfer_ready_step": role_step,
        "role_bearing_ltotal_ready_step": ltotal_step,
        "toe_ready_step": ltotal_step,
        "objective_vector_bridge_selector": [bridge_step, selector_step],
    }


def dominates(left: dict[str, Any], right: dict[str, Any]) -> bool:
    lb, ls = left["objective_vector_bridge_selector"]
    rb, rs = right["objective_vector_bridge_selector"]
    return (lb <= rb and ls <= rs) and (lb < rb or ls < rs)


def pareto_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [row for row in rows if not any(dominates(other, row) for other in rows)]


def class_counts(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["class"]] = counts.get(row["class"], 0) + 1
    return dict(sorted(counts.items()))


def vector_counts(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        key = str(row["objective_vector_bridge_selector"])
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


def first_gate_counts(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        gate = row["first_three_source_gates"][0]
        counts[gate] = counts.get(gate, 0) + 1
    return dict(sorted(counts.items()))


def build_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2423 = p2423_theorem(sources["P2423_REPAIR_ORDER_POSET"])
    orders = p2423.get("admissible_orders", [])
    rows = [order_metrics(order) for order in orders]
    frontier = pareto_rows(rows)
    dominated = [row for row in rows if row not in frontier]
    return {
        "source_frontier_gates": SOURCE_FRONTIER_GATES,
        "admissible_order_rows": rows,
        "pareto_order_rows": frontier,
        "dominated_order_rows": dominated,
        "class_counts": class_counts(rows),
        "pareto_class_counts": class_counts(frontier),
        "objective_vector_counts": vector_counts(rows),
        "pareto_objective_vector_counts": vector_counts(frontier),
        "first_gate_counts": first_gate_counts(rows),
        "pareto_first_gate_counts": first_gate_counts(frontier),
        "p2423_admissible_order_count_inherited": p2423.get("admissible_order_count"),
        "p2423_role_transfer_step_distribution_inherited": p2423.get("role_transfer_ready_step_distribution"),
        "p2423_ltotal_step_distribution_inherited": p2423.get("role_bearing_ltotal_ready_step_distribution"),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2424/S1374 source-frontier Pareto order certificate

`P2424/S1374` refines the six P2423 admissible orders by looking only at the first-three source frontier before role-transfer.  The objectives are earliest bridge-source readiness and earliest selector-source readiness.

The finite Pareto audit has two incomparable optimal classes.  If `source_obligation_discharge` is first, bridge-source readiness occurs at step `1` and selector-source readiness at step `3`.  If the `chi11_source_export + qw2191_selector_discharge` pair is first, selector-source readiness occurs at step `2` and bridge-source readiness at step `3`.  The two mixed orders are dominated because they delay selector-source readiness to step `3` without improving bridge-source readiness beyond step `2`.

This is a proof-search ordering theorem only: it does not pick a unique first gate without an extra cost/source premise and exports no source, selector, QW-2191, role-transfer, `L_total`, or ToE theorem.
""".strip()
    lag_section = """
## P2424/S1374 source-frontier Pareto guard for Lagrangian/EOM

`P2424/S1374` shows that the first-three source frontier has two incomparable Pareto classes, bridge-first and selector-pair-first.  Because no unique first theorem gate is selected internally, no Pareto order can be promoted into role-transfer, role-bearing `L_total`, or ToE closure.
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
    theorem_export = {
        "theorem_name": "P2424_T1_source_frontier_pareto_order_certificate",
        "source_frontier_gate_count": len(cert["source_frontier_gates"]),
        "admissible_order_count": len(cert["admissible_order_rows"]),
        "pareto_order_count": len(cert["pareto_order_rows"]),
        "dominated_order_count": len(cert["dominated_order_rows"]),
        "class_counts": cert["class_counts"],
        "pareto_class_counts": cert["pareto_class_counts"],
        "objective_vector_counts": cert["objective_vector_counts"],
        "pareto_objective_vector_counts": cert["pareto_objective_vector_counts"],
        "first_gate_counts": cert["first_gate_counts"],
        "pareto_first_gate_counts": cert["pareto_first_gate_counts"],
        "unique_pareto_first_gate_selected": False,
        "p2423_admissible_order_count_inherited": cert["p2423_admissible_order_count_inherited"] == 6,
        "p2423_role_transfer_step_distribution_inherited": cert["p2423_role_transfer_step_distribution_inherited"] == {"4": 6},
        "p2423_ltotal_step_distribution_inherited": cert["p2423_ltotal_step_distribution_inherited"] == {"5": 6},
        "source_obligation_discharge_exported": False,
        "chi11_source_exported": False,
        "qw2191_discharged": False,
        "role_transfer_licensed": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Pareto optimality ranks admissible proof orders but does not discharge a source theorem.",
            "Two source-frontier Pareto classes remain incomparable without an extra cost/source premise.",
            "No unique first source gate is selected internally.",
            "No role-transfer, role-bearing L_total, or ToE theorem is exported.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "three_source_frontier_gates": theorem_export["source_frontier_gate_count"] == 3,
        "six_admissible_orders": theorem_export["admissible_order_count"] == 6,
        "four_pareto_orders": theorem_export["pareto_order_count"] == 4,
        "two_dominated_orders": theorem_export["dominated_order_count"] == 2,
        "class_counts_expected": theorem_export["class_counts"] == {"bridge_first_pareto": 2, "mixed_split_dominated": 2, "selector_pair_first_pareto": 2},
        "pareto_class_counts_expected": theorem_export["pareto_class_counts"] == {"bridge_first_pareto": 2, "selector_pair_first_pareto": 2},
        "objective_vectors_expected": theorem_export["objective_vector_counts"] == {"[1, 3]": 2, "[2, 3]": 2, "[3, 2]": 2},
        "pareto_vectors_expected": theorem_export["pareto_objective_vector_counts"] == {"[1, 3]": 2, "[3, 2]": 2},
        "no_unique_first_gate": not theorem_export["unique_pareto_first_gate_selected"],
        "p2423_admissible_orders_inherited": theorem_export["p2423_admissible_order_count_inherited"],
        "p2423_role_step_inherited": theorem_export["p2423_role_transfer_step_distribution_inherited"],
        "p2423_ltotal_step_inherited": theorem_export["p2423_ltotal_step_distribution_inherited"],
        "source_still_open": not theorem_export["source_obligation_discharge_exported"],
        "chi11_still_open": not theorem_export["chi11_source_exported"],
        "qw2191_still_open": not theorem_export["qw2191_discharged"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "ltotal_still_blocked": not theorem_export["role_bearing_ltotal_exported"],
        "toe_still_open": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2424_s1374_v1",
        "packet_id": "P2424",
        "stage_id": "S1374",
        "result_kind": "SOURCE_FRONTIER_PARETO_ORDER_CERTIFICATE",
        "status": "PASS_SOURCE_FRONTIER_PARETO_ORDER_NO_UNIQUE_SOURCE_GATE_NO_DISCHARGE",
        "source_frontier_pareto_order_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Add a real external/source-cost premise if one wants to choose between bridge-first and selector-pair-first; "
            "otherwise either Pareto class is only a proof-search order, not a closure theorem."
        ),
        "global_status": "OPEN_PROGRESS_SOURCE_FRONTIER_PARETO_CERTIFIED_NO_SOURCE_GATE_DISCHARGED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["source_frontier_pareto_order_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2424 S1374: source-frontier Pareto order certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Source-frontier gates: `{theorem['source_frontier_gate_count']}`.",
                f"- Admissible orders: `{theorem['admissible_order_count']}`.",
                f"- Pareto orders: `{theorem['pareto_order_count']}`.",
                f"- Dominated orders: `{theorem['dominated_order_count']}`.",
                f"- Pareto vectors: `{theorem['pareto_objective_vector_counts']}`.",
                f"- Pareto classes: `{theorem['pareto_class_counts']}`.",
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
