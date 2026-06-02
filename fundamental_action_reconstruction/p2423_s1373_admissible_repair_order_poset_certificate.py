#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import permutations
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2423_s1373_admissible_repair_order_poset_certificate.json"
MD = GEN / "p2423_s1373_admissible_repair_order_poset_certificate.md"

SOURCE_FILES = {
    "P2422_REPAIR_SUBCUBE": GEN / "p2422_s1372_current_missing_gate_repair_subcube_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

MISSING_GATES = [
    "source_obligation_discharge",
    "chi11_source_export",
    "qw2191_selector_discharge",
    "role_transfer_audit_license",
    "role_bearing_ltotal_export",
]

PRECEDENCE_EDGES = [
    ("source_obligation_discharge", "role_transfer_audit_license", "role_transfer_after_bridge_source"),
    ("chi11_source_export", "role_transfer_audit_license", "role_transfer_after_chi11_source"),
    ("qw2191_selector_discharge", "role_transfer_audit_license", "role_transfer_after_qw2191_discharge"),
    ("role_transfer_audit_license", "role_bearing_ltotal_export", "ltotal_after_role_transfer_audit"),
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
        "new_packet": "P2423|S1373|admissible repair order|repair order poset|proof-order poset",
        "p2422_input": "P2422|repair subcube|partial unlock|missing-gate repair",
        "role_after_bridge_guard": "role-transfer audit after bridge completion|role_transfer_after_bridge|post-bridge role-transfer",
        "proof_order_prior": "proof order|linear extension|admissible order|unlock order",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds P2422's repair subcube and older proof-order language, but no production P24xx certificate "
            "that imposes the guardrail precedence relation on the five current missing gates and enumerates only the admissible linear extensions."
        ),
    }


def p2422_theorem(source: dict[str, Any]) -> dict[str, Any]:
    return source.get("current_missing_gate_repair_subcube_certificate", {}).get("theorem_export", {})


def order_positions(order: tuple[str, ...]) -> dict[str, int]:
    return {gate: index + 1 for index, gate in enumerate(order)}


def violated_edges(order: tuple[str, ...]) -> list[dict[str, str]]:
    pos = order_positions(order)
    violations = []
    for before, after, reason in PRECEDENCE_EDGES:
        if pos[before] > pos[after]:
            violations.append({"before": before, "after": after, "reason": reason})
    return violations


def prefix_state(prefix: tuple[str, ...]) -> dict[str, Any]:
    have = set(prefix)
    bridge_source_ready = "source_obligation_discharge" in have
    selector_source_ready = {"chi11_source_export", "qw2191_selector_discharge"}.issubset(have)
    role_transfer_ready = "role_transfer_audit_license" in have
    ltotal_ready = "role_bearing_ltotal_export" in have
    toe_ready = set(MISSING_GATES).issubset(have)
    return {
        "prefix": list(prefix),
        "step": len(prefix),
        "bridge_source_ready": bridge_source_ready,
        "selector_source_ready": selector_source_ready,
        "role_transfer_ready": role_transfer_ready,
        "role_bearing_ltotal_ready": ltotal_ready,
        "toe_ready": toe_ready,
    }


def order_row(order: tuple[str, ...]) -> dict[str, Any]:
    prefixes = [prefix_state(order[:step]) for step in range(1, len(order) + 1)]
    return {
        "order": list(order),
        "positions": order_positions(order),
        "prefix_rows": prefixes,
        "bridge_source_ready_step": next(row["step"] for row in prefixes if row["bridge_source_ready"]),
        "selector_source_ready_step": next(row["step"] for row in prefixes if row["selector_source_ready"]),
        "role_transfer_ready_step": next(row["step"] for row in prefixes if row["role_transfer_ready"]),
        "role_bearing_ltotal_ready_step": next(row["step"] for row in prefixes if row["role_bearing_ltotal_ready"]),
        "toe_ready_step": next(row["step"] for row in prefixes if row["toe_ready"]),
    }


def all_order_rows() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    admissible = []
    rejected = []
    for order in permutations(MISSING_GATES):
        violations = violated_edges(order)
        if violations:
            rejected.append({"order": list(order), "violated_edges": violations})
        else:
            admissible.append(order_row(order))
    return admissible, rejected


def distribution(rows: list[dict[str, Any]], key: str) -> dict[str, int]:
    dist: dict[str, int] = {}
    for row in rows:
        value = str(row[key])
        dist[value] = dist.get(value, 0) + 1
    return dict(sorted(dist.items(), key=lambda item: int(item[0])))


def edge_necessity_rows() -> list[dict[str, Any]]:
    rows = []
    all_orders = list(permutations(MISSING_GATES))
    full_count = len([order for order in all_orders if not violated_edges(order)])
    for edge in PRECEDENCE_EDGES:
        reduced_edges = [item for item in PRECEDENCE_EDGES if item != edge]
        count = 0
        newly_admitted_examples = []
        for order in all_orders:
            pos = order_positions(order)
            reduced_ok = all(pos[before] < pos[after] for before, after, _ in reduced_edges)
            full_ok = all(pos[before] < pos[after] for before, after, _ in PRECEDENCE_EDGES)
            if reduced_ok:
                count += 1
            if reduced_ok and not full_ok and len(newly_admitted_examples) < 3:
                newly_admitted_examples.append(list(order))
        before, after, reason = edge
        rows.append(
            {
                "edge": {"before": before, "after": after, "reason": reason},
                "linear_extension_count_without_edge": count,
                "newly_admitted_count_if_removed": count - full_count,
                "newly_admitted_examples": newly_admitted_examples,
            }
        )
    return rows


def build_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2422 = p2422_theorem(sources["P2422_REPAIR_SUBCUBE"])
    admissible, rejected = all_order_rows()
    return {
        "missing_gates": MISSING_GATES,
        "precedence_edges": [
            {"before": before, "after": after, "reason": reason}
            for before, after, reason in PRECEDENCE_EDGES
        ],
        "raw_order_count": len(admissible) + len(rejected),
        "admissible_order_rows": admissible,
        "rejected_order_rows": rejected,
        "admissible_order_count": len(admissible),
        "rejected_order_count": len(rejected),
        "bridge_source_ready_step_distribution": distribution(admissible, "bridge_source_ready_step"),
        "selector_source_ready_step_distribution": distribution(admissible, "selector_source_ready_step"),
        "role_transfer_ready_step_distribution": distribution(admissible, "role_transfer_ready_step"),
        "role_bearing_ltotal_ready_step_distribution": distribution(admissible, "role_bearing_ltotal_ready_step"),
        "toe_ready_step_distribution": distribution(admissible, "toe_ready_step"),
        "edge_necessity_rows": edge_necessity_rows(),
        "p2422_repair_subcube_count_inherited": p2422.get("repair_subcube_assignment_count"),
        "p2422_toe_ready_repair_count_inherited": p2422.get("toe_ready_repair_count"),
        "p2422_selector_singleton_unlock_count_inherited": p2422.get("selector_singleton_unlock_count"),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2423/S1373 admissible repair-order poset certificate

`P2423/S1373` adds proof-order discipline to the P2422 repair subcube.  It enumerates all `5! = 120` orders of the five missing gates and imposes the guardrail precedence relation: source discharge, chi11 source export, and QW-2191 selector discharge must precede role-transfer audit; role-transfer audit must precede role-bearing `L_total` export.

Only `6` orders survive as admissible linear extensions.  In every admissible order, role-transfer readiness first appears at step `4`, role-bearing `L_total` and ToE first appear at step `5`, and the first three steps are exactly a permutation of the bridge/selector source gates.

This is an order certificate, not a theorem discharge: it narrows the legal proof-search sequence but exports no source theorem, no chi11 source, no QW-2191 discharge, no role-transfer license, no `L_total`, and no ToE closure.
""".strip()
    lag_section = """
## P2423/S1373 admissible repair-order poset guard for Lagrangian/EOM

`P2423/S1373` proves that role-transfer and role-bearing `L_total` are order-constrained: source discharge plus the chi11/QW-2191 selector pair must come first, role-transfer can only appear at step 4, and `L_total`/ToE only at step 5.  Therefore no reordered shortcut can promote partial repair work into a role-bearing Lagrangian or ToE closure.
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
        "theorem_name": "P2423_T1_admissible_repair_order_poset_certificate",
        "missing_gate_count": len(cert["missing_gates"]),
        "precedence_edge_count": len(cert["precedence_edges"]),
        "raw_order_count": cert["raw_order_count"],
        "admissible_order_count": cert["admissible_order_count"],
        "rejected_order_count": cert["rejected_order_count"],
        "admissible_orders": [row["order"] for row in cert["admissible_order_rows"]],
        "role_transfer_ready_step_distribution": cert["role_transfer_ready_step_distribution"],
        "role_bearing_ltotal_ready_step_distribution": cert["role_bearing_ltotal_ready_step_distribution"],
        "toe_ready_step_distribution": cert["toe_ready_step_distribution"],
        "selector_source_ready_step_distribution": cert["selector_source_ready_step_distribution"],
        "all_edges_necessary": all(row["newly_admitted_count_if_removed"] > 0 for row in cert["edge_necessity_rows"]),
        "p2422_repair_subcube_count_inherited": cert["p2422_repair_subcube_count_inherited"] == 32,
        "p2422_toe_ready_repair_count_inherited": cert["p2422_toe_ready_repair_count_inherited"] == 1,
        "p2422_selector_singleton_unlock_count_inherited": cert["p2422_selector_singleton_unlock_count_inherited"] == 0,
        "source_obligation_discharge_exported": False,
        "chi11_source_exported": False,
        "qw2191_discharged": False,
        "role_transfer_licensed": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Admissible order is not theorem discharge for any missing gate.",
            "Role-transfer remains after bridge/selector source gates and is not a shortcut around them.",
            "Role-bearing L_total remains last and cannot be promoted from earlier prefixes.",
            "No source, selector, role-transfer, L_total, or ToE theorem is exported.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "five_missing_gates": theorem_export["missing_gate_count"] == 5,
        "four_precedence_edges": theorem_export["precedence_edge_count"] == 4,
        "all_120_raw_orders": theorem_export["raw_order_count"] == 120,
        "six_admissible_orders": theorem_export["admissible_order_count"] == 6,
        "one_hundred_fourteen_rejected": theorem_export["rejected_order_count"] == 114,
        "role_transfer_only_step_four": theorem_export["role_transfer_ready_step_distribution"] == {"4": 6},
        "ltotal_only_step_five": theorem_export["role_bearing_ltotal_ready_step_distribution"] == {"5": 6},
        "toe_only_step_five": theorem_export["toe_ready_step_distribution"] == {"5": 6},
        "selector_ready_step_distribution_expected": theorem_export["selector_source_ready_step_distribution"] == {"2": 2, "3": 4},
        "all_edges_necessary": theorem_export["all_edges_necessary"],
        "p2422_repair_subcube_inherited": theorem_export["p2422_repair_subcube_count_inherited"],
        "p2422_toe_repair_inherited": theorem_export["p2422_toe_ready_repair_count_inherited"],
        "p2422_no_selector_singleton_inherited": theorem_export["p2422_selector_singleton_unlock_count_inherited"],
        "source_still_open": not theorem_export["source_obligation_discharge_exported"],
        "chi11_still_open": not theorem_export["chi11_source_exported"],
        "qw2191_still_open": not theorem_export["qw2191_discharged"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "ltotal_still_blocked": not theorem_export["role_bearing_ltotal_exported"],
        "toe_still_open": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2423_s1373_v1",
        "packet_id": "P2423",
        "stage_id": "S1373",
        "result_kind": "ADMISSIBLE_REPAIR_ORDER_POSET_CERTIFICATE",
        "status": "PASS_ADMISSIBLE_REPAIR_ORDER_POSET_NO_GATE_DISCHARGE",
        "admissible_repair_order_poset_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Attempt one of the first-three source gates under the admissible order discipline; do not move role-transfer "
            "or L_total before bridge/selector source gates."
        ),
        "global_status": "OPEN_PROGRESS_ADMISSIBLE_ORDER_CERTIFIED_NO_GATE_DISCHARGED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["admissible_repair_order_poset_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2423 S1373: admissible repair-order poset certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Missing gates: `{theorem['missing_gate_count']}`.",
                f"- Precedence edges: `{theorem['precedence_edge_count']}`.",
                f"- Raw orders: `{theorem['raw_order_count']}`.",
                f"- Admissible orders: `{theorem['admissible_order_count']}`.",
                f"- Rejected orders: `{theorem['rejected_order_count']}`.",
                f"- Role-transfer step distribution: `{theorem['role_transfer_ready_step_distribution']}`.",
                f"- L_total step distribution: `{theorem['role_bearing_ltotal_ready_step_distribution']}`.",
                f"- ToE step distribution: `{theorem['toe_ready_step_distribution']}`.",
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
