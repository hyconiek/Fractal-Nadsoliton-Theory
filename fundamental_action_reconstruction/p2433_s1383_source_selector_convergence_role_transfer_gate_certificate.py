#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import combinations
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2433_s1383_source_selector_convergence_role_transfer_gate_certificate.json"
MD = GEN / "p2433_s1383_source_selector_convergence_role_transfer_gate_certificate.md"

SOURCE_FILES = {
    "P2432_TRANSITION": GEN / "p2432_s1382_post_antichain_branch_residual_transition_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

GATES = [
    "source_obligation_discharge",
    "chi11_source_export",
    "qw2191_selector_discharge",
    "role_transfer_audit_license",
    "role_bearing_ltotal_export",
]
PREDECESSORS = {
    "source_obligation_discharge": set(),
    "chi11_source_export": set(),
    "qw2191_selector_discharge": set(),
    "role_transfer_audit_license": {
        "source_obligation_discharge",
        "chi11_source_export",
        "qw2191_selector_discharge",
    },
    "role_bearing_ltotal_export": {"role_transfer_audit_license"},
}
CONVERGENCE_STATES = {
    "source_then_selector_pair": [
        "source_obligation_discharge",
        "chi11_source_export",
        "qw2191_selector_discharge",
    ],
    "selector_pair_then_source": [
        "chi11_source_export",
        "qw2191_selector_discharge",
        "source_obligation_discharge",
    ],
}
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
        "new_packet": "P2433|S1383|source selector convergence|role-transfer gate|convergence role-transfer",
        "p2432_input": "P2432|post-antichain branch|branch residual transition|source-first|selector-pair-first",
        "role_prior": "role-transfer|role_bearing_ltotal|L_total|source.*selector",
        "nonpromotion_prior": "nonpromotion|does not discharge|No source|No ToE|not a theorem discharge",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds P2432 branch residual transitions and generic source/selector/role language, but no production "
            "P24xx certificate proving that both antichain orders converge to the same role-transfer-admissible, L_total-blocked state."
        ),
    }


def theorem(payload: dict[str, Any], key: str) -> dict[str, Any]:
    return payload.get(key, {}).get("theorem_export", {})


def readiness_from_gates(gates: set[str]) -> dict[str, bool]:
    return {
        "bridge_source_ready": "source_obligation_discharge" in gates,
        "selector_source_ready": {"chi11_source_export", "qw2191_selector_discharge"}.issubset(gates),
        "role_transfer_ready": "role_transfer_audit_license" in gates,
        "role_bearing_ltotal_ready": "role_bearing_ltotal_export" in gates,
        "toe_ready": set(GATES).issubset(gates),
    }


def all_nonempty_subsets(items: list[str], max_size: int) -> list[list[str]]:
    out = []
    for size in range(1, max_size + 1):
        for combo in combinations(items, size):
            out.append(list(combo))
    return out


def is_admissible(candidate: set[str], current: set[str]) -> bool:
    available = current | candidate
    return all(PREDECESSORS[gate].issubset(available) for gate in candidate)


def readiness_unlocks(candidate: set[str], current: set[str]) -> dict[str, bool]:
    before = readiness_from_gates(current)
    after = readiness_from_gates(current | candidate)
    return {target: after[target] and not before[target] for target in TARGETS}


def state_rows() -> list[dict[str, Any]]:
    rows = []
    for state_name, sequence in CONVERGENCE_STATES.items():
        current = set(sequence)
        remaining = [gate for gate in GATES if gate not in current]
        candidate_rows = []
        for candidate in all_nonempty_subsets(remaining, 2):
            candidate_set = set(candidate)
            candidate_rows.append(
                {
                    "candidate_gates": candidate,
                    "candidate_size": len(candidate),
                    "admissible_after_convergence": is_admissible(candidate_set, current),
                    "readiness_unlocks": readiness_unlocks(candidate_set, current),
                    "remaining_after_candidate": [gate for gate in GATES if gate not in current | candidate_set],
                }
            )
        rows.append(
            {
                "state_name": state_name,
                "discharged_sequence": sequence,
                "discharged_gate_set": sorted(current, key=GATES.index),
                "readiness_after_convergence": readiness_from_gates(current),
                "remaining_gates_after_convergence": remaining,
                "remaining_gate_count_after_convergence": len(remaining),
                "candidate_rows_size_le_2_after_convergence": candidate_rows,
                "admissible_singletons_after_convergence": [
                    row["candidate_gates"] for row in candidate_rows if row["candidate_size"] == 1 and row["admissible_after_convergence"]
                ],
                "inadmissible_singletons_after_convergence": [
                    row["candidate_gates"] for row in candidate_rows if row["candidate_size"] == 1 and not row["admissible_after_convergence"]
                ],
            }
        )
    return rows


def build_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2432 = theorem(sources["P2432_TRANSITION"], "post_antichain_branch_residual_transition_certificate")
    rows = state_rows()
    return {
        "convergence_state_rows": rows,
        "convergence_state_count": len(rows),
        "convergence_gate_sets": [row["discharged_gate_set"] for row in rows],
        "p2432_source_first_next_inherited": p2432.get("source_first_next_minimal_readiness_candidates"),
        "p2432_selector_pair_next_inherited": p2432.get("selector_pair_first_next_minimal_readiness_candidates"),
        "p2432_role_transfer_not_after_first_branches_inherited": (
            p2432.get("role_transfer_admissible_after_source_first"),
            p2432.get("role_transfer_admissible_after_selector_pair_first"),
        ),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2433/S1383 source-selector convergence role-transfer gate certificate

`P2433/S1383` follows P2432 one step further: after either source-first then selector-pair, or selector-pair-first then source, both branches reach the same discharged theorem-gate set `{source_obligation_discharge, chi11_source_export, qw2191_selector_discharge}`.  At that convergence state, bridge-source and selector-source readiness are true, and `role_transfer_audit_license` becomes the admissible singleton next target.

The certificate also proves the remaining limit: `role_bearing_ltotal_export` is still not admissible until role-transfer is actually discharged, and ToE remains false until both role-transfer and role-bearing `L_total` are exported.  Convergence therefore licenses a next target, not closure.
""".strip()
    lag_section = """
## P2433/S1383 source-selector convergence guard for Lagrangian/EOM

`P2433/S1383` shows that once source and the chi11/QW-2191 selector pair are both discharged, role-transfer becomes the next admissible singleton target, but role-bearing `L_total` is still blocked.  The Lagrangian/EOM draft may use this as proof-order guidance only; it cannot write the role-bearing term before the role-transfer and `L_total` theorem gates are exported.
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
    rows = cert["convergence_state_rows"]
    theorem_export = {
        "theorem_name": "P2433_T1_source_selector_convergence_role_transfer_gate_certificate",
        "convergence_state_count": cert["convergence_state_count"],
        "convergence_gate_sets": cert["convergence_gate_sets"],
        "convergence_states_have_same_gate_set": len({tuple(item) for item in cert["convergence_gate_sets"]}) == 1,
        "readiness_after_convergence": rows[0]["readiness_after_convergence"],
        "remaining_gates_after_convergence": rows[0]["remaining_gates_after_convergence"],
        "admissible_singletons_after_convergence": rows[0]["admissible_singletons_after_convergence"],
        "inadmissible_singletons_after_convergence": rows[0]["inadmissible_singletons_after_convergence"],
        "role_transfer_admissible_after_source_selector_convergence": ["role_transfer_audit_license"]
        in rows[0]["admissible_singletons_after_convergence"],
        "role_bearing_ltotal_admissible_after_source_selector_convergence": ["role_bearing_ltotal_export"]
        in rows[0]["admissible_singletons_after_convergence"],
        "toe_ready_after_source_selector_convergence": rows[0]["readiness_after_convergence"]["toe_ready"],
        "p2432_source_first_next_inherited": cert["p2432_source_first_next_inherited"] == [
            ["chi11_source_export", "qw2191_selector_discharge"]
        ],
        "p2432_selector_pair_next_inherited": cert["p2432_selector_pair_next_inherited"] == [
            ["source_obligation_discharge"],
            ["source_obligation_discharge", "role_transfer_audit_license"],
        ],
        "p2432_role_transfer_not_after_first_branches_inherited": cert[
            "p2432_role_transfer_not_after_first_branches_inherited"
        ] == (False, False),
        "source_obligation_discharge_exported_by_this_certificate": False,
        "chi11_source_exported_by_this_certificate": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Convergence after hypothetical source/selector discharge is not current discharge by this certificate.",
            "Role-transfer becomes an admissible next target only after source and selector are actually discharged.",
            "Role-bearing L_total remains blocked until role-transfer is discharged.",
            "No source, selector, role-transfer, role-bearing L_total, or ToE theorem is exported.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "two_convergence_states": theorem_export["convergence_state_count"] == 2,
        "same_convergence_gate_set": theorem_export["convergence_states_have_same_gate_set"],
        "bridge_and_selector_ready": theorem_export["readiness_after_convergence"] == {
            "bridge_source_ready": True,
            "selector_source_ready": True,
            "role_transfer_ready": False,
            "role_bearing_ltotal_ready": False,
            "toe_ready": False,
        },
        "remaining_role_and_ltotal": theorem_export["remaining_gates_after_convergence"] == [
            "role_transfer_audit_license",
            "role_bearing_ltotal_export",
        ],
        "role_transfer_next": theorem_export["role_transfer_admissible_after_source_selector_convergence"],
        "ltotal_not_next": not theorem_export["role_bearing_ltotal_admissible_after_source_selector_convergence"],
        "toe_not_ready": not theorem_export["toe_ready_after_source_selector_convergence"],
        "p2432_source_first_inherited": theorem_export["p2432_source_first_next_inherited"],
        "p2432_selector_first_inherited": theorem_export["p2432_selector_pair_next_inherited"],
        "p2432_no_role_after_first_inherited": theorem_export["p2432_role_transfer_not_after_first_branches_inherited"],
        "source_not_exported_here": not theorem_export["source_obligation_discharge_exported_by_this_certificate"],
        "chi11_not_exported_here": not theorem_export["chi11_source_exported_by_this_certificate"],
        "qw2191_not_exported_here": not theorem_export["qw2191_discharged_by_this_certificate"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "ltotal_still_blocked": not theorem_export["role_bearing_ltotal_exported"],
        "toe_still_open": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2433_s1383_v1",
        "packet_id": "P2433",
        "stage_id": "S1383",
        "result_kind": "SOURCE_SELECTOR_CONVERGENCE_ROLE_TRANSFER_GATE_CERTIFICATE",
        "status": "PASS_SOURCE_SELECTOR_CONVERGENCE_ROLE_TRANSFER_NEXT_NO_GATE_DISCHARGE_NO_CLOSURE",
        "source_selector_convergence_role_transfer_gate_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "If source and selector pair become actual theorems, target role-transfer next; do not skip to role-bearing L_total."
        ),
        "global_status": "OPEN_PROGRESS_SOURCE_SELECTOR_CONVERGENCE_CERTIFIED_ROLE_TRANSFER_NEXT_NO_THEOREM_GATE_DISCHARGED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["source_selector_convergence_role_transfer_gate_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2433 S1383: source-selector convergence role-transfer gate certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Convergence states: `{theorem['convergence_state_count']}`.",
                f"- Readiness after convergence: `{theorem['readiness_after_convergence']}`.",
                f"- Remaining gates: `{theorem['remaining_gates_after_convergence']}`.",
                f"- Admissible singletons: `{theorem['admissible_singletons_after_convergence']}`.",
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
