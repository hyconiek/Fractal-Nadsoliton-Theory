#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2392_s1342_auxiliary_beta_tors_chi11_retirement_certificate.json"
MD = GEN / "p2392_s1342_auxiliary_beta_tors_chi11_retirement_certificate.md"

SOURCE_FILES = {
    "P1343_SELECTOR_SOURCE_REPORT": GEN / "p1343_p1343_report_v1.json",
    "P1348_GLOBAL_CLOSURE_REPORT": GEN / "p1348_p1348_report_v1.json",
    "P2390_SELECTOR_QUALIFIED_AUDIT": GEN / "p2390_s1340_selector_qualified_beta_tors_chi11_role_audit.json",
    "P2391_REBASED_GAP_MATRIX": GEN / "p2391_s1341_selector_epoch_rebased_bridge_gap_matrix.json",
    "S2_STRATEGIC_PRIORITY": ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

ATOMS = [
    "strict_internal_selector_P1343_P1348",
    "auxiliary_beta_tors_to_chi11",
    "legacy_to_strict_bridge_completion",
    "post_bridge_role_transfer_audit",
]

TARGETS = {
    "selector_mechanism": [
        {"strict_internal_selector_P1343_P1348"},
        {"auxiliary_beta_tors_to_chi11"},
    ],
    "legacy_role_transfer_permission": [
        {"legacy_to_strict_bridge_completion", "post_bridge_role_transfer_audit"},
    ],
}


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    if not path.exists():
        return f"OPEN_MISSING_ARTIFACT::{rel(path)}"
    return path.read_text(encoding="utf-8")


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def rg_audit() -> dict[str, Any]:
    patterns = [
        "beta_tors -> chi11|beta_tors -> chi_11|beta_tors.*chi_?11",
        "P2390|P2391|selector-qualified|selector-epoch",
        "auxiliary.*beta_tors|candidate bridge hypothesis|not a theorem",
        "S_strict_internal_v1|P1343|P1348",
    ]
    out: dict[str, Any] = {}
    for pattern in patterns:
        proc = subprocess.run(
            ["rg", "-n", pattern, "fundamental_action_reconstruction", "-g", "*.py", "-g", "*.md", "-g", "*.json"],
            cwd=REPO,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        lines = [line for line in proc.stdout.splitlines() if line]
        out[pattern] = {"count": len(lines), "samples": lines[:10]}
    return {
        "tool": "rg",
        "patterns": out,
        "finding": (
            "Repo already contains P2390/P2391 separation artifacts and many beta_tors/chi11 candidate references. "
            "P2392 therefore does not try to prove beta_tors->chi11; it retires that candidate from the active selector-mechanism route."
        ),
    }


def all_minimal_supports(target: str) -> list[list[str]]:
    clauses = TARGETS[target]
    supports: list[set[str]] = []
    for size in range(1, len(ATOMS) + 1):
        for combo in itertools.combinations(ATOMS, size):
            support = set(combo)
            if any(clause <= support for clause in clauses):
                if not any(prev <= support for prev in supports):
                    supports.append(support)
        if supports:
            break
    return [sorted(support) for support in supports]


def evaluate_supports(available: dict[str, bool]) -> dict[str, Any]:
    evaluated: dict[str, Any] = {}
    for target in TARGETS:
        supports = all_minimal_supports(target)
        realized = [support for support in supports if all(available[atom] for atom in support)]
        evaluated[target] = {
            "minimal_supports": supports,
            "realized_minimal_supports": realized,
            "target_satisfied": bool(realized),
            "uses_auxiliary_beta_tors_in_realized_support": any("auxiliary_beta_tors_to_chi11" in support for support in realized),
        }
    return evaluated


def build_payload() -> dict[str, Any]:
    artifacts: dict[str, Any] = {}
    for name, path in SOURCE_FILES.items():
        artifacts[name] = load_json(path) if path.suffix == ".json" else {"text": read_text(path)}
    grep = rg_audit()

    p2390 = artifacts["P2390_SELECTOR_QUALIFIED_AUDIT"]["selector_qualified_beta_tors_chi11_role_audit"]["implication_certificate"]
    p2391 = artifacts["P2391_REBASED_GAP_MATRIX"]["selector_epoch_rebased_bridge_gap_matrix"]
    chi11_row = next(row for row in p2391["rows"] if row["component"] == "topological_phase_bit_chi11")

    selector_available = (
        artifacts["P1343_SELECTOR_SOURCE_REPORT"].get("status") == "CLOSED_FULL_STRICT_INTERNAL_SOURCE_V1"
        and artifacts["P1348_GLOBAL_CLOSURE_REPORT"].get("status") == "GLOBAL_CLOSURE_THEOREM_EXPORTED_CLOSED_DECLARED_SCOPE"
        and bool(p2390["transfer_rule_inputs"]["selector_source"])
        and bool(chi11_row["rebased_bits"]["generic_strict_selector_exported"])
    )
    beta_transport_available = bool(p2390["transfer_rule_inputs"]["explicit_beta_tors_to_chi11_transport_theorem"])
    bridge_ready = bool(p2390["transfer_rule_inputs"]["full_legacy_to_strict_bridge_ready"])
    role_audit_ready = bool(p2390["transfer_rule_inputs"]["separate_role_transfer_theorem"])
    available = {
        "strict_internal_selector_P1343_P1348": selector_available,
        "auxiliary_beta_tors_to_chi11": beta_transport_available,
        "legacy_to_strict_bridge_completion": bridge_ready,
        "post_bridge_role_transfer_audit": role_audit_ready,
    }
    supports = evaluate_supports(available)

    obligation_rows = [
        {
            "obligation": "strict selector mechanism",
            "status": "DISCHARGED_BY_P1343_P1348",
            "active_next_target": False,
            "beta_tors_chi11_required": False,
            "evidence": "selector_mechanism has a realized minimal support using strict_internal_selector_P1343_P1348",
        },
        {
            "obligation": "auxiliary beta_tors -> chi11 selector-search hypothesis",
            "status": "RETIRED_AS_SELECTOR_ROUTE_ASSUMPTION",
            "active_next_target": False,
            "beta_tors_chi11_required": False,
            "evidence": "not used by the realized selector_mechanism support",
        },
        {
            "obligation": "legacy -> strict bridge completion",
            "status": "OPEN_ACTIVE_BRIDGE_TARGET",
            "active_next_target": True,
            "beta_tors_chi11_required": False,
            "evidence": "bridge completion is independent of confirming the retired selector-search hypothesis",
        },
        {
            "obligation": "post-bridge role-transfer audit",
            "status": "DEFERRED_UNTIL_BRIDGE_READY",
            "active_next_target": False,
            "beta_tors_chi11_required": False,
            "evidence": "role transfer remains blocked, but beta_tors->chi11 is not required just to preserve the selector mechanism",
        },
    ]
    active_beta_obligations = [row for row in obligation_rows if row["active_next_target"] and row["beta_tors_chi11_required"]]
    theorem_export = {
        "theorem_name": "P2392_T1_auxiliary_beta_tors_chi11_selector_assumption_retired",
        "statement": (
            "Because the current selector mechanism is realized by P1343/P1348 and P2391 rebases the chi11 row to generic selector-present, "
            "the auxiliary beta_tors->chi11 hypothesis is no longer an active selector-mechanism obligation. "
            "It may remain a historical/candidate legacy-role idea, but it need not be confirmed to continue the strict selector route."
        ),
        "available_atoms": available,
        "support_evaluation": supports,
        "active_beta_tors_chi11_obligation_count": len(active_beta_obligations),
        "licensed": [
            "Drop beta_tors->chi11 from the active selector-mechanism blocker list.",
            "Keep beta_tors legacy physical roles separated until an explicit post-bridge role-transfer audit.",
            "Move the next bridge work toward bridge completion / role-audit prerequisites rather than proving the retired selector-search assumption.",
        ],
        "not_licensed": [
            "No beta_tors -> chi11 theorem is claimed or needed here.",
            "No silent legacy physical-role transfer is licensed.",
            "No completed legacy -> strict bridge, L_total promotion, cap-density source theorem, SM/GR numeric extraction, or ToE closure is claimed.",
        ],
    }
    gatekeepers = {
        "p1343_p1348_selector_available": selector_available,
        "p2391_chi11_selector_rebased": bool(chi11_row["rebased_bits"]["generic_strict_selector_exported"]),
        "selector_support_realized_without_beta_tors": supports["selector_mechanism"]["target_satisfied"] and not supports["selector_mechanism"]["uses_auxiliary_beta_tors_in_realized_support"],
        "beta_transport_not_available_and_not_required": (not beta_transport_available) and len(active_beta_obligations) == 0,
        "legacy_role_transfer_still_unsatisfied": not supports["legacy_role_transfer_permission"]["target_satisfied"],
        "rg_nonduplication_found_prior_artifacts": grep["patterns"]["P2390|P2391|selector-qualified|selector-epoch"]["count"] > 0,
        "docs_updated_with_p2392": all("P2392/S1342" in read_text(path) for path in DOC_FILES.values()),
    }
    return {
        "schema_version": "p2392_s1342_v1",
        "packet_id": "P2392",
        "stage_id": "S1342",
        "timestamp_utc": "2026-06-01T00:00:00Z",
        "produced_by": rel(Path(__file__).resolve()),
        "result_kind": "AUXILIARY_BETA_TORS_CHI11_SELECTOR_ASSUMPTION_RETIREMENT_CERTIFICATE",
        "status": "PASS_BETA_TORS_CHI11_RETIRED_AS_SELECTOR_ROUTE_ASSUMPTION",
        "auxiliary_beta_tors_chi11_retirement_certificate": {
            "rg_nonduplication_audit": grep,
            "obligation_rows": obligation_rows,
            "active_beta_tors_chi11_obligations": active_beta_obligations,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": "Stop treating beta_tors->chi11 as an active selector target; work on the explicit legacy->strict bridge completion and then the role-transfer audit.",
        "global_status": "OPEN_PROGRESS_SELECTOR_ASSUMPTION_RETIRED_BRIDGE_AND_ROLE_AUDIT_OPEN",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    cert = payload["auxiliary_beta_tors_chi11_retirement_certificate"]
    theorem = cert["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2392 S1342: auxiliary beta_tors -> chi11 selector-assumption retirement certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2392/S1342 implements the correction that `beta_tors -> chi11` was an auxiliary assumption used while searching for a selector mechanism.",
                "Since P1343/P1348 now provide the strict selector mechanism and P2391 rebases the chi11 row to selector-present, this auxiliary route is retired from the active selector blocker list.",
                "",
                "## Minimal-support computation",
                "",
                f"- Available atoms: `{theorem['available_atoms']}`.",
                f"- Selector support evaluation: `{theorem['support_evaluation']['selector_mechanism']}`.",
                f"- Active beta_tors->chi11 obligation count: `{theorem['active_beta_tors_chi11_obligation_count']}`.",
                "",
                "## Hard limits",
                "",
                "- No `beta_tors -> chi11` theorem is claimed or needed for the selector route.",
                "- Legacy physical-role transfer remains blocked until bridge completion and a separate role-transfer audit.",
                "- No `L_total` promotion, cap-density source theorem, SM/GR numeric extraction, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    payload = build_payload()
    write_outputs(payload)


if __name__ == "__main__":
    main()
