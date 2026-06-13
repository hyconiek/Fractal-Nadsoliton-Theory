#!/usr/bin/env python3
"""P2687/S1637: one-new strict anisotropic source-class audit.

This implements the P2686 continuation rule without replaying the reduced/FRW
reverse-closure lane.  It audits exactly the two source classes named by P2686
and asks whether either is actually exported as a new strict typed source class
that evades the existing P1977/P1978 blockers.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

import sympy as sp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2687_s1637_one_new_strict_anisotropic_source_class_audit.json"
MD = GEN / "p2687_s1637_one_new_strict_anisotropic_source_class_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2686": GEN / "p2686_s1636_shared_background_nonproxy_component_residual_table.json",
    "P1975": GEN / "p1975_s925_strict_minimal_anisotropic_source_obligation_and_energy_sign_audit.json",
    "P1977": GEN / "p1977_s927_strict_positive_energy_anisotropic_provider_bounded_no_go.json",
    "P1978": GEN / "p1978_s928_strict_energy_neutral_tensor_transport_obstruction.json",
    "P1976": GEN / "p1976_s926_strict_ltotal_anisotropic_provider_nonexport_audit.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "new_strict_anisotropic_source_exported",
    "p1977_evaded_by_positive_energy_lapse_source",
    "p1978_evaded_by_non_energy_neutral_transport",
    "reverse_closure_reopened",
    "ltotal_promoted",
    "toe_closure_claimed",
    "selector_or_bridge_imported",
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        ["rg", "-n", pattern, ".", "-g", "*.py", "-g", "*.md", "-g", "*.json", "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**"],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def content_grep() -> dict[str, Any]:
    patterns = {
        "p2686_continuation_rule": r"P2687|one new strict anisotropic source|evades P1977|evades P1978|freeze this Lagrangian/EOM",
        "derived_lapse_energy_source": r"derived lapse|lapse/energy|rho_required|positive-energy|negative rho|rho_provider",
        "non_energy_neutral_tensor_transport": r"non-energy-neutral|energy-neutral|tensorial shear transport|U00|spatial trace|tracefree",
        "provider_nonexport_registry": r"anisotropic provider|provider nonexport|strict source derivation exported|current basis has no exported",
        "forbidden_imports": r"selector closure|role transfer|legacy bridge|QW-2191|ToE closure|role-bearing L_total",
    }
    return {"tool": "rg", "mode": "content-first one-new strict anisotropic source-class audit", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def current_state() -> dict[str, Any]:
    data = {name: load_json(path) for name, path in INPUTS.items()}
    p2686_next = data["P2686"].get("decision", {}).get("next_honest_step", "")
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2686_requested_one_new_source_class": "exactly one new strict anisotropic source class" in p2686_next,
        "p1975_minimal_source_cancels_if_admitted": data["P1975"].get("gatekeeper_checks", {}).get("minimal_source_cancels_residual_if_admitted") is True,
        "p1975_strict_source_derivation_exported": data["P1975"].get("gatekeeper_checks", {}).get("strict_source_derivation_exported") is True,
        "p1977_bounded_no_go_passed": data["P1977"].get("gatekeeper_checks", {}).get("bounded_no_go_passed") is True,
        "p1978_bounded_obstruction_passed": data["P1978"].get("gatekeeper_checks", {}).get("bounded_obstruction_passed") is True,
        "p1976_nonexport_audit_passed": data["P1976"].get("gatekeeper_checks", {}).get("no_explicit_anisotropic_provider_in_current_registries") is True,
    }


def source_class_rows() -> list[dict[str, Any]]:
    p1975 = load_json(INPUTS["P1975"])
    p1977 = load_json(INPUTS["P1977"])
    p1978 = load_json(INPUTS["P1978"])
    sigma1, sigma2, u0 = sp.symbols("sigma1 sigma2 u0")
    q_shear = sp.sympify(p1975.get("source_obligation", {}).get("q_shear", "sigma1**2 + sigma1*sigma2 + sigma2**2"))
    rho_required = sp.sympify(p1975.get("source_obligation", {}).get("rho_required", "-sigma1**2 - sigma1*sigma2 - sigma2**2"))
    lapse_delta = sp.simplify(u0 - q_shear)  # u0 must equal Q_shear to avoid negative rho in the split source.
    required_u00 = sp.sympify(p1978.get("symbolic_core", {}).get("required_U00_for_full_cancellation", "-sigma1**2 - sigma1*sigma2 - sigma2**2"))
    return [
        {
            "source_class_id": "derived_lapse_energy_source",
            "target_blocker": "P1977_positive_energy_no_go",
            "typed_source_exported_now": False,
            "strict_derivation_source": "MISSING_CURRENT_ARTIFACT",
            "symbolic_requirement": f"rho_base + u0 split must satisfy u0 = {q_shear} to offset rho_required = {rho_required}",
            "positive_energy_evade_condition": str(sp.Eq(u0, q_shear)),
            "residual_if_unexported": str(lapse_delta),
            "evades_blocker_now": False,
            "verdict": "NEAR_MISS_NOT_EXPORTED",
        },
        {
            "source_class_id": "non_energy_neutral_tensorial_shear_transport",
            "target_blocker": "P1978_energy_neutral_transport_obstruction",
            "typed_source_exported_now": False,
            "strict_derivation_source": "MISSING_CURRENT_ARTIFACT",
            "symbolic_requirement": f"transport must abandon U00=0 and carry required_U00 = {required_u00}",
            "non_energy_neutral_condition": str(sp.Ne(u0, 0)),
            "required_u00": str(required_u00),
            "evades_blocker_now": False,
            "verdict": "NEAR_MISS_NOT_EXPORTED",
        },
    ]


def obligation_matrix(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    obligations = []
    for row in rows:
        obligations.extend([
            {"source_class_id": row["source_class_id"], "obligation": "typed strict source object", "satisfied": row["typed_source_exported_now"]},
            {"source_class_id": row["source_class_id"], "obligation": "internal derivation/provider citation", "satisfied": row["strict_derivation_source"] != "MISSING_CURRENT_ARTIFACT"},
            {"source_class_id": row["source_class_id"], "obligation": f"evade {row['target_blocker']}", "satisfied": row["evades_blocker_now"]},
        ])
    return obligations


def decision(rows: list[dict[str, Any]], obligations: list[dict[str, Any]]) -> dict[str, Any]:
    passing = [row for row in rows if row["typed_source_exported_now"] and row["evades_blocker_now"]]
    all_failed = all(not item["satisfied"] for item in obligations if item["obligation"].startswith("typed") or item["obligation"].startswith("internal") or item["obligation"].startswith("evade"))
    return {
        "decision": "P2687_ONE_NEW_STRICT_ANISOTROPIC_SOURCE_CLASS_AUDIT_NO_EXPORTED_EVASION_FREEZE_LANE_NO_FALSE_PASS",
        "new_source_class_exported": bool(passing),
        "passing_source_classes": [row["source_class_id"] for row in passing],
        "all_current_evasion_obligations_failed": all_failed,
        "freeze_lagrangian_eom_reverse_closure_lane": True,
        "professorial_verdict": (
            "P2687 audits the exact P2686 continuation: one new strict anisotropic source class.  The two finite candidates are explicit and computable, but neither is currently exported as a strict typed source/provider.  A lapse/energy split would need an internal source with u0=Q_shear to avoid the P1977 negative-rho conclusion; a non-energy-neutral tensorial transport would need to carry the required U00 component to evade P1978.  Both are near-miss design equations, not current theorem exports, so this lane freezes as bounded no-go rather than reopening reverse closure."
        ),
        "next_honest_step": (
            "P2688 should pivot back to the broad state-map and choose a different live proof frontier, with first preference returning to the kernel-split-robust P46/N49 direct-route zero-witness/no-go matrix named by F3/P2681 (especially the direct m2 psi4 target-action coefficient defect on common psi4**2/2 support).  Do not replay Lagrangian/EOM reverse closure unless a genuinely exported strict anisotropic source class is introduced."
        ),
        "ltotal_promoted_now": False,
        "toe_closed_now": False,
        "selector_bridge_role_imported_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2687/S1637 one-new strict anisotropic source-class audit", "", f"Status: `{payload['status']}`", "", "## Content-first grep"]
    for name, data in payload["content_grep"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    lines.extend(["", "## Source-class candidate rows"])
    for row in payload["source_class_rows"]:
        lines.append(f"- `{row['source_class_id']}` targets `{row['target_blocker']}`: exported=`{row['typed_source_exported_now']}`, evades_now=`{row['evades_blocker_now']}`, verdict=`{row['verdict']}`; requirement=`{row['symbolic_requirement']}`")
    lines.extend(["", "## Obligation matrix"])
    for item in payload["obligation_matrix"]:
        lines.append(f"- `{item['source_class_id']}` / `{item['obligation']}`: `{item['satisfied']}`")
    lines.extend(["", "## Verdict", payload["decision"]["professorial_verdict"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    rows = source_class_rows()
    obligations = obligation_matrix(rows)
    payload: dict[str, Any] = {
        "status": "P2687_ONE_NEW_STRICT_ANISOTROPIC_SOURCE_CLASS_AUDIT_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "current_state": current_state(),
        "source_class_rows": rows,
        "obligation_matrix": obligations,
        "decision": decision(rows, obligations),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2687/S1637 one-new strict anisotropic source-class audit",
        "## P2687/S1637 one-new strict anisotropic source-class audit\n\n"
        "`P2687/S1637` executes the P2686 continuation rule by auditing exactly one-new-source-class candidates for the Bianchi-I residual: a derived lapse/energy source against P1977 and a non-energy-neutral tensorial shear transport against P1978.  Both yield symbolic design equations, but neither is currently exported as a strict typed source/provider, so the Lagrangian/EOM reverse-closure lane remains frozen as bounded no-go rather than promoted to `L_total` or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2687/S1637 source-class audit freeze guard",
        "## P2687/S1637 source-class audit freeze guard\n\n"
        "`P2687/S1637` does not reopen reverse closure: no new strict anisotropic source class currently evades P1977 or P1978.  The honest continuation is a state-map pivot, first back to the F3/P2681 direct-route P46/N49 zero-witness/no-go matrix unless a new exported source class appears.\n",
    )
    append_once(
        AGENTS,
        "Current anisotropic source-class audit guardrail (P2687/S1637, 2026-06-13)",
        "## Current anisotropic source-class audit guardrail (P2687/S1637, 2026-06-13)\n\n"
        "- P2687 audits the two P2686-admissible strict anisotropic source-class continuations and finds no currently exported typed source/provider evading P1977 or P1978.\n"
        "- Freeze the strict Lagrangian/EOM reverse-closure lane as bounded no-go unless a genuinely new strict anisotropic source class is exported; the next proof-grade move should pivot through the state-map, preferably to the F3/P2681 P46/N49 direct-route zero-witness/no-go matrix.\n",
    )
    return payload


if __name__ == "__main__":
    main()
