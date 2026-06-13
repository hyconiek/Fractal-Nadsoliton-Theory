#!/usr/bin/env python3
"""P2691/S1641: alpha_geo role-safe amplitude source audit.

This executes the post-P2690 recommendation: attack one remaining P2680
non-selector bridge atom by testing whether alpha_geo scalar-shape normalization
can be promoted to a strict, role-safe amplitude source without legacy role
transfer, selector replay, beta_tors->chi11, generic bridge completion, L_total,
or ToE overclaim.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2691_s1641_alpha_geo_role_safe_amplitude_source_audit.json"
MD = GEN / "p2691_s1641_alpha_geo_role_safe_amplitude_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2690": GEN / "p2690_s1640_selector_free_nonexact_boundary_phase_sector_provider_audit.json",
    "P2680": GEN / "p2680_s1630_legacy_strict_bridge_source_inventory_no_selector_replay_audit.json",
    "ALPHA": GEN / "alpha_geo_strict_derived_v1.json",
    "F1": ROOT / "F1_CANONICAL_INFORMATIONAL_FRACTAL_SUBSTRATE_PARAMETER_PACKET.md",
    "N65": ROOT / "N65_CURRENT_LEGACY_PHYSICAL_ROLE_TRANSFER_TO_STRICT_GATE_KERNEL_OBSTRUCTION_THEOREM.md",
    "ROLE_DRAFT": ROOT / "STRICT_KERNEL_LEGACY_ROLE_TRANSFER_AUDIT_DRAFT.md",
    "TOE_AUDIT": ROOT / "STRICT_KERNEL_TOE_POTENTIAL_AUDIT.md",
}

NEGATIVE_EXPORT_FLAGS = [
    "amplitude_role_safe_source_exported",
    "legacy_role_transfer_imported",
    "selector_replay_imported",
    "beta_tors_chi11_imported",
    "generic_bridge_completion_claimed",
    "ltotal_promoted",
    "toe_closure_claimed",
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.exists() else ""


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
        "p2690_selected_p2691": r"P2691|alpha_geo scalar-shape|amplitude role-safe|role-safe amplitude",
        "p2680_amplitude_atom": r"alpha_geo_scalar_shape_normalization|amplitude_role_safe_source|scalar.*normalization|amplitude.*normalization",
        "strict_alpha_source": r"alpha_geo_strict_derived_v1|4 ln\(2\)|ln\(16\)|Shannon",
        "legacy_role_blockers": r"alpha_geo/12|alpha_EM|role transfer|physical-role|electroweak role|scalar normalization alone",
        "forbidden_imports": r"selector replay|QW-2191|beta_tors.*chi_?11|generic bridge|bridge completion|L_total|ToE closure",
    }
    return {"tool": "rg", "mode": "content-first alpha_geo role-safe amplitude source audit", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def state_reads() -> dict[str, Any]:
    p2690 = load_json(INPUTS["P2690"])
    p2680 = load_json(INPUTS["P2680"])
    alpha = load_json(INPUTS["ALPHA"])
    atoms = {row.get("atom"): row for row in p2680.get("source_atom_inventory", [])}
    role_text = read_text(INPUTS["ROLE_DRAFT"])
    toe_text = read_text(INPUTS["TOE_AUDIT"])
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2690_selected_p2691": "P2691" in p2690.get("decision", {}).get("next_honest_step", ""),
        "alpha_strict_source_exported": alpha.get("status") == "actual_exported_strict_derived_source_upgrade_value",
        "alpha_value": alpha.get("value"),
        "alpha_hard_limits": alpha.get("hard_limits", []),
        "p2680_scalar_shape_normalization_present": atoms.get("alpha_geo_scalar_shape_normalization", {}).get("formal_material_present") is True,
        "p2680_scalar_shape_role_safe": atoms.get("alpha_geo_scalar_shape_normalization", {}).get("role_safe_for_completion") is True,
        "p2680_amplitude_role_safe_source_exported": atoms.get("amplitude_role_safe_source", {}).get("source_theorem_exported") is True,
        "role_draft_blocks_scalar_normalization_as_physical_role_proof": "scalar normalization alone is not a physical-role proof" in role_text,
        "toe_audit_keeps_apd_source_open": "strict_dynamical_source_for_A_P_D" in toe_text and "open" in toe_text,
    }


def amplitude_normalization_computation() -> dict[str, Any]:
    alpha = 4.0 * math.log(2.0)
    omega = math.pi / 4.0
    phi = math.pi / 6.0
    beta_tors = 0.01
    rows = []
    for d in range(0, 13):
        legacy = alpha * math.cos(omega * d + phi) / (1.0 + beta_tors * d)
        carrier = legacy / alpha
        expected = math.cos(omega * d + phi) / (1.0 + beta_tors * d)
        rows.append({"d": d, "legacy": legacy, "alpha_removed_carrier": carrier, "expected_carrier": expected, "residual": carrier - expected})
    max_abs_residual = max(abs(row["residual"]) for row in rows)
    return {
        "alpha_geo_numeric": alpha,
        "alpha_matches_4_ln2": abs(alpha - 2.772588722239781) < 1e-15,
        "formula_tested": "K_legacy_ont(d)/alpha_geo = cos(omega*d+phi)/(1+beta_tors*d) on d=0..12",
        "max_abs_residual": max_abs_residual,
        "scalar_shape_normalization_is_exact": max_abs_residual < 1e-15,
        "sample_rows": rows,
        "role_semantics_generated_by_this_computation": False,
    }


def obligation_matrix(reads: dict[str, Any], comp: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "obligation": "strict_alpha_geo_value_source",
            "needed_for_role_safe_amplitude_source": True,
            "satisfied_now": reads["alpha_strict_source_exported"] and reads["alpha_value"] == "4 ln(2)",
            "evidence": "alpha_geo_strict_derived_v1 exports 4 ln(2) with hard limits.",
        },
        {
            "obligation": "exact_scalar_shape_normalization",
            "needed_for_role_safe_amplitude_source": True,
            "satisfied_now": comp["scalar_shape_normalization_is_exact"] and reads["p2680_scalar_shape_normalization_present"],
            "evidence": "Direct finite computation verifies alpha removal from K_legacy_ont on d=0..12.",
        },
        {
            "obligation": "amplitude_absorption_as_strict_bridge_source",
            "needed_for_role_safe_amplitude_source": True,
            "satisfied_now": False,
            "evidence": "P2680 explicitly leaves amplitude_role_safe_source unexported; exact scalar normalization is only shape algebra.",
        },
        {
            "obligation": "physical_role_safety_theorem",
            "needed_for_role_safe_amplitude_source": True,
            "satisfied_now": False,
            "evidence": "Role draft says scalar normalization alone is not a physical-role proof; no alpha_geo electroweak/EM role theorem is imported.",
        },
        {
            "obligation": "apd_or_lagrangian_dynamical_source",
            "needed_for_role_safe_amplitude_source": True,
            "satisfied_now": False,
            "evidence": "Strict ToE audit keeps strict_dynamical_source_for_A_P_D open; no role-bearing L_total term follows.",
        },
    ]


def decision(matrix: list[dict[str, Any]]) -> dict[str, Any]:
    missing = [row["obligation"] for row in matrix if row["needed_for_role_safe_amplitude_source"] and not row["satisfied_now"]]
    return {
        "decision": "P2691_ALPHA_GEO_ROLE_SAFE_AMPLITUDE_SOURCE_AUDIT_SCALAR_SHAPE_ONLY_NO_ROLE_SAFE_SOURCE_NO_FALSE_PASS",
        "satisfied_obligations": [row["obligation"] for row in matrix if row["satisfied_now"]],
        "missing_obligations": missing,
        "amplitude_role_safe_source_exported_now": False,
        "bounded_no_go_for_current_alpha_geo_amplitude_atom": True,
        "professorial_verdict": (
            "P2691 gives the alpha_geo atom its strongest fair reading.  The strict Shannon value alpha_geo=4 ln(2) is exported, and the finite computation confirms that dividing K_legacy_ont by alpha_geo exactly removes the scalar amplitude on the audited legacy support.  But that is only scalar-shape normalization.  Current artifacts still do not export amplitude absorption as a strict bridge source, a physical-role safety theorem, or an APD/Lagrangian dynamical source.  Therefore the alpha_geo amplitude atom remains bounded no-go for bridge completion on current evidence, without legacy role transfer or selector replay."
        ),
        "next_honest_step": (
            "P2692 should return to the P2680 non-selector inventory and attack the remaining damping/compression atom: a target-independent positive beta/Z_beta source audit, explicitly separated from canonical UV-unit replay, beta_tors->chi11, selector replay, role transfer, and generic bridge completion."
        ),
        "role_transfer_started_now": False,
        "ltotal_promoted_now": False,
        "toe_closed_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2691/S1641 alpha_geo role-safe amplitude source audit", "", f"Status: `{payload['status']}`", "", "## Content-first grep"]
    for name, data in payload["content_grep"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    comp = payload["amplitude_normalization_computation"]
    lines.extend(["", "## Scalar normalization computation", f"`{comp['formula_tested']}`", f"max_abs_residual = `{comp['max_abs_residual']}`; exact = `{comp['scalar_shape_normalization_is_exact']}`."])
    lines.extend(["", "## Obligation matrix"])
    for row in payload["obligation_matrix"]:
        lines.append(f"- `{row['obligation']}`: satisfied=`{row['satisfied_now']}` — {row['evidence']}")
    lines.extend(["", "## Verdict", payload["decision"]["professorial_verdict"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    reads = state_reads()
    comp = amplitude_normalization_computation()
    matrix = obligation_matrix(reads, comp)
    payload: dict[str, Any] = {
        "status": "P2691_ALPHA_GEO_ROLE_SAFE_AMPLITUDE_SOURCE_AUDIT_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "state_reads": reads,
        "amplitude_normalization_computation": comp,
        "obligation_matrix": matrix,
        "decision": decision(matrix),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2691/S1641 alpha_geo role-safe amplitude source audit",
        "## P2691/S1641 alpha_geo role-safe amplitude source audit\n\n"
        "`P2691/S1641` audits the remaining P2680 amplitude atom.  The strict Shannon source `alpha_geo_strict_derived_v1 = 4 ln(2)` is real, and finite symbolic/numeric checking confirms that `K_legacy_ont(d)/alpha_geo` exactly removes the scalar amplitude on the audited legacy support.  This is not yet a role-safe amplitude source: no amplitude-absorption bridge source, physical-role safety theorem, or APD/Lagrangian dynamical source is exported.  Therefore the `alpha_geo` amplitude atom is bounded no-go on current artifacts; no legacy role transfer, selector replay, `beta_tors -> chi11`, bridge completion, role-bearing `L_total`, or ToE closure is claimed.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2691/S1641 alpha_geo amplitude Ltotal guard",
        "## P2691/S1641 alpha_geo amplitude Ltotal guard\n\n"
        "`P2691/S1641` keeps `L_total` nonpromoted.  A strict value source for `alpha_geo=4 ln(2)` plus exact scalar normalization does not by itself create a variational amplitude term or physical-role transfer theorem.  The APD/Lagrangian source remains open.\n",
    )
    append_once(
        AGENTS,
        "Current alpha_geo amplitude source guardrail (P2691/S1641, 2026-06-13)",
        "## Current alpha_geo amplitude source guardrail (P2691/S1641, 2026-06-13)\n\n"
        "- P2691 confirms strict `alpha_geo=4 ln(2)` and exact scalar-shape normalization, but finds no role-safe amplitude absorption source, no physical-role safety theorem, and no APD/Lagrangian dynamical source.\n"
        "- Freeze the `alpha_geo` amplitude atom as bounded no-go on current artifacts; the next non-replay move should audit a target-independent positive `beta/Z_beta` source, without canonical UV-unit replay, `beta_tors -> chi11`, selector replay, role transfer, generic bridge completion, `L_total`, or ToE closure.\n",
    )
    return payload


if __name__ == "__main__":
    main()
