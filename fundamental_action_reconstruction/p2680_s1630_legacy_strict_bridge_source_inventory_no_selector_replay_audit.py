#!/usr/bin/env python3
"""P2680/S1630: legacy->strict bridge-source inventory without selector replay.

P2679 forbids reopening selector, tau_src->pair12, and beta_tors->chi11 lanes
without a genuinely new object.  This packet follows the admissible pivot:
a bridge-source inventory for amplitude/normalization and damping/compression
atoms in the legacy -> strict kernel completion problem.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2680_s1630_legacy_strict_bridge_source_inventory_no_selector_replay_audit.json"
MD = GEN / "p2680_s1630_legacy_strict_bridge_source_inventory_no_selector_replay_audit.md"
P2679 = GEN / "p2679_s1629_reopen_repetition_gate_and_bridge_pivot_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NEGATIVE_EXPORT_FLAGS = [
    "full_legacy_strict_bridge_completed",
    "amplitude_role_safe_source_exported",
    "positive_beta_renormalization_source_exported",
    "canonical_length_uv_unit_exported",
    "selector_replay_used",
    "beta_tors_chi11_reopened",
    "role_transfer_theorem_exported",
    "role_bearing_ltotal_reenabled",
    "q_w_2191_discharged",
    "toe_closure_claimed",
]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:60]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "legacy_strict_kernel_content": "K_legacy_ont|K_strict_gate|legacy -> strict|legacy-to-strict|completion bridge",
        "amplitude_normalization_content": "alpha_geo|amplitude.*normalization|scalar.*normalization|amplitude absorption|amplitude",
        "damping_compression_content": "beta_tors\\*d|linear torsion|d\\^eta|eta = 9/5|eta=9/5|nonlinear compression|damping/compression",
        "positive_beta_source_content": "Z_beta|positive_beta|positive beta|beta renormalization|beta source|target-independent.*beta",
        "canonical_unit_content": "canonical length|UV unit|unit source|normalization orbit|beta=1|gauge-fixed working normalization",
        "forbidden_replay_content": "selector replay|tau_src -> pair12|beta_tors -> chi_11|chi11|QW-2191|role transfer|ToE closure",
    }
    return {
        "tool": "rg",
        "mode": "content-first legacy->strict bridge-source inventory excluding selector/tau_src/beta_tors->chi11 replay",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def upstream_consistency() -> dict[str, Any]:
    p2679 = load_json(P2679)
    decision = p2679.get("closure_decision", {})
    return {
        "p2679_blocks_old_lanes": decision.get("old_lanes_reopened_now") is False,
        "p2679_no_new_reopening_evidence": decision.get("new_reopening_evidence_after_p2678") is False,
        "p2679_allowed_bridge_source_pivot": "legacy->strict bridge source audit excluding selector/tau_src/beta_tors->chi11 repeats" in decision.get("allowed_next_lanes", []),
        "p2679_hash": sha256_file(P2679),
    }


def source_atom_inventory() -> list[dict[str, Any]]:
    return [
        {
            "atom": "legacy_kernel_source_layer_visible",
            "bridge_component": "input discipline",
            "formal_material_present": True,
            "source_theorem_exported": True,
            "role_safe_for_completion": True,
            "notes": "K_legacy_ont and canonical D_f/alpha_geo/beta_tors layer remain visible as bridge input, not final strict identity.",
        },
        {
            "atom": "strict_kernel_target_layer_visible",
            "bridge_component": "target discipline",
            "formal_material_present": True,
            "source_theorem_exported": True,
            "role_safe_for_completion": True,
            "notes": "K_strict_gate is available as operational strict target, not silently identical to the legacy kernel.",
        },
        {
            "atom": "alpha_geo_scalar_shape_normalization",
            "bridge_component": "amplitude/normalization",
            "formal_material_present": True,
            "source_theorem_exported": True,
            "role_safe_for_completion": False,
            "notes": "Scalar shape normalization exists, but role-safe amplitude absorption/physical-role transfer is not exported.",
        },
        {
            "atom": "amplitude_role_safe_source",
            "bridge_component": "amplitude/normalization",
            "formal_material_present": False,
            "source_theorem_exported": False,
            "role_safe_for_completion": False,
            "notes": "No theorem says legacy alpha_geo roles survive or transform into strict Lagrangian/value roles.",
        },
        {
            "atom": "fractal_pushforward_linear_to_power_damping_shape",
            "bridge_component": "damping/compression",
            "formal_material_present": True,
            "source_theorem_exported": True,
            "role_safe_for_completion": False,
            "notes": "The formal q(d)=d^(9/5) shape route can change linear damping form into power damping shape, but it is not the full coefficient source.",
        },
        {
            "atom": "target_independent_positive_beta_or_z_beta_source",
            "bridge_component": "damping/compression",
            "formal_material_present": False,
            "source_theorem_exported": False,
            "role_safe_for_completion": False,
            "notes": "The missing source atom is a target-independent positive beta/Z_beta coefficient source, not a normalization fit.",
        },
        {
            "atom": "canonical_length_or_uv_unit_source",
            "bridge_component": "damping/compression",
            "formal_material_present": False,
            "source_theorem_exported": False,
            "role_safe_for_completion": False,
            "notes": "The beta=1 normalization orbit is understood, but no canonical length/UV unit theorem fixes the gauge physically.",
        },
        {
            "atom": "selector_phase_orientation_source",
            "bridge_component": "phase/topological selector",
            "formal_material_present": False,
            "source_theorem_exported": False,
            "role_safe_for_completion": False,
            "notes": "Deliberately not reopened by P2680; P2679 forbids selector replay without a new object.",
        },
    ]


def bridge_completion_lattice(atoms: list[dict[str, Any]]) -> dict[str, Any]:
    obligations = [
        "legacy_input_layer_visible",
        "strict_target_layer_visible",
        "amplitude_shape_normalization_present",
        "amplitude_role_safe_source_exported",
        "damping_power_shape_map_present",
        "positive_beta_z_beta_source_exported",
        "canonical_length_uv_unit_exported",
        "selector_replay_excluded_and_guards_preserved",
    ]
    current = {
        obligations[0]: True,
        obligations[1]: True,
        obligations[2]: True,
        obligations[3]: False,
        obligations[4]: True,
        obligations[5]: False,
        obligations[6]: False,
        obligations[7]: True,
    }
    pass_count = 0
    rows = []
    for bits in itertools.product([False, True], repeat=len(obligations)):
        state = dict(zip(obligations, bits))
        passes = all(state.values())
        pass_count += int(passes)
        if passes or state == current:
            rows.append({"state": state, "passes_bridge_source_gate": passes})
    missing = [key for key, value in current.items() if not value]
    return {
        "obligations": obligations,
        "total_states": 2 ** len(obligations),
        "passing_states": pass_count,
        "current_state": current,
        "distinguished_rows": rows,
        "missing_current_obligations": missing,
        "hamming_distance_to_pass": len(missing),
        "full_bridge_completed_now": False,
    }


def component_readiness_matrix(atoms: list[dict[str, Any]]) -> list[dict[str, Any]]:
    components = ["amplitude/normalization", "damping/compression", "phase/topological selector", "role transfer"]
    rows = []
    for component in components:
        related = [atom for atom in atoms if atom["bridge_component"] == component]
        if component == "role transfer":
            rows.append({
                "component": component,
                "ready_now": False,
                "reason": "Role transfer is downstream of full bridge completion and remains forbidden by guardrail.",
                "missing_atoms": ["full bridge completion", "claim-by-claim role-transfer theorem"],
            })
            continue
        ready = bool(related) and all(atom["source_theorem_exported"] and atom["role_safe_for_completion"] for atom in related)
        rows.append({
            "component": component,
            "ready_now": ready,
            "reason": "all source atoms role-safe" if ready else "one or more source atoms remain missing or not role-safe",
            "missing_atoms": [atom["atom"] for atom in related if not (atom["source_theorem_exported"] and atom["role_safe_for_completion"])],
        })
    return rows


def closure_decision(lattice: dict[str, Any], readiness: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "decision": "P2680_LEGACY_STRICT_BRIDGE_SOURCE_INVENTORY__AMPLITUDE_AND_DAMPING_SOURCES_STILL_INCOMPLETE_NO_SELECTOR_REPLAY",
        "professorial_verdict": (
            "P2680 follows the P2679 pivot and does not reopen selector, tau_src->pair12, or beta_tors->chi11.  The bridge-source inventory finds real formal material: legacy/strict kernel layers are visible, scalar alpha_geo shape normalization exists, and a fractal pushforward can supply the damping power-shape route.  The proof-grade bridge still fails because three non-selector source atoms are missing or not role-safe: amplitude role-safe source, target-independent positive beta/Z_beta source, and canonical length/UV unit source."
        ),
        "next_honest_step": (
            "The next honest proof-grade move is not selector work.  Choose one missing non-selector atom and run a construction-or-no-go audit, with highest leverage on the target-independent positive beta/Z_beta source or the canonical length/UV unit source.  Only after those source atoms are exported should role-transfer auditing be attempted."
        ),
        "ready_components": [row["component"] for row in readiness if row["ready_now"]],
        "blocked_components": [row["component"] for row in readiness if not row["ready_now"]],
        "hamming_distance_to_bridge_source_pass": lattice["hamming_distance_to_pass"],
        "full_bridge_completed_now": False,
        "selector_replay_used": False,
        "beta_tors_chi11_reopened_now": False,
        "role_transfer_allowed_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2680/S1630 legacy->strict bridge-source inventory without selector replay",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first repo grep",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    lines.extend(["", "## Source atom inventory"])
    for atom in payload["source_atom_inventory"]:
        lines.append(f"- `{atom['atom']}` ({atom['bridge_component']}): formal=`{atom['formal_material_present']}`, source_export=`{atom['source_theorem_exported']}`, role_safe=`{atom['role_safe_for_completion']}` — {atom['notes']}")
    lat = payload["bridge_completion_lattice"]
    lines.extend([
        "", "## Bridge-source lattice",
        f"Total states: `{lat['total_states']}`; passing states: `{lat['passing_states']}`.",
        f"Current Hamming distance to bridge-source pass: `{lat['hamming_distance_to_pass']}`.",
        f"Missing current obligations: `{lat['missing_current_obligations']}`.",
        "", "## Component readiness",
    ])
    for row in payload["component_readiness_matrix"]:
        lines.append(f"- `{row['component']}`: ready_now=`{row['ready_now']}`; missing={row['missing_atoms']} — {row['reason']}")
    lines.extend([
        "", "## Verdict", payload["closure_decision"]["professorial_verdict"],
        f"Decision: `{payload['closure_decision']['decision']}`.",
        "", "## Next honest step", payload["closure_decision"]["next_honest_step"],
        "", "## Negative exports",
    ])
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    atoms = source_atom_inventory()
    lattice = bridge_completion_lattice(atoms)
    readiness = component_readiness_matrix(atoms)
    payload: dict[str, Any] = {
        "status": "P2680_LEGACY_STRICT_BRIDGE_SOURCE_INVENTORY_NO_SELECTOR_REPLAY_AUDIT_NO_FALSE_PASS",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "upstream_consistency": upstream_consistency(),
        "source_atom_inventory": atoms,
        "bridge_completion_lattice": lattice,
        "component_readiness_matrix": readiness,
        "closure_decision": closure_decision(lattice, readiness),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2680/S1630 legacy-strict bridge-source inventory guard",
        "## P2680/S1630 legacy-strict bridge-source inventory guard\n\n"
        "`P2680/S1630` follows the P2679 pivot without reopening selector, `tau_src -> pair12`, or `beta_tors -> chi_11`.  The non-selector bridge-source inventory finds real formal material: the legacy/strict kernel layers are visible, scalar `alpha_geo` shape normalization exists, and a fractal pushforward can supply the linear-to-power damping shape route.  The full bridge still fails because no role-safe amplitude source, no target-independent positive `beta`/`Z_beta` source, and no canonical length/UV unit source are exported.  Therefore bridge completion, role transfer, O4/O5, `QW-2191` discharge, role-bearing `L_total`, and ToE closure remain blocked.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2680/S1630 bridge-source inventory Ltotal guard",
        "## P2680/S1630 bridge-source inventory Ltotal guard\n\n"
        "`P2680/S1630` keeps `L_total` closed while auditing non-selector bridge-source atoms.  Formal amplitude normalization and damping-shape maps do not become role-bearing variational terms until the missing role-safe amplitude source, positive `beta`/`Z_beta` source, and canonical length/UV unit source are exported and followed by a role-transfer theorem.\n",
    )
    return payload


if __name__ == "__main__":
    main()
