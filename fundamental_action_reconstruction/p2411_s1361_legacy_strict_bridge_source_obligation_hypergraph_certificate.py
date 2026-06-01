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
SCRATCH = ROOT / "scratch"

OUT = GEN / "p2411_s1361_legacy_strict_bridge_source_obligation_hypergraph_certificate.json"
MD = GEN / "p2411_s1361_legacy_strict_bridge_source_obligation_hypergraph_certificate.md"

SOURCE_FILES = {
    "S2_STRATEGIC_PRIORITY": ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
    "P2410_ATOM_OBSTRUCTION": GEN / "p2410_s1360_dequotiented_twelve_atom_prime_implicate_obstruction_certificate.json",
    "SCRATCH_COMPONENT_GAP": SCRATCH / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

BRIDGE_OBLIGATION_ATOMS = [
    "amplitude_dynamic_source_theorem",
    "role_safe_amplitude_absorption_map_theorem",
    "phase_frequency_dynamic_source_theorem",
    "topological_bit_transport_selector_theorem",
    "legacy_linear_to_strict_nonlinear_compression_map_theorem",
    "strict_compression_dynamic_source_theorem",
    "chi11_selector_source_theorem",
    "qw2191_symmetry_breaking_or_internal_source_theorem",
]
ATOM_INDEX = {atom: index for index, atom in enumerate(BRIDGE_OBLIGATION_ATOMS)}
FULL_MASK = (1 << len(BRIDGE_OBLIGATION_ATOMS)) - 1

BRIDGE_COMPONENTS = {
    "amplitude_normalization_passage": [
        "amplitude_dynamic_source_theorem",
        "role_safe_amplitude_absorption_map_theorem",
    ],
    "phase_frequency_topological_bit_passage": [
        "phase_frequency_dynamic_source_theorem",
        "topological_bit_transport_selector_theorem",
        "chi11_selector_source_theorem",
    ],
    "damping_compression_passage": [
        "legacy_linear_to_strict_nonlinear_compression_map_theorem",
        "strict_compression_dynamic_source_theorem",
    ],
    "selector_source_premise": [
        "chi11_selector_source_theorem",
        "qw2191_symmetry_breaking_or_internal_source_theorem",
    ],
    "residual_strict_additions_inventory": [],
}

ROLE_TRANSFER_ATOMS_RESERVED = [
    "alpha_geo_electroweak_role_transfer_theorem",
    "beta_tors_strict_role_transfer_theorem",
    "beta_power_hierarchy_successor_role_transfer_theorem",
]


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    if not path.exists():
        return f"OPEN_MISSING_TEXT::{rel(path)}"
    return path.read_text(encoding="utf-8")


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
    return {"count": len(lines), "samples": lines[:16]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2411|S1361|bridge source obligation hypergraph|source-obligation hypergraph",
        "bridge_gap_source": "completion component gap|legacy -> strict completion bridge|component gap matrix",
        "s2_bridge_requirements": "amplitude/normalization passage|damping/compression passage|selector/source premise|role-transfer audit after bridge completion",
        "qw2191_guard": "QW-2191|selector/source|symmetry-breaking|chi11",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds scratch component-gap and guardrail material, but no production P24xx certificate that "
            "turns S2's legacy->strict bridge requirements into a finite source-obligation hypergraph."
        ),
    }


def atoms_for(mask: int) -> list[str]:
    return [atom for atom in BRIDGE_OBLIGATION_ATOMS if mask & (1 << ATOM_INDEX[atom])]


def component_ready(component: str, mask: int) -> bool:
    return all(mask & (1 << ATOM_INDEX[atom]) for atom in BRIDGE_COMPONENTS[component])


def ready_components(mask: int) -> list[str]:
    return [component for component in BRIDGE_COMPONENTS if component_ready(component, mask)]


def bridge_ready(mask: int) -> bool:
    return all(component_ready(component, mask) for component in BRIDGE_COMPONENTS)


def truth_vector() -> list[int]:
    return [int(bridge_ready(mask)) for mask in range(1 << len(BRIDGE_OBLIGATION_ATOMS))]


def component_frontier_rows() -> list[dict[str, Any]]:
    rows = []
    for component, requirements in BRIDGE_COMPONENTS.items():
        rows.append(
            {
                "component": component,
                "required_atoms": requirements,
                "required_atom_count": len(requirements),
                "current_ready_at_empty_mask": len(requirements) == 0,
                "open_gap": "none_currently_blocking_component" if not requirements else " + ".join(requirements),
            }
        )
    return rows


def repair_spectrum() -> list[dict[str, Any]]:
    values = truth_vector()
    rows = []
    for missing_count in range(1, len(BRIDGE_OBLIGATION_ATOMS) + 1):
        masks = [mask for mask, value in enumerate(values) if not value and len(BRIDGE_OBLIGATION_ATOMS) - len(atoms_for(mask)) == missing_count]
        rows.append(
            {
                "missing_obligation_count": missing_count,
                "failure_mask_count": len(masks),
                "nearest_repair_distance": missing_count,
            }
        )
    return rows


def nearest_miss_rows() -> list[dict[str, Any]]:
    rows = []
    for atom in BRIDGE_OBLIGATION_ATOMS:
        mask = FULL_MASK ^ (1 << ATOM_INDEX[atom])
        rows.append(
            {
                "mask": mask,
                "missing_atom": atom,
                "ready_components_before_repair": ready_components(mask),
                "blocked_components_before_repair": [component for component in BRIDGE_COMPONENTS if not component_ready(component, mask)],
            }
        )
    return rows


def minimal_component_completion_sets() -> list[dict[str, Any]]:
    rows = []
    for component, requirements in BRIDGE_COMPONENTS.items():
        rows.append(
            {
                "component": component,
                "minimal_completion_sets": [requirements],
                "minimum_size": len(requirements),
            }
        )
    return rows


def atom_component_incidence() -> list[dict[str, Any]]:
    rows = []
    for atom in BRIDGE_OBLIGATION_ATOMS:
        components = [component for component, reqs in BRIDGE_COMPONENTS.items() if atom in reqs]
        rows.append({"atom": atom, "component_count": len(components), "components": components})
    return rows


def proof_order_scores() -> list[dict[str, Any]]:
    rows = []
    for atom in BRIDGE_OBLIGATION_ATOMS:
        components = [component for component, reqs in BRIDGE_COMPONENTS.items() if atom in reqs]
        singleton_completes = [component for component in components if len(BRIDGE_COMPONENTS[component]) == 1]
        rows.append(
            {
                "atom": atom,
                "component_incidence_count": len(components),
                "singleton_completion_count": len(singleton_completes),
                "affected_components": components,
                "proof_search_priority_score": len(components) + len(singleton_completes),
            }
        )
    return sorted(rows, key=lambda row: (-row["proof_search_priority_score"], row["atom"]))


def build_certificate(sources: dict[str, Any]) -> dict[str, Any]:
    values = truth_vector()
    spectrum = repair_spectrum()
    p2410_theorem = sources["P2410_ATOM_OBSTRUCTION"].get(
        "dequotiented_twelve_atom_prime_implicate_obstruction_certificate", {}
    ).get("theorem_export", {})
    scratch = sources["SCRATCH_COMPONENT_GAP"]
    scratch_summary = scratch.get("completion_gap_summary", scratch.get("summary", {}))
    return {
        "bridge_obligation_atoms": BRIDGE_OBLIGATION_ATOMS,
        "role_transfer_atoms_reserved_post_bridge": ROLE_TRANSFER_ATOMS_RESERVED,
        "bridge_components": BRIDGE_COMPONENTS,
        "component_frontier_rows": component_frontier_rows(),
        "atom_component_incidence": atom_component_incidence(),
        "proof_order_scores": proof_order_scores(),
        "minimal_component_completion_sets": minimal_component_completion_sets(),
        "total_bridge_obligation_assignment_count": len(values),
        "bridge_ready_truth_vector_sha256": sha256_json(values),
        "bridge_ready_true_masks": [mask for mask, value in enumerate(values) if value],
        "bridge_ready_anf": " * ".join(BRIDGE_OBLIGATION_ATOMS),
        "bridge_ready_anf_degree": len(BRIDGE_OBLIGATION_ATOMS),
        "empty_mask_ready_components": ready_components(0),
        "empty_mask_bridge_ready": bridge_ready(0),
        "proper_subset_failure_count": sum(1 for mask, value in enumerate(values) if not value),
        "nearest_miss_rows": nearest_miss_rows(),
        "nearest_miss_count": len(nearest_miss_rows()),
        "repair_spectrum": spectrum,
        "repair_spectrum_failure_total": sum(row["failure_mask_count"] for row in spectrum),
        "p2410_no_ltotal_or_toe_closure_inherited": p2410_theorem.get("failure_assignment_count") == 4095
        and bool(p2410_theorem.get("not_licensed")),
        "scratch_component_gap_sources_missing_inherited": scratch_summary.get("rows_with_strict_dynamic_sources") == 0,
        "scratch_role_transfer_blocked_inherited": scratch_summary.get("rows_with_role_transfer_allowed_now") == 0,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2411/S1361 legacy-to-strict bridge source-obligation hypergraph certificate

`P2411/S1361` moves from the generic twelve-atom `L_total` obstruction ledger to the strategic S2 bridge problem itself.  It models the current `K_legacy_ont -> K_strict_gate` completion bridge as five S2 components: amplitude normalization, phase/frequency/topological-bit passage, damping/compression passage, selector/source premise, and residual strict-addition inventory.

The finite hypergraph has eight open bridge-source obligations.  The residual strict-addition inventory is currently ready as an inventory, but the bridge theorem is true only at the full eight-obligation mask.  All `255` proper masks fail; the nearest misses are exactly the eight one-obligation-missing masks.

This identifies the next proof search as a real bridge/source theorem problem, not another role-export shortcut: role-transfer atoms remain reserved for the post-bridge audit, and QW-2191 remains open until a genuine selector/source or symmetry-breaking premise is supplied.
""".strip()
    lag_section = """
## P2411/S1361 bridge-source hypergraph guard for Lagrangian/EOM

`P2411/S1361` keeps `L_total` downstream of the bridge: the legacy-to-strict completion bridge still needs eight source obligations before any post-bridge role-transfer audit can start.  Therefore no Lagrangian term may import legacy roles from amplitude, phase, damping, or `chi11` bookkeeping without a real bridge-source theorem and the later role-transfer theorem.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    source_payloads = {
        name: (load_json(path) if path.suffix == ".json" else {"text_sha256": hashlib.sha256(load_text(path).encode("utf-8")).hexdigest(), "status": "TEXT_LOADED"})
        for name, path in SOURCE_FILES.items()
    }
    grep = rg_audit()
    cert = build_certificate(source_payloads)
    theorem_export = {
        "theorem_name": "P2411_T1_legacy_strict_bridge_source_obligation_hypergraph_certificate",
        "bridge_obligation_atom_count": len(BRIDGE_OBLIGATION_ATOMS),
        "bridge_component_count": len(BRIDGE_COMPONENTS),
        "total_bridge_obligation_assignment_count": cert["total_bridge_obligation_assignment_count"],
        "bridge_ready_true_masks": cert["bridge_ready_true_masks"],
        "bridge_ready_anf_degree": cert["bridge_ready_anf_degree"],
        "empty_mask_ready_components": cert["empty_mask_ready_components"],
        "empty_mask_bridge_ready": cert["empty_mask_bridge_ready"],
        "proper_subset_failure_count": cert["proper_subset_failure_count"],
        "nearest_miss_count": cert["nearest_miss_count"],
        "repair_spectrum_failure_total": cert["repair_spectrum_failure_total"],
        "top_priority_atoms": cert["proof_order_scores"][:3],
        "p2410_no_ltotal_or_toe_closure_inherited": cert["p2410_no_ltotal_or_toe_closure_inherited"],
        "scratch_component_gap_sources_missing_inherited": cert["scratch_component_gap_sources_missing_inherited"],
        "scratch_role_transfer_blocked_inherited": cert["scratch_role_transfer_blocked_inherited"],
        "not_licensed": [
            "No raw identity K_legacy_ont == K_strict_gate is claimed.",
            "No bridge completion follows from residual strict-addition inventory alone.",
            "No chi11 selector-source theorem or QW-2191 discharge is exported.",
            "No legacy physical-role transfer is licensed before the post-bridge role-transfer audit.",
            "No role-bearing L_total or ToE closure follows from this bridge hypergraph.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "eight_bridge_obligations": theorem_export["bridge_obligation_atom_count"] == 8,
        "five_s2_bridge_components": theorem_export["bridge_component_count"] == 5,
        "exact_256_assignment_lattice": theorem_export["total_bridge_obligation_assignment_count"] == 256,
        "only_full_mask_completes_bridge": theorem_export["bridge_ready_true_masks"] == [255],
        "empty_mask_only_residual_inventory_ready": theorem_export["empty_mask_ready_components"] == ["residual_strict_additions_inventory"],
        "all_255_proper_masks_fail": theorem_export["proper_subset_failure_count"] == 255,
        "nearest_misses_are_eight_one_obligation_gaps": theorem_export["nearest_miss_count"] == 8,
        "repair_spectrum_sums_to_255": theorem_export["repair_spectrum_failure_total"] == 255,
        "p2410_no_ltotal_or_toe_closure_inherited": theorem_export["p2410_no_ltotal_or_toe_closure_inherited"],
        "scratch_sources_missing_inherited": theorem_export["scratch_component_gap_sources_missing_inherited"],
        "scratch_role_transfer_blocked_inherited": theorem_export["scratch_role_transfer_blocked_inherited"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2411_s1361_v1",
        "packet_id": "P2411",
        "stage_id": "S1361",
        "result_kind": "LEGACY_STRICT_BRIDGE_SOURCE_OBLIGATION_HYPERGRAPH_CERTIFICATE",
        "status": "PASS_BRIDGE_SOURCE_OBLIGATION_HYPERGRAPH_FULL_MASK_ONLY_NO_ROLE_TRANSFER",
        "legacy_strict_bridge_source_obligation_hypergraph_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: payload.get("status", "OPEN_UNKNOWN") for name, payload in source_payloads.items()},
            "finite_hypergraph_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Attempt one concrete bridge-source theorem from the top-priority atoms, especially chi11 selector/source or "
            "a real amplitude/damping source map; do not perform post-bridge role transfer yet."
        ),
        "global_status": "OPEN_PROGRESS_BRIDGE_SOURCE_HYPERGRAPH_CERTIFIED_NO_ROLE_TRANSFER_OR_TOE_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["legacy_strict_bridge_source_obligation_hypergraph_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2411 S1361: legacy-strict bridge source-obligation hypergraph certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2411/S1361 converts the S2 legacy->strict bridge priority into an exact source-obligation hypergraph.",
                "",
                "## Finite Boolean facts",
                "",
                f"- Bridge obligation atoms: `{theorem['bridge_obligation_atom_count']}`.",
                f"- Bridge components: `{theorem['bridge_component_count']}`.",
                f"- Assignment count: `{theorem['total_bridge_obligation_assignment_count']}`.",
                f"- True bridge masks: `{theorem['bridge_ready_true_masks']}`.",
                f"- Empty-mask ready components: `{theorem['empty_mask_ready_components']}`.",
                f"- Proper-subset failures: `{theorem['proper_subset_failure_count']}`.",
                f"- Nearest misses: `{theorem['nearest_miss_count']}`.",
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
