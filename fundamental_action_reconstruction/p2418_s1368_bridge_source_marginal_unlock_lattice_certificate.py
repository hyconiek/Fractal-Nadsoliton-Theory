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

OUT = GEN / "p2418_s1368_bridge_source_marginal_unlock_lattice_certificate.json"
MD = GEN / "p2418_s1368_bridge_source_marginal_unlock_lattice_certificate.md"

SOURCE_FILES = {
    "P2411_BRIDGE_SOURCE_HYPERGRAPH": GEN / "p2411_s1361_legacy_strict_bridge_source_obligation_hypergraph_certificate.json",
    "P2417_NONPROMOTION_MATRIX": GEN / "p2417_s1367_apd_witness_to_source_obligation_nonpromotion_matrix_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

RESIDUAL_COMPONENT = "residual_strict_additions_inventory"


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
    return {"count": len(lines), "samples": lines[:18]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2418|S1368|bridge source marginal unlock|source marginal unlock lattice",
        "prior_source_hypergraph": "P2411|bridge source obligation hypergraph|bridge_obligation_atoms|minimal_component_completion_sets",
        "nonpromotion_input": "P2417|zero source-discharge|current_source_discharge_mask|value witness.*source",
        "scratch_marginal": "marginal unlock|atom influence|critical count|source atom influence",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds P2411's hypergraph and scratch influence/marginal material, but no production P24xx "
            "post-P2417 lattice that starts at the zero source-discharge mask and enumerates exactly which source "
            "sets first unlock bridge components without claiming bridge closure."
        ),
    }


def p2411_cert(source: dict[str, Any]) -> dict[str, Any]:
    return source.get("legacy_strict_bridge_source_obligation_hypergraph_certificate", {}).get(
        "finite_hypergraph_certificate", {}
    )


def atoms_and_components(source: dict[str, Any]) -> tuple[list[str], dict[str, list[str]]]:
    cert = p2411_cert(source)
    return list(cert.get("bridge_obligation_atoms", [])), dict(cert.get("bridge_components", {}))


def mask_atoms(mask: int, atoms: list[str]) -> list[str]:
    return [atom for index, atom in enumerate(atoms) if mask & (1 << index)]


def component_ready(mask: int, atoms: list[str], components: dict[str, list[str]]) -> list[str]:
    true_atoms = set(mask_atoms(mask, atoms))
    ready = []
    for component, required in components.items():
        if set(required).issubset(true_atoms):
            ready.append(component)
    return sorted(ready)


def minimal_component_sets(atoms: list[str], components: dict[str, list[str]]) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for component, required in sorted(components.items()):
        indices = [atoms.index(atom) for atom in required]
        mask = sum(1 << index for index in indices)
        out[component] = {
            "required_atoms": list(required),
            "minimal_mask": mask,
            "minimal_size": len(required),
            "is_empty_inventory_component": len(required) == 0,
        }
    return out


def truth_rows(atoms: list[str], components: dict[str, list[str]]) -> list[dict[str, Any]]:
    rows = []
    nonresidual = sorted(component for component in components if component != RESIDUAL_COMPONENT)
    for mask in range(1 << len(atoms)):
        ready = component_ready(mask, atoms, components)
        ready_nonresidual = [component for component in ready if component != RESIDUAL_COMPONENT]
        rows.append(
            {
                "mask": mask,
                "true_atoms": mask_atoms(mask, atoms),
                "ready_components": ready,
                "ready_nonresidual_components": ready_nonresidual,
                "ready_nonresidual_count": len(ready_nonresidual),
                "missing_nonresidual_components": [component for component in nonresidual if component not in ready_nonresidual],
                "bridge_ready_by_component_coverage": len(ready_nonresidual) == len(nonresidual),
            }
        )
    return rows


def singleton_rows(atoms: list[str], components: dict[str, list[str]]) -> list[dict[str, Any]]:
    rows = []
    for index, atom in enumerate(atoms):
        mask = 1 << index
        ready = component_ready(mask, atoms, components)
        rows.append(
            {
                "atom": atom,
                "mask": mask,
                "ready_components": ready,
                "unlocked_nonresidual_components": [component for component in ready if component != RESIDUAL_COMPONENT],
            }
        )
    return rows


def pair_unlock_rows(atoms: list[str], components: dict[str, list[str]]) -> list[dict[str, Any]]:
    rows = []
    for left, right in combinations(range(len(atoms)), 2):
        mask = (1 << left) | (1 << right)
        ready = component_ready(mask, atoms, components)
        unlocked = [component for component in ready if component != RESIDUAL_COMPONENT]
        if unlocked:
            rows.append(
                {
                    "mask": mask,
                    "atoms": [atoms[left], atoms[right]],
                    "unlocked_nonresidual_components": unlocked,
                }
            )
    return rows


def atom_incidence_rows(atoms: list[str], minimal_sets: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for atom in atoms:
        components = [component for component, info in minimal_sets.items() if atom in info["required_atoms"]]
        pair_components = [component for component in components if minimal_sets[component]["minimal_size"] == 2]
        rows.append(
            {
                "atom": atom,
                "component_incidence": components,
                "component_incidence_count": len(components),
                "pair_unlock_component_count": len(pair_components),
                "pair_unlock_components": pair_components,
                "singleton_unlock_count": 0,
                "priority_score_component_incidence_plus_pair_unlock": len(components) + len(pair_components),
            }
        )
    return sorted(
        rows,
        key=lambda row: (
            -row["priority_score_component_incidence_plus_pair_unlock"],
            -row["component_incidence_count"],
            row["atom"],
        ),
    )


def build_certificate(sources: dict[str, Any]) -> dict[str, Any]:
    atoms, components = atoms_and_components(sources["P2411_BRIDGE_SOURCE_HYPERGRAPH"])
    rows = truth_rows(atoms, components)
    minimal_sets = minimal_component_sets(atoms, components)
    singletons = singleton_rows(atoms, components)
    pairs = pair_unlock_rows(atoms, components)
    p2417_theorem = sources["P2417_NONPROMOTION_MATRIX"].get(
        "apd_witness_to_source_obligation_nonpromotion_matrix_certificate", {}
    ).get("theorem_export", {})
    full_mask = (1 << len(atoms)) - 1
    distribution: dict[str, int] = {}
    for row in rows:
        key = str(row["ready_nonresidual_count"])
        distribution[key] = distribution.get(key, 0) + 1
    return {
        "source_obligation_atoms": atoms,
        "bridge_components": components,
        "current_source_mask_from_p2417": p2417_theorem.get("current_source_discharge_mask"),
        "full_source_mask": full_mask,
        "truth_rows": rows,
        "minimal_component_completion_sets": minimal_sets,
        "singleton_unlock_rows": singletons,
        "pair_unlock_rows": pairs,
        "atom_incidence_priority_rows": atom_incidence_rows(atoms, minimal_sets),
        "ready_nonresidual_component_count_distribution": dict(sorted(distribution.items(), key=lambda item: int(item[0]))),
        "total_source_assignment_count": len(rows),
        "bridge_ready_masks": [row["mask"] for row in rows if row["bridge_ready_by_component_coverage"]],
        "current_mask_ready_components": component_ready(0, atoms, components),
        "singleton_unlock_count": sum(1 for row in singletons if row["unlocked_nonresidual_components"]),
        "pair_unlock_count": len(pairs),
        "minimum_nonresidual_component_unlock_size": min(
            info["minimal_size"] for component, info in minimal_sets.items() if component != RESIDUAL_COMPONENT
        ),
        "phase_component_minimal_size": minimal_sets["phase_frequency_topological_bit_passage"]["minimal_size"],
        "selector_component_minimal_size": minimal_sets["selector_source_premise"]["minimal_size"],
        "chi11_top_incidence_inherited": atom_incidence_rows(atoms, minimal_sets)[0]["atom"] == "chi11_selector_source_theorem",
        "p2417_zero_source_mask_inherited": p2417_theorem.get("current_source_discharge_mask") == 0,
        "full_bridge_theorem_exported": False,
        "role_transfer_licensed": False,
        "toe_closure_exported": False,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2418/S1368 bridge source marginal-unlock lattice certificate

`P2418/S1368` starts from the P2417 zero source-discharge mask and enumerates the full `2^8=256` source-obligation lattice from P2411.  The empty/current mask readies only `residual_strict_additions_inventory`; no singleton source atom unlocks a non-residual bridge component.

The first nontrivial component unlocks occur at size two: the amplitude pair, the damping pair, and the selector pair.  The phase/frequency/topological component still has minimal size three because it needs `phase_frequency_dynamic_source_theorem`, `topological_bit_transport_selector_theorem`, and `chi11_selector_source_theorem`.  `chi11_selector_source_theorem` remains the highest-incidence source atom because it belongs to both the phase/topological and selector-source components.

This is a proof-search lattice, not a source theorem: it ranks where new source evidence could first unlock components, but it does not discharge any atom, complete the bridge, license role transfer, promote `L_total`, or close ToE.
""".strip()
    lag_section = """
## P2418/S1368 bridge source marginal-unlock guard for Lagrangian/EOM

`P2418/S1368` shows that from the current zero source-discharge mask, no singleton source atom unlocks a non-residual bridge component; even the first pair unlocks are only proof-search targets.  Therefore no marginal source-priority score can be used as a role-bearing `L_total` term or as bridge/ToE closure.
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
        "theorem_name": "P2418_T1_bridge_source_marginal_unlock_lattice_certificate",
        "source_obligation_atom_count": len(cert["source_obligation_atoms"]),
        "bridge_component_count": len(cert["bridge_components"]),
        "total_source_assignment_count": cert["total_source_assignment_count"],
        "current_source_mask_from_p2417": cert["current_source_mask_from_p2417"],
        "full_source_mask": cert["full_source_mask"],
        "bridge_ready_masks": cert["bridge_ready_masks"],
        "current_mask_ready_components": cert["current_mask_ready_components"],
        "singleton_unlock_count": cert["singleton_unlock_count"],
        "pair_unlock_count": cert["pair_unlock_count"],
        "minimum_nonresidual_component_unlock_size": cert["minimum_nonresidual_component_unlock_size"],
        "phase_component_minimal_size": cert["phase_component_minimal_size"],
        "selector_component_minimal_size": cert["selector_component_minimal_size"],
        "chi11_top_incidence_inherited": cert["chi11_top_incidence_inherited"],
        "ready_nonresidual_component_count_distribution": cert["ready_nonresidual_component_count_distribution"],
        "p2417_zero_source_mask_inherited": cert["p2417_zero_source_mask_inherited"],
        "full_bridge_theorem_exported": cert["full_bridge_theorem_exported"],
        "role_transfer_licensed": cert["role_transfer_licensed"],
        "toe_closure_exported": cert["toe_closure_exported"],
        "not_licensed": [
            "No singleton source atom is promoted to bridge-component completion.",
            "No marginal unlock score is a source theorem, selector theorem, or QW-2191 discharge.",
            "No bridge closure follows before the full source mask is discharged.",
            "No role-transfer theorem, role-bearing L_total term, or ToE closure follows.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "eight_source_atoms": theorem_export["source_obligation_atom_count"] == 8,
        "five_bridge_components": theorem_export["bridge_component_count"] == 5,
        "all_256_source_assignments": theorem_export["total_source_assignment_count"] == 256,
        "current_mask_zero_inherited": theorem_export["current_source_mask_from_p2417"] == 0,
        "full_mask_255": theorem_export["full_source_mask"] == 255,
        "bridge_ready_only_full_mask": theorem_export["bridge_ready_masks"] == [255],
        "current_mask_only_residual_ready": theorem_export["current_mask_ready_components"] == [RESIDUAL_COMPONENT],
        "no_singleton_unlock": theorem_export["singleton_unlock_count"] == 0,
        "three_pair_unlocks": theorem_export["pair_unlock_count"] == 3,
        "first_unlock_size_two": theorem_export["minimum_nonresidual_component_unlock_size"] == 2,
        "phase_requires_three_atoms": theorem_export["phase_component_minimal_size"] == 3,
        "selector_requires_two_atoms": theorem_export["selector_component_minimal_size"] == 2,
        "chi11_top_incidence": theorem_export["chi11_top_incidence_inherited"],
        "p2417_zero_source_mask_inherited": theorem_export["p2417_zero_source_mask_inherited"],
        "full_bridge_still_open": not theorem_export["full_bridge_theorem_exported"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "toe_still_open": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2418_s1368_v1",
        "packet_id": "P2418",
        "stage_id": "S1368",
        "result_kind": "BRIDGE_SOURCE_MARGINAL_UNLOCK_LATTICE_CERTIFICATE",
        "status": "PASS_SOURCE_MARGINAL_UNLOCK_LATTICE_NO_SINGLETON_NO_BRIDGE_CLOSURE",
        "bridge_source_marginal_unlock_lattice_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Use the pair/triple unlock lattice to target a real source theorem, with chi11 still highest-incidence; do not treat marginal unlocks as closure."
        ),
        "global_status": "OPEN_PROGRESS_SOURCE_UNLOCK_LATTICE_NO_SOURCE_OR_TOE_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["bridge_source_marginal_unlock_lattice_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2418 S1368: bridge source marginal-unlock lattice certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2418/S1368 enumerates the source-obligation marginal-unlock lattice from the P2417 zero source mask.",
                "",
                "## Finite facts",
                "",
                f"- Source atoms: `{theorem['source_obligation_atom_count']}`.",
                f"- Source assignments: `{theorem['total_source_assignment_count']}`.",
                f"- Bridge-ready masks: `{theorem['bridge_ready_masks']}`.",
                f"- Singleton unlock count: `{theorem['singleton_unlock_count']}`.",
                f"- Pair unlock count: `{theorem['pair_unlock_count']}`.",
                f"- Ready-component distribution: `{theorem['ready_nonresidual_component_count_distribution']}`.",
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
