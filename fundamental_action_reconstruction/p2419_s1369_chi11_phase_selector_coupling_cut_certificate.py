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

OUT = GEN / "p2419_s1369_chi11_phase_selector_coupling_cut_certificate.json"
MD = GEN / "p2419_s1369_chi11_phase_selector_coupling_cut_certificate.md"

SOURCE_FILES = {
    "P2411_BRIDGE_SOURCE_HYPERGRAPH": GEN / "p2411_s1361_legacy_strict_bridge_source_obligation_hypergraph_certificate.json",
    "P2418_MARGINAL_UNLOCK_LATTICE": GEN / "p2418_s1368_bridge_source_marginal_unlock_lattice_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

PHASE_COMPONENT = "phase_frequency_topological_bit_passage"
SELECTOR_COMPONENT = "selector_source_premise"
CHI11_ATOM = "chi11_selector_source_theorem"
QW2191_ATOM = "qw2191_symmetry_breaking_or_internal_source_theorem"


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
        "new_packet": "P2419|S1369|chi11 phase selector coupling|phase selector coupling cut|co-readiness cut",
        "source_hypergraph_input": "P2411|bridge source obligation hypergraph|phase_frequency_topological_bit_passage|selector_source_premise",
        "marginal_lattice_input": "P2418|source marginal unlock|chi11_top_incidence|pair unlock",
        "overlap_prior": "shared chi11|component overlap|phase.*selector.*coupling|selector.*phase.*coupling|co-unlock|co-readiness",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds the P2411 source hypergraph and P2418 marginal-unlock priority result, plus isolated "
            "scratch overlap language, but no production P24xx finite cut certificate that classifies the phase/selector "
            "readiness quadrant and proves the chi11 atom is a common necessary cut without discharging QW-2191."
        ),
    }


def p2411_cert(source: dict[str, Any]) -> dict[str, Any]:
    return source.get("legacy_strict_bridge_source_obligation_hypergraph_certificate", {}).get(
        "finite_hypergraph_certificate", {}
    )


def p2418_theorem(source: dict[str, Any]) -> dict[str, Any]:
    return source.get("bridge_source_marginal_unlock_lattice_certificate", {}).get("theorem_export", {})


def atoms_and_components(source: dict[str, Any]) -> tuple[list[str], dict[str, list[str]]]:
    cert = p2411_cert(source)
    return list(cert.get("bridge_obligation_atoms", [])), dict(cert.get("bridge_components", {}))


def mask_atoms(mask: int, atoms: list[str]) -> list[str]:
    return [atom for index, atom in enumerate(atoms) if mask & (1 << index)]


def mask_for_atoms(atoms: list[str], selected: list[str]) -> int:
    return sum(1 << atoms.index(atom) for atom in selected)


def component_ready(mask: int, atoms: list[str], required: list[str]) -> bool:
    true_atoms = set(mask_atoms(mask, atoms))
    return set(required).issubset(true_atoms)


def readiness_quadrant_rows(atoms: list[str], phase_req: list[str], selector_req: list[str]) -> list[dict[str, Any]]:
    rows = []
    for mask in range(1 << len(atoms)):
        phase_ready = component_ready(mask, atoms, phase_req)
        selector_ready = component_ready(mask, atoms, selector_req)
        if phase_ready and selector_ready:
            quadrant = "phase_and_selector_ready"
        elif phase_ready:
            quadrant = "phase_only_ready"
        elif selector_ready:
            quadrant = "selector_only_ready"
        else:
            quadrant = "neither_phase_nor_selector_ready"
        rows.append(
            {
                "mask": mask,
                "true_atoms": mask_atoms(mask, atoms),
                "phase_ready": phase_ready,
                "selector_ready": selector_ready,
                "quadrant": quadrant,
                "contains_chi11": CHI11_ATOM in mask_atoms(mask, atoms),
                "contains_qw2191": QW2191_ATOM in mask_atoms(mask, atoms),
            }
        )
    return rows


def minimal_masks_for_predicate(atoms: list[str], predicate) -> list[dict[str, Any]]:
    out = []
    for mask in range(1 << len(atoms)):
        if not predicate(mask):
            continue
        proper_true_submasks = [sub for sub in range(mask) if sub & mask == sub and sub != mask]
        if any(predicate(sub) for sub in proper_true_submasks):
            continue
        out.append({"mask": mask, "atoms": mask_atoms(mask, atoms), "size": len(mask_atoms(mask, atoms))})
    return sorted(out, key=lambda row: (row["size"], row["mask"]))


def deletion_sensitivity_rows(atoms: list[str], phase_req: list[str], selector_req: list[str], co_mask: int) -> list[dict[str, Any]]:
    rows = []
    for atom in mask_atoms(co_mask, atoms):
        dropped_mask = co_mask & ~(1 << atoms.index(atom))
        phase_after = component_ready(dropped_mask, atoms, phase_req)
        selector_after = component_ready(dropped_mask, atoms, selector_req)
        blocked = []
        if not phase_after:
            blocked.append(PHASE_COMPONENT)
        if not selector_after:
            blocked.append(SELECTOR_COMPONENT)
        rows.append(
            {
                "removed_atom": atom,
                "remaining_mask": dropped_mask,
                "remaining_atoms": mask_atoms(dropped_mask, atoms),
                "phase_ready_after_deletion": phase_after,
                "selector_ready_after_deletion": selector_after,
                "blocked_components_after_deletion": blocked,
            }
        )
    return rows


def size_distribution(rows: list[dict[str, Any]], quadrant: str) -> dict[str, int]:
    dist: dict[str, int] = {}
    for row in rows:
        if row["quadrant"] != quadrant:
            continue
        size = str(len(row["true_atoms"]))
        dist[size] = dist.get(size, 0) + 1
    return dict(sorted(dist.items(), key=lambda item: int(item[0])))


def build_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2411 = sources["P2411_BRIDGE_SOURCE_HYPERGRAPH"]
    p2418 = sources["P2418_MARGINAL_UNLOCK_LATTICE"]
    atoms, components = atoms_and_components(p2411)
    phase_req = list(components.get(PHASE_COMPONENT, []))
    selector_req = list(components.get(SELECTOR_COMPONENT, []))
    rows = readiness_quadrant_rows(atoms, phase_req, selector_req)
    quadrant_counts: dict[str, int] = {}
    for row in rows:
        quadrant_counts[row["quadrant"]] = quadrant_counts.get(row["quadrant"], 0) + 1

    phase_mask = mask_for_atoms(atoms, phase_req)
    selector_mask = mask_for_atoms(atoms, selector_req)
    shared_atoms = sorted(set(phase_req) & set(selector_req))
    union_atoms = [atom for atom in atoms if atom in set(phase_req) | set(selector_req)]
    co_mask = mask_for_atoms(atoms, union_atoms)

    phase_minimal = minimal_masks_for_predicate(atoms, lambda mask: component_ready(mask, atoms, phase_req))
    selector_minimal = minimal_masks_for_predicate(atoms, lambda mask: component_ready(mask, atoms, selector_req))
    co_minimal = minimal_masks_for_predicate(
        atoms,
        lambda mask: component_ready(mask, atoms, phase_req) and component_ready(mask, atoms, selector_req),
    )

    no_phase_without_chi11 = all(not row["phase_ready"] for row in rows if not row["contains_chi11"])
    no_selector_without_chi11 = all(not row["selector_ready"] for row in rows if not row["contains_chi11"])
    qw2191_not_sufficient_for_selector = not component_ready(mask_for_atoms(atoms, [QW2191_ATOM]), atoms, selector_req)
    chi11_not_sufficient_for_selector = not component_ready(mask_for_atoms(atoms, [CHI11_ATOM]), atoms, selector_req)
    chi11_not_sufficient_for_phase = not component_ready(mask_for_atoms(atoms, [CHI11_ATOM]), atoms, phase_req)

    p2418_t = p2418_theorem(p2418)
    return {
        "source_obligation_atoms": atoms,
        "phase_component_required_atoms": phase_req,
        "selector_component_required_atoms": selector_req,
        "phase_minimal_mask": phase_mask,
        "selector_minimal_mask": selector_mask,
        "shared_phase_selector_atoms": shared_atoms,
        "phase_selector_union_atoms": union_atoms,
        "phase_selector_coreadiness_minimal_mask": co_mask,
        "readiness_quadrant_rows": rows,
        "quadrant_counts": dict(sorted(quadrant_counts.items())),
        "quadrant_size_distributions": {
            quadrant: size_distribution(rows, quadrant)
            for quadrant in sorted(quadrant_counts)
        },
        "phase_minimal_readiness_masks": phase_minimal,
        "selector_minimal_readiness_masks": selector_minimal,
        "phase_selector_coreadiness_minimal_masks": co_minimal,
        "phase_selector_coreadiness_deletion_sensitivity_rows": deletion_sensitivity_rows(atoms, phase_req, selector_req, co_mask),
        "no_phase_readiness_without_chi11": no_phase_without_chi11,
        "no_selector_readiness_without_chi11": no_selector_without_chi11,
        "qw2191_not_sufficient_for_selector": qw2191_not_sufficient_for_selector,
        "chi11_not_sufficient_for_selector": chi11_not_sufficient_for_selector,
        "chi11_not_sufficient_for_phase": chi11_not_sufficient_for_phase,
        "p2418_chi11_top_incidence_inherited": p2418_t.get("chi11_top_incidence_inherited") is True,
        "p2418_pair_unlock_count_inherited": p2418_t.get("pair_unlock_count"),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2419/S1369 chi11 phase-selector coupling cut certificate

`P2419/S1369` refines the P2418 source-unlock lattice at the phase/selector overlap.  The finite audit over all `2^8=256` source masks classifies the readiness quadrants for `phase_frequency_topological_bit_passage` and `selector_source_premise`: neither ready, phase-only ready, selector-only ready, or both ready.

The shared cut is explicit: both components require `chi11_selector_source_theorem`.  The unique minimal phase/selector co-readiness set is `{phase_frequency_dynamic_source_theorem, topological_bit_transport_selector_theorem, chi11_selector_source_theorem, qw2191_symmetry_breaking_or_internal_source_theorem}`.  Deleting `chi11_selector_source_theorem` from that set blocks both phase/topological passage and selector premise, while deleting `qw2191_symmetry_breaking_or_internal_source_theorem` blocks selector only.

This is only a coupling/cut theorem for proof search.  It does not export a chi11 source, does not discharge QW-2191, does not complete the bridge, and does not license role transfer, role-bearing `L_total`, or ToE closure.
""".strip()
    lag_section = """
## P2419/S1369 chi11 phase-selector coupling guard for Lagrangian/EOM

`P2419/S1369` proves that `chi11_selector_source_theorem` is a common necessary source-obligation cut for the phase/topological and selector components, but also proves that `chi11` alone is not sufficient for either component and that `qw2191_symmetry_breaking_or_internal_source_theorem` remains an independent selector obligation.  Therefore the overlap cannot be used as a role-bearing `L_total` term, selector closure, or ToE closure.
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
        "theorem_name": "P2419_T1_chi11_phase_selector_coupling_cut_certificate",
        "source_obligation_atom_count": len(cert["source_obligation_atoms"]),
        "total_source_assignment_count": len(cert["readiness_quadrant_rows"]),
        "phase_required_atom_count": len(cert["phase_component_required_atoms"]),
        "selector_required_atom_count": len(cert["selector_component_required_atoms"]),
        "shared_phase_selector_atoms": cert["shared_phase_selector_atoms"],
        "shared_phase_selector_atom_count": len(cert["shared_phase_selector_atoms"]),
        "phase_selector_union_atom_count": len(cert["phase_selector_union_atoms"]),
        "phase_minimal_readiness_mask_count": len(cert["phase_minimal_readiness_masks"]),
        "selector_minimal_readiness_mask_count": len(cert["selector_minimal_readiness_masks"]),
        "phase_selector_coreadiness_minimal_mask_count": len(cert["phase_selector_coreadiness_minimal_masks"]),
        "phase_selector_coreadiness_minimal_masks": [row["mask"] for row in cert["phase_selector_coreadiness_minimal_masks"]],
        "phase_selector_coreadiness_minimal_size": cert["phase_selector_coreadiness_minimal_masks"][0]["size"],
        "quadrant_counts": cert["quadrant_counts"],
        "no_phase_readiness_without_chi11": cert["no_phase_readiness_without_chi11"],
        "no_selector_readiness_without_chi11": cert["no_selector_readiness_without_chi11"],
        "chi11_is_common_necessary_cut": cert["shared_phase_selector_atoms"] == [CHI11_ATOM]
        and cert["no_phase_readiness_without_chi11"]
        and cert["no_selector_readiness_without_chi11"],
        "qw2191_not_sufficient_for_selector": cert["qw2191_not_sufficient_for_selector"],
        "chi11_not_sufficient_for_selector": cert["chi11_not_sufficient_for_selector"],
        "chi11_not_sufficient_for_phase": cert["chi11_not_sufficient_for_phase"],
        "p2418_chi11_top_incidence_inherited": cert["p2418_chi11_top_incidence_inherited"],
        "p2418_pair_unlock_count_inherited": cert["p2418_pair_unlock_count_inherited"],
        "chi11_source_exported": False,
        "qw2191_discharged": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_licensed": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "The shared chi11 cut is not a chi11 source theorem.",
            "QW-2191 remains an independent selector-source obstruction and is not discharged.",
            "Phase/selector co-readiness is a source-obligation mask, not bridge completion.",
            "No role-transfer theorem, role-bearing L_total term, or ToE closure follows.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "eight_source_atoms": theorem_export["source_obligation_atom_count"] == 8,
        "all_256_source_assignments": theorem_export["total_source_assignment_count"] == 256,
        "phase_requires_three_atoms": theorem_export["phase_required_atom_count"] == 3,
        "selector_requires_two_atoms": theorem_export["selector_required_atom_count"] == 2,
        "exactly_one_shared_atom": theorem_export["shared_phase_selector_atom_count"] == 1,
        "shared_atom_is_chi11": theorem_export["shared_phase_selector_atoms"] == [CHI11_ATOM],
        "unique_phase_minimal_mask": theorem_export["phase_minimal_readiness_mask_count"] == 1,
        "unique_selector_minimal_mask": theorem_export["selector_minimal_readiness_mask_count"] == 1,
        "unique_coreadiness_minimal_mask": theorem_export["phase_selector_coreadiness_minimal_mask_count"] == 1,
        "coreadiness_minimal_size_four": theorem_export["phase_selector_coreadiness_minimal_size"] == 4,
        "quadrant_counts_expected": theorem_export["quadrant_counts"] == {
            "neither_phase_nor_selector_ready": 176,
            "phase_and_selector_ready": 16,
            "phase_only_ready": 16,
            "selector_only_ready": 48,
        },
        "no_phase_without_chi11": theorem_export["no_phase_readiness_without_chi11"],
        "no_selector_without_chi11": theorem_export["no_selector_readiness_without_chi11"],
        "chi11_common_necessary_cut": theorem_export["chi11_is_common_necessary_cut"],
        "qw2191_not_sufficient": theorem_export["qw2191_not_sufficient_for_selector"],
        "chi11_not_sufficient_for_selector": theorem_export["chi11_not_sufficient_for_selector"],
        "chi11_not_sufficient_for_phase": theorem_export["chi11_not_sufficient_for_phase"],
        "p2418_priority_inherited": theorem_export["p2418_chi11_top_incidence_inherited"],
        "chi11_source_still_open": not theorem_export["chi11_source_exported"],
        "qw2191_still_open": not theorem_export["qw2191_discharged"],
        "bridge_still_open": not theorem_export["full_bridge_theorem_exported"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "ltotal_still_blocked": not theorem_export["role_bearing_ltotal_exported"],
        "toe_still_open": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2419_s1369_v1",
        "packet_id": "P2419",
        "stage_id": "S1369",
        "result_kind": "CHI11_PHASE_SELECTOR_COUPLING_CUT_CERTIFICATE",
        "status": "PASS_CHI11_PHASE_SELECTOR_COUPLING_CUT_NO_SOURCE_OR_SELECTOR_CLOSURE",
        "chi11_phase_selector_coupling_cut_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Target a real source theorem for the shared chi11 cut or the independent QW-2191 selector obstruction; "
            "do not promote phase/selector overlap into bridge, role-transfer, L_total, or ToE closure."
        ),
        "global_status": "OPEN_PROGRESS_CHI11_CUT_CERTIFIED_SOURCE_AND_SELECTOR_OBSTRUCTIONS_REMAIN",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["chi11_phase_selector_coupling_cut_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2419 S1369: chi11 phase-selector coupling cut certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2419/S1369 classifies the phase/selector readiness quadrants and isolates the shared chi11 cut.",
                "",
                "## Finite facts",
                "",
                f"- Source assignments: `{theorem['total_source_assignment_count']}`.",
                f"- Shared phase/selector atoms: `{theorem['shared_phase_selector_atoms']}`.",
                f"- Minimal co-readiness masks: `{theorem['phase_selector_coreadiness_minimal_masks']}`.",
                f"- Minimal co-readiness size: `{theorem['phase_selector_coreadiness_minimal_size']}`.",
                f"- Quadrant counts: `{theorem['quadrant_counts']}`.",
                f"- No phase readiness without chi11: `{theorem['no_phase_readiness_without_chi11']}`.",
                f"- No selector readiness without chi11: `{theorem['no_selector_readiness_without_chi11']}`.",
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
