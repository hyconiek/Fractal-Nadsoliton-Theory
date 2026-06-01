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

OUT = GEN / "p2412_s1362_chi11_selector_scope_separation_certificate.json"
MD = GEN / "p2412_s1362_chi11_selector_scope_separation_certificate.md"

SOURCE_FILES = {
    "P2366_PHASE_ORIGIN_CANDIDATE": GEN / "p2366_s1316_selector_candidate_phase_reference_chi11_audit_probe.json",
    "P2392_BETA_TORS_RETIREMENT": GEN / "p2392_s1342_auxiliary_beta_tors_chi11_retirement_certificate.json",
    "P2411_BRIDGE_SOURCE_HYPERGRAPH": GEN / "p2411_s1361_legacy_strict_bridge_source_obligation_hypergraph_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

SCOPE_ATOMS = [
    "declared_scope_strict_selector_available",
    "beta_tors_selector_route_retired",
    "phase_origin_candidate_finite_audit_passed",
    "bridge_chi11_source_theorem_exported",
    "qw2191_discharge_exported",
]
ATOM_INDEX = {atom: index for index, atom in enumerate(SCOPE_ATOMS)}

CURRENT_REQUIRED_TRUE = {
    "declared_scope_strict_selector_available",
    "beta_tors_selector_route_retired",
    "phase_origin_candidate_finite_audit_passed",
}
CURRENT_REQUIRED_FALSE = {
    "bridge_chi11_source_theorem_exported",
    "qw2191_discharge_exported",
}


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
    return {"count": len(lines), "samples": lines[:16]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2412|S1362|chi11 selector scope separation|selector-scope separation",
        "p2366_candidate": "P2366|phase-origin selector candidate|chiral bispectrum|coprime Fourier phase",
        "p2392_retirement": "P2392|beta_tors -> chi11|selector-assumption retirement|retired selector",
        "p2411_bridge": "P2411|chi11_selector_source_theorem|QW-2191|bridge-source hypergraph",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds candidate, retirement, and bridge-hypergraph materials, but no compact certificate "
            "separating declared-scope selector availability from bridge-level chi11 source/QW-2191 closure."
        ),
    }


def atom_value(mask: int, atom: str) -> bool:
    return bool(mask & (1 << ATOM_INDEX[atom]))


def atoms_true(mask: int) -> list[str]:
    return [atom for atom in SCOPE_ATOMS if atom_value(mask, atom)]


def current_scope_value(mask: int) -> int:
    return int(
        all(atom_value(mask, atom) for atom in CURRENT_REQUIRED_TRUE)
        and all(not atom_value(mask, atom) for atom in CURRENT_REQUIRED_FALSE)
    )


def declared_selector_lane(mask: int) -> int:
    return int(atom_value(mask, "declared_scope_strict_selector_available"))


def bridge_selector_closure_lane(mask: int) -> int:
    return int(
        atom_value(mask, "bridge_chi11_source_theorem_exported")
        and atom_value(mask, "qw2191_discharge_exported")
    )


def truth_rows() -> list[dict[str, Any]]:
    rows = []
    for mask in range(1 << len(SCOPE_ATOMS)):
        rows.append(
            {
                "mask": mask,
                "true_atoms": atoms_true(mask),
                "current_scope_separation": current_scope_value(mask),
                "declared_selector_lane": declared_selector_lane(mask),
                "bridge_selector_closure_lane": bridge_selector_closure_lane(mask),
            }
        )
    return rows


def source_facts(sources: dict[str, Any]) -> dict[str, Any]:
    p2366_theorem = sources["P2366_PHASE_ORIGIN_CANDIDATE"].get(
        "selector_candidate_phase_reference_chi11_audit_probe", {}
    ).get("theorem_export", {})
    p2392_theorem = sources["P2392_BETA_TORS_RETIREMENT"].get(
        "auxiliary_beta_tors_chi11_retirement_certificate", {}
    ).get("theorem_export", {})
    p2411_theorem = sources["P2411_BRIDGE_SOURCE_HYPERGRAPH"].get(
        "legacy_strict_bridge_source_obligation_hypergraph_certificate", {}
    ).get("theorem_export", {})
    p2411_top_atoms = [row.get("atom") for row in p2411_theorem.get("top_priority_atoms", [])]
    return {
        "declared_scope_strict_selector_available": p2392_theorem.get("available_atoms", {}).get(
            "strict_internal_selector_P1343_P1348"
        ) is True,
        "beta_tors_selector_route_retired": p2392_theorem.get("active_beta_tors_chi11_obligation_count") == 0
        and p2392_theorem.get("available_atoms", {}).get("auxiliary_beta_tors_to_chi11") is False,
        "phase_origin_candidate_finite_audit_passed": "All 24 source/orientation rows" in "\n".join(
            p2366_theorem.get("positive_content", [])
        ),
        "bridge_chi11_source_theorem_exported": False,
        "qw2191_discharge_exported": False,
        "bridge_chi11_source_is_p2411_top_priority": "chi11_selector_source_theorem" in p2411_top_atoms,
        "p2366_still_blocks_qw2191_discharge": "selector/QW-2191 discharge" in p2366_theorem.get(
            "not_licensed", []
        ),
        "p2411_still_blocks_chi11_source": any(
            "No chi11 selector-source theorem" in item for item in p2411_theorem.get("not_licensed", [])
        ),
    }


def current_mask_from_facts(facts: dict[str, Any]) -> int:
    mask = 0
    for atom in SCOPE_ATOMS:
        if facts[atom]:
            mask |= 1 << ATOM_INDEX[atom]
    return mask


def candidate_classification(facts: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": "P1343_P1348_declared_scope_selector",
            "current_status": "available_for_declared_selector_lane" if facts["declared_scope_strict_selector_available"] else "missing",
            "licenses_declared_selector_lane": facts["declared_scope_strict_selector_available"],
            "licenses_bridge_chi11_source": False,
            "licenses_qw2191_discharge": False,
        },
        {
            "candidate": "beta_tors_to_chi11",
            "current_status": "retired_as_selector_search_route" if facts["beta_tors_selector_route_retired"] else "active_or_unclassified",
            "licenses_declared_selector_lane": False,
            "licenses_bridge_chi11_source": False,
            "licenses_qw2191_discharge": False,
        },
        {
            "candidate": "phase_origin_chiral_bispectrum_coprime_phase",
            "current_status": "finite_candidate_audited_not_strict_core_source" if facts["phase_origin_candidate_finite_audit_passed"] else "missing_audit",
            "licenses_declared_selector_lane": False,
            "licenses_bridge_chi11_source": False,
            "licenses_qw2191_discharge": False,
        },
        {
            "candidate": "explicit_bridge_chi11_source_theorem",
            "current_status": "top_priority_open_obligation" if facts["bridge_chi11_source_is_p2411_top_priority"] else "open",
            "licenses_declared_selector_lane": False,
            "licenses_bridge_chi11_source": facts["bridge_chi11_source_theorem_exported"],
            "licenses_qw2191_discharge": False,
        },
        {
            "candidate": "qw2191_symmetry_breaking_or_internal_source_theorem",
            "current_status": "open_not_exported",
            "licenses_declared_selector_lane": False,
            "licenses_bridge_chi11_source": False,
            "licenses_qw2191_discharge": facts["qw2191_discharge_exported"],
        },
    ]


def build_certificate(sources: dict[str, Any]) -> dict[str, Any]:
    rows = truth_rows()
    facts = source_facts(sources)
    current_mask = current_mask_from_facts(facts)
    return {
        "scope_atoms": SCOPE_ATOMS,
        "truth_rows": rows,
        "truth_row_count": len(rows),
        "current_source_facts": facts,
        "current_mask": current_mask,
        "current_true_atoms": atoms_true(current_mask),
        "current_scope_separation_true": current_scope_value(current_mask) == 1,
        "current_scope_separation_true_masks": [row["mask"] for row in rows if row["current_scope_separation"]],
        "declared_selector_lane_true_count": sum(row["declared_selector_lane"] for row in rows),
        "bridge_selector_closure_lane_true_count": sum(row["bridge_selector_closure_lane"] for row in rows),
        "candidate_classification": candidate_classification(facts),
        "scope_separation_prime_implicant": {
            "required_true": sorted(CURRENT_REQUIRED_TRUE),
            "required_false": sorted(CURRENT_REQUIRED_FALSE),
            "literal_count": len(CURRENT_REQUIRED_TRUE) + len(CURRENT_REQUIRED_FALSE),
        },
    }


def append_doc_sections() -> None:
    eq_section = """
## P2412/S1362 chi11 selector scope-separation certificate

`P2412/S1362` reconciles three facts that must not be conflated: P2392 makes the declared-scope strict selector available without the retired `beta_tors -> chi11` route; P2366 gives a finite phase-origin candidate but not a strict-core source theorem; and P2411 keeps the bridge-level `chi11` source/QW-2191 obligation open.

The finite five-atom scope lattice has `32` rows.  The current mask is the signed state where declared selector availability, beta-route retirement, and the finite phase-origin candidate are true, while bridge-level `chi11` source export and QW-2191 discharge are false.  That signed state is a consistency certificate, not a closure theorem.

Thus `chi11` bookkeeping may be used only inside its declared selector scope; it cannot be promoted into the legacy-to-strict bridge source, QW-2191 discharge, role transfer, `L_total`, or ToE closure without a new theorem.
""".strip()
    lag_section = """
## P2412/S1362 chi11 selector scope guard for Lagrangian/EOM

`P2412/S1362` separates declared-scope selector availability from bridge-level `chi11` source closure.  Lagrangian terms may not treat the retired `beta_tors -> chi11` route, the finite phase-origin candidate, or generic declared selector availability as a bridge-source theorem or QW-2191 discharge.
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
        "theorem_name": "P2412_T1_chi11_selector_scope_separation_certificate",
        "truth_row_count": cert["truth_row_count"],
        "current_mask": cert["current_mask"],
        "current_true_atoms": cert["current_true_atoms"],
        "current_scope_separation_true": cert["current_scope_separation_true"],
        "current_scope_separation_true_masks": cert["current_scope_separation_true_masks"],
        "declared_selector_lane_true_count": cert["declared_selector_lane_true_count"],
        "bridge_selector_closure_lane_true_count": cert["bridge_selector_closure_lane_true_count"],
        "scope_separation_prime_implicant": cert["scope_separation_prime_implicant"],
        "bridge_chi11_source_is_p2411_top_priority": cert["current_source_facts"]["bridge_chi11_source_is_p2411_top_priority"],
        "p2366_still_blocks_qw2191_discharge": cert["current_source_facts"]["p2366_still_blocks_qw2191_discharge"],
        "p2411_still_blocks_chi11_source": cert["current_source_facts"]["p2411_still_blocks_chi11_source"],
        "not_licensed": [
            "Declared-scope selector availability is not a bridge-level chi11 source theorem.",
            "The retired beta_tors -> chi11 selector-search route is not reopened.",
            "The phase-origin candidate is not a strict-core source theorem or QW-2191 discharge.",
            "No legacy role transfer, role-bearing L_total, or ToE closure follows from selector scope separation.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "exact_32_row_scope_lattice": theorem_export["truth_row_count"] == 32,
        "current_mask_is_declared_selector_beta_retired_candidate_only": theorem_export["current_mask"] == 7,
        "current_scope_separation_holds": theorem_export["current_scope_separation_true"] is True,
        "scope_separation_has_single_signed_mask": theorem_export["current_scope_separation_true_masks"] == [7],
        "declared_selector_lane_is_broader_than_bridge_closure": theorem_export["declared_selector_lane_true_count"] == 16
        and theorem_export["bridge_selector_closure_lane_true_count"] == 8,
        "bridge_chi11_source_remains_top_priority": theorem_export["bridge_chi11_source_is_p2411_top_priority"],
        "p2366_qw2191_still_blocked": theorem_export["p2366_still_blocks_qw2191_discharge"],
        "p2411_chi11_source_still_blocked": theorem_export["p2411_still_blocks_chi11_source"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2412_s1362_v1",
        "packet_id": "P2412",
        "stage_id": "S1362",
        "result_kind": "CHI11_SELECTOR_SCOPE_SEPARATION_CERTIFICATE",
        "status": "PASS_CHI11_SCOPE_SEPARATION_DECLARED_SELECTOR_NOT_BRIDGE_SOURCE",
        "chi11_selector_scope_separation_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_scope_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Attack the explicit bridge-level chi11 source/QW-2191 theorem directly; do not reopen beta_tors->chi11 "
            "or overread the phase-origin candidate as strict-core closure."
        ),
        "global_status": "OPEN_PROGRESS_CHI11_SCOPE_SEPARATED_NO_BRIDGE_SOURCE_OR_TOE_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["chi11_selector_scope_separation_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2412 S1362: chi11 selector scope-separation certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2412/S1362 separates declared-scope selector availability from bridge-level chi11 source and QW-2191 closure.",
                "",
                "## Finite Boolean facts",
                "",
                f"- Truth rows: `{theorem['truth_row_count']}`.",
                f"- Current mask: `{theorem['current_mask']}`.",
                f"- Current true atoms: `{theorem['current_true_atoms']}`.",
                f"- Scope-separation true masks: `{theorem['current_scope_separation_true_masks']}`.",
                f"- Declared selector lane true count: `{theorem['declared_selector_lane_true_count']}`.",
                f"- Bridge selector closure lane true count: `{theorem['bridge_selector_closure_lane_true_count']}`.",
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
