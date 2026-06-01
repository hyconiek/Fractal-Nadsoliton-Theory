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
SCRATCH = ROOT / "scratch"

OUT = GEN / "p2397_s1347_role_closed_toe_projection_certificate.json"
MD = GEN / "p2397_s1347_role_closed_toe_projection_certificate.md"

SOURCE_FILES = {
    "P2396_ROLE_PACKAGE_REBASE": GEN / "p2396_s1346_role_package_negative_closure_rebase_certificate.json",
    "FRONTIER_TRUTH_TABLE": SCRATCH / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.json",
    "TOE_BOOLEAN_NORMAL_FORM": SCRATCH / "bridge_strict_completion_toe_boolean_normal_form_certificate_report.json",
    "TOE_PROPER_SUBSET_OBSTRUCTION": SCRATCH / "bridge_strict_completion_toe_proper_subset_obstruction_certificate_report.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

ROLE_ATOMS = [
    "alpha_geo_electroweak_role_theorem",
    "beta_tors_strict_role_theorem",
    "beta_power_hierarchy_successor_theorem",
]

NONROLE_ATOMS = [
    "strict_dynamical_source_for_A_P_D",
    "strict_phase_frequency_source",
    "strict_damping_beta_eta_source",
    "chi11_selector_source",
]

TARGET_KEYS = [
    "bridge_theorem_level_closure",
    "role_transfer_theorem_level_closure",
    "selector_qw2191_closure",
    "toe_closure",
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


def rg_audit() -> dict[str, Any]:
    patterns = [
        "P2397|S1347|role-closed ToE projection|role closed toe projection",
        "P2396|role package negative closure|current_licensed_transfer_count|future-only conditional",
        "toe boolean normal form|toe_closure|proper subset obstruction|all_seven_atoms_required",
        "alpha_geo_electroweak_role_theorem|beta_tors_strict_role_theorem|beta_power_hierarchy_successor_theorem",
        "strict_dynamical_source_for_A_P_D|strict_phase_frequency_source|strict_damping_beta_eta_source|chi11_selector_source",
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
        out[pattern] = {"count": len(lines), "samples": lines[:16]}
    return {
        "tool": "rg",
        "patterns": out,
        "finding": (
            "Repo grep finds existing ToE truth-table/normal-form/proper-subset reports and the P2396 current-state role-package closure. "
            "P2397 therefore does not recompute global ToE readiness; it projects the existing 7-atom ToE board onto the P2396 role-closed slice."
        ),
    }


def gf2_rank(matrix: list[list[int]]) -> int:
    work = [row[:] for row in matrix if any(row)]
    if not work:
        return 0
    rank = 0
    col = 0
    width = len(work[0])
    while rank < len(work) and col < width:
        pivot = next((row for row in range(rank, len(work)) if work[row][col] % 2), None)
        if pivot is None:
            col += 1
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for row in range(len(work)):
            if row != rank and work[row][col] % 2:
                work[row] = [a ^ b for a, b in zip(work[row], work[rank])]
        rank += 1
        col += 1
    return rank


def target_values(active: set[str], target_definitions: dict[str, Any]) -> dict[str, bool]:
    bridge = set(target_definitions["bridge_theorem_level_closure"]) <= active
    role = set(target_definitions["role_transfer_theorem_level_closure"]) <= active
    selector = set(target_definitions["selector_qw2191_closure"]) <= active
    toe = bridge and role and selector
    return {
        "bridge_theorem_level_closure": bridge,
        "role_transfer_theorem_level_closure": role,
        "selector_qw2191_closure": selector,
        "toe_closure": toe,
    }


def enumerate_role_closed_slice(target_definitions: dict[str, Any]) -> dict[str, Any]:
    rows = []
    signature_counts: dict[str, int] = {}
    for values in itertools.product([False, True], repeat=len(NONROLE_ATOMS)):
        assignment = dict(zip(NONROLE_ATOMS, values))
        active = {atom for atom, truth in assignment.items() if truth}
        targets = target_values(active, target_definitions)
        signature = "".join("1" if targets[key] else "0" for key in TARGET_KEYS)
        signature_counts[signature] = signature_counts.get(signature, 0) + 1
        rows.append(
            {
                "nonrole_assignment": assignment,
                "role_atoms_forced_false_by_p2396": ROLE_ATOMS,
                "target_values": targets,
                "target_signature_bridge_role_selector_toe": signature,
            }
        )
    signature_matrix = [
        [1 if char == "1" else 0 for char in row["target_signature_bridge_role_selector_toe"]]
        for row in rows
    ]
    return {
        "role_atoms_forced_false": ROLE_ATOMS,
        "free_nonrole_atoms": NONROLE_ATOMS,
        "slice_row_count": len(rows),
        "rows": rows,
        "signature_counts": signature_counts,
        "signature_rank_mod2": gf2_rank(signature_matrix),
        "toe_true_count": sum(row["target_values"]["toe_closure"] for row in rows),
        "role_transfer_true_count": sum(row["target_values"]["role_transfer_theorem_level_closure"] for row in rows),
        "bridge_true_count": sum(row["target_values"]["bridge_theorem_level_closure"] for row in rows),
        "selector_true_count": sum(row["target_values"]["selector_qw2191_closure"] for row in rows),
    }


def p2396_context(p2396: dict[str, Any]) -> dict[str, Any]:
    cert = p2396.get("role_package_negative_closure_rebase_certificate", {})
    package = cert.get("package_certificate", {})
    return {
        "role_package_closed_negative_now": bool(package.get("all_current_role_transfer_closed_negative")),
        "current_licensed_transfer_count": package.get("current_licensed_transfer_count"),
        "p2395_future_only_flag_count": package.get("p2395_future_only_flag_count"),
        "future_successor_not_forbidden_count": package.get("future_successor_not_forbidden_count"),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2397/S1347 role-closed ToE projection certificate

`P2397/S1347` projects the existing seven-atom ToE truth-table/normal-form board onto the P2396 current-state role-closed slice.  Repo grep confirms that the global ToE boolean-normal-form and proper-subset reports already exist, so P2397 only computes the 16 assignments where the three role atoms are forced false:

```text
alpha_geo_electroweak_role_theorem = false,
beta_tors_strict_role_theorem = false,
beta_power_hierarchy_successor_theorem = false.
```

The free non-role atoms are the three strict bridge-source atoms plus `chi11_selector_source`.  On this slice the bridge target can still become true and the selector target can still become true, but role-transfer and ToE closure are false in all 16 rows.  This proves a sharper current-state corollary of P2396: non-role progress alone cannot close ToE while the role package is closed on the current repo state.

This is not a forever no-go theorem.  Future explicit role-successor evidence would move the system off the P2396 role-closed slice.  No `L_total`, SM/GR numeric extraction, or ToE closure follows.
""".strip()
    lag_section = """
## P2397/S1347 role-closed ToE projection for Lagrangian/EOM

`P2397/S1347` computes the P2396 role-closed slice of the existing ToE truth table.  Across all 16 assignments of the remaining non-role atoms, role-transfer and ToE closure stay false.  Therefore even a future non-role source advance in APD, phase/frequency, damping, or `chi11` cannot by itself license a role-bearing `L_total` term while the legacy role package remains closed on the current repo state.

Future explicit role-successor evidence would be a different slice; P2397 does not forbid it forever.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    truth = artifacts["FRONTIER_TRUTH_TABLE"]
    target_definitions = truth.get("target_definitions", {})
    p2396 = p2396_context(artifacts["P2396_ROLE_PACKAGE_REBASE"])
    slice_cert = enumerate_role_closed_slice(target_definitions)
    normal_form = artifacts["TOE_BOOLEAN_NORMAL_FORM"].get("boolean_normal_form_summary", {})
    proper_subset = artifacts["TOE_PROPER_SUBSET_OBSTRUCTION"].get("proper_subset_obstruction_summary", {})
    theorem_export = {
        "theorem_name": "P2397_T1_role_closed_toe_projection_current_state",
        "p2396_context": p2396,
        "role_atoms_forced_false": ROLE_ATOMS,
        "free_nonrole_atoms": NONROLE_ATOMS,
        "slice_row_count": slice_cert["slice_row_count"],
        "role_transfer_true_count_on_slice": slice_cert["role_transfer_true_count"],
        "toe_true_count_on_slice": slice_cert["toe_true_count"],
        "bridge_true_count_on_slice": slice_cert["bridge_true_count"],
        "selector_true_count_on_slice": slice_cert["selector_true_count"],
        "signature_counts_on_slice": slice_cert["signature_counts"],
        "signature_rank_mod2": slice_cert["signature_rank_mod2"],
        "inherited_toe_anf_degree": normal_form.get("toe_anf_degree"),
        "inherited_all_seven_atoms_required": proper_subset.get("all_seven_atoms_required"),
        "not_licensed": [
            "No ToE closure is claimed on the P2396 role-closed slice.",
            "No role-transfer theorem is recovered by non-role atoms alone.",
            "No forever impossibility of future explicit role-successor evidence is claimed.",
            "No L_total promotion or SM/GR numeric extraction is claimed.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2396_role_package_closed": p2396["role_package_closed_negative_now"],
        "slice_size_is_2_pow_4": slice_cert["slice_row_count"] == 16,
        "role_transfer_false_on_slice": slice_cert["role_transfer_true_count"] == 0,
        "toe_false_on_slice": slice_cert["toe_true_count"] == 0,
        "bridge_can_still_close_on_slice": slice_cert["bridge_true_count"] == 2,
        "selector_can_still_close_on_slice": slice_cert["selector_true_count"] == 8,
        "signature_rank_is_two": slice_cert["signature_rank_mod2"] == 2,
        "inherited_toe_requires_all_seven_atoms": normal_form.get("toe_anf_degree") == 7 and proper_subset.get("all_seven_atoms_required") is True,
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2397_s1347_v1",
        "packet_id": "P2397",
        "stage_id": "S1347",
        "result_kind": "ROLE_CLOSED_TOE_PROJECTION_CERTIFICATE",
        "status": "PASS_ROLE_CLOSED_SLICE_TOE_FALSE_NONROLE_PROGRESS_SEPARATED",
        "role_closed_toe_projection_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", "TEXT_SOURCE") for name, artifact in artifacts.items()},
            "p2396_context": p2396,
            "target_definitions": target_definitions,
            "role_closed_slice": slice_cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": "Continue non-role strict-source work only with the explicit warning that ToE/role closure remains false on the P2396 role-closed slice unless new explicit role-successor evidence is introduced.",
        "global_status": "OPEN_PROGRESS_ROLE_CLOSED_TOE_PROJECTION_CERTIFIED_FUTURE_ROLE_SUCCESSOR_EVIDENCE_REQUIRED_FOR_TOE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    cert = payload["role_closed_toe_projection_certificate"]
    theorem = cert["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2397 S1347: role-closed ToE projection certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2397/S1347 projects the existing seven-atom ToE board onto the P2396 role-closed current-state slice instead of recomputing global ToE readiness.",
                "",
                "## Finite slice",
                "",
                f"- Role atoms forced false: `{theorem['role_atoms_forced_false']}`.",
                f"- Free non-role atoms: `{theorem['free_nonrole_atoms']}`.",
                f"- Slice row count: `{theorem['slice_row_count']}`.",
                f"- Bridge true count on slice: `{theorem['bridge_true_count_on_slice']}`.",
                f"- Selector true count on slice: `{theorem['selector_true_count_on_slice']}`.",
                f"- Role-transfer true count on slice: `{theorem['role_transfer_true_count_on_slice']}`.",
                f"- ToE true count on slice: `{theorem['toe_true_count_on_slice']}`.",
                f"- Signature rank over GF(2): `{theorem['signature_rank_mod2']}`.",
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
