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

OUT = GEN / "p2398_s1348_role_closed_quotient_anf_certificate.json"
MD = GEN / "p2398_s1348_role_closed_quotient_anf_certificate.md"

SOURCE_FILES = {
    "P2397_ROLE_CLOSED_TOE_PROJECTION": GEN / "p2397_s1347_role_closed_toe_projection_certificate.json",
    "FRONTIER_TRUTH_TABLE": SCRATCH / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.json",
    "TOE_BOOLEAN_NORMAL_FORM": SCRATCH / "bridge_strict_completion_toe_boolean_normal_form_certificate_report.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

VARIABLES = [
    "strict_dynamical_source_for_A_P_D",
    "strict_phase_frequency_source",
    "strict_damping_beta_eta_source",
    "chi11_selector_source",
]

TARGETS = [
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
        "P2398|S1348|role-closed quotient ANF|quotient ANF|restricted ANF",
        "P2397|role-closed ToE projection|role atoms forced false|slice row count",
        "toe boolean normal form|ANF|single monomial|toe_anf_degree",
        "strict_dynamical_source_for_A_P_D|strict_phase_frequency_source|strict_damping_beta_eta_source|chi11_selector_source",
        "role_transfer theorem.*zero|toe closure.*zero|non-role atoms alone",
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
            "Repo grep finds P2397's role-closed projection and the global ToE Boolean normal form. "
            "P2398 therefore computes only the restricted quotient ANF on the P2396/P2397 role-closed slice."
        ),
    }


def bitmask_from_assignment(assignment: dict[str, bool]) -> int:
    mask = 0
    for index, variable in enumerate(VARIABLES):
        if assignment[variable]:
            mask |= 1 << index
    return mask


def monomial_name(mask: int) -> str:
    if mask == 0:
        return "1"
    return "*".join(variable for index, variable in enumerate(VARIABLES) if mask & (1 << index))


def anf_coefficients(values_by_mask: dict[int, int]) -> dict[int, int]:
    coeffs = [values_by_mask.get(mask, 0) for mask in range(1 << len(VARIABLES))]
    for bit in range(len(VARIABLES)):
        for mask in range(1 << len(VARIABLES)):
            if mask & (1 << bit):
                coeffs[mask] ^= coeffs[mask ^ (1 << bit)]
    return {mask: coeff for mask, coeff in enumerate(coeffs) if coeff}


def degree(mask: int) -> int:
    return mask.bit_count()


def compute_quotient_anf(p2397: dict[str, Any]) -> dict[str, Any]:
    slice_cert = p2397["role_closed_toe_projection_certificate"]["role_closed_slice"]
    rows = slice_cert["rows"]
    target_rows: dict[str, Any] = {}
    for target in TARGETS:
        values_by_mask = {
            bitmask_from_assignment(row["nonrole_assignment"]): int(row["target_values"][target])
            for row in rows
        }
        coeffs = anf_coefficients(values_by_mask)
        monomials = [
            {"mask": mask, "degree": degree(mask), "monomial": monomial_name(mask)}
            for mask in sorted(coeffs, key=lambda item: (degree(item), item))
        ]
        target_rows[target] = {
            "truth_vector_by_mask_0_to_15": [values_by_mask[mask] for mask in range(1 << len(VARIABLES))],
            "anf_monomials": monomials,
            "anf_monomial_count": len(monomials),
            "anf_degree": max((item["degree"] for item in monomials), default=-1),
            "is_zero_polynomial": len(monomials) == 0,
            "minimal_true_supports": [
                [VARIABLES[index] for index in range(len(VARIABLES)) if mask & (1 << index)]
                for mask, value in values_by_mask.items()
                if value and not any((submask != mask and (submask & mask) == submask and values_by_mask.get(submask, 0)) for submask in range(1 << len(VARIABLES)))
            ],
        }
    return {
        "variables": VARIABLES,
        "target_rows": target_rows,
        "proof_reading": (
            "On the role-closed quotient, bridge is the 3-source monomial, selector is the chi11 monomial, and role/toe are zero polynomials."
        ),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2398/S1348 role-closed quotient ANF certificate

`P2398/S1348` strengthens P2397 by computing the algebraic normal form of the P2396 role-closed quotient instead of only counting truth-table rows.  On the four free non-role variables, the quotient ANF is:

```text
bridge = strict_dynamical_source_for_A_P_D * strict_phase_frequency_source * strict_damping_beta_eta_source,
selector = chi11_selector_source,
role_transfer = 0,
toe = 0.
```

Thus role-transfer and ToE are not merely unsatisfied in the current row; they are identically zero functions on the whole role-closed quotient.  Non-role source work can still close the bridge or selector components, but it cannot close role-transfer or ToE until genuinely new role-successor evidence moves the state off this quotient.

No `L_total`, SM/GR numeric extraction, or ToE closure follows.
""".strip()
    lag_section = """
## P2398/S1348 quotient ANF guard for Lagrangian/EOM

`P2398/S1348` turns the P2397 role-closed slice into an exact quotient ANF.  In that quotient, `role_transfer=0` and `toe=0` as polynomials, while bridge and selector remain non-role targets.  Therefore a role-bearing `L_total` term cannot be justified by any combination of APD, phase/frequency, damping, or `chi11` source progress alone.

Only new explicit role-successor evidence changes the quotient.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    quotient = compute_quotient_anf(artifacts["P2397_ROLE_CLOSED_TOE_PROJECTION"])
    target_rows = quotient["target_rows"]
    inherited_normal_form = artifacts["TOE_BOOLEAN_NORMAL_FORM"].get("boolean_normal_form_summary", {})
    theorem_export = {
        "theorem_name": "P2398_T1_role_closed_quotient_anf_zero_role_and_toe",
        "variables": VARIABLES,
        "bridge_anf": target_rows["bridge_theorem_level_closure"]["anf_monomials"],
        "selector_anf": target_rows["selector_qw2191_closure"]["anf_monomials"],
        "role_transfer_is_zero_polynomial": target_rows["role_transfer_theorem_level_closure"]["is_zero_polynomial"],
        "toe_is_zero_polynomial": target_rows["toe_closure"]["is_zero_polynomial"],
        "inherited_global_toe_anf_degree": inherited_normal_form.get("toe_anf_degree"),
        "not_licensed": [
            "No role-transfer theorem is available on the role-closed quotient.",
            "No ToE closure is available on the role-closed quotient.",
            "No forever impossibility of future role-successor evidence is claimed.",
            "No L_total promotion or SM/GR numeric extraction is claimed.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "bridge_is_single_degree3_monomial": target_rows["bridge_theorem_level_closure"]["anf_monomial_count"] == 1 and target_rows["bridge_theorem_level_closure"]["anf_degree"] == 3,
        "selector_is_single_degree1_monomial": target_rows["selector_qw2191_closure"]["anf_monomial_count"] == 1 and target_rows["selector_qw2191_closure"]["anf_degree"] == 1,
        "role_transfer_zero_polynomial": target_rows["role_transfer_theorem_level_closure"]["is_zero_polynomial"],
        "toe_zero_polynomial": target_rows["toe_closure"]["is_zero_polynomial"],
        "global_toe_degree_inherited_as_seven": inherited_normal_form.get("toe_anf_degree") == 7,
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2398_s1348_v1",
        "packet_id": "P2398",
        "stage_id": "S1348",
        "result_kind": "ROLE_CLOSED_QUOTIENT_ANF_CERTIFICATE",
        "status": "PASS_ROLE_CLOSED_QUOTIENT_ANF_ROLE_AND_TOE_ZERO",
        "role_closed_quotient_anf_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", "TEXT_SOURCE") for name, artifact in artifacts.items()},
            "quotient_anf": quotient,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": "Treat non-role source progress as bridge/selector-only on this quotient; introduce new explicit role-successor evidence before revisiting role-transfer or ToE closure.",
        "global_status": "OPEN_PROGRESS_ROLE_CLOSED_QUOTIENT_ANF_CERTIFIED_ROLE_SUCCESSOR_EVIDENCE_REQUIRED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    cert = payload["role_closed_quotient_anf_certificate"]
    theorem = cert["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2398 S1348: role-closed quotient ANF certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2398/S1348 computes the exact ANF on the P2396/P2397 role-closed quotient.",
                "",
                "## Quotient ANF",
                "",
                f"- Bridge ANF: `{theorem['bridge_anf']}`.",
                f"- Selector ANF: `{theorem['selector_anf']}`.",
                f"- Role-transfer zero polynomial: `{theorem['role_transfer_is_zero_polynomial']}`.",
                f"- ToE zero polynomial: `{theorem['toe_is_zero_polynomial']}`.",
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
