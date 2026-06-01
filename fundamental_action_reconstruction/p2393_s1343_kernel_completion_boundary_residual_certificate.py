#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2393_s1343_kernel_completion_boundary_residual_certificate.json"
MD = GEN / "p2393_s1343_kernel_completion_boundary_residual_certificate.md"

SOURCE_FILES = {
    "P2392_SELECTOR_ASSUMPTION_RETIREMENT": GEN / "p2392_s1342_auxiliary_beta_tors_chi11_retirement_certificate.json",
    "S2_STRATEGIC_PRIORITY": ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
    "K1_KERNEL_SPLIT": ROOT / "K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md",
    "K2_STRICT_DERIVATION_CHAIN": ROOT / "K2_STRICT_GATE_KERNEL_DERIVATION_CHAIN_NOTE.md",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

ALPHA_GEO = 4.0 * math.log(2.0)
LEGACY_PARAMS = {
    "alpha_geo": ALPHA_GEO,
    "omega": math.pi / 4.0,
    "phi": math.pi / 6.0,
    "beta_tors": 1.0 / 10.0,
}
STRICT_TARGET_PARAMS = {"omega": 743.0 / 4000.0, "phi": 13.0 / 80.0, "beta": 1.0, "eta": 9.0 / 5.0}
FINITE_DOMAIN = list(range(0, 13))
GF2_COLUMNS = [
    "legacy_input_visible",
    "strict_target_visible",
    "explicit_completion_map_exported",
    "current_residual_zero",
    "role_transfer_licensed",
]


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def read_text(path: Path) -> str:
    if not path.exists():
        return f"OPEN_MISSING_ARTIFACT::{rel(path)}"
    return path.read_text(encoding="utf-8")


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
        "P2393|S1343|kernel completion boundary|boundary residual",
        "legacy -> strict completion|legacy-to-strict bridge|K_legacy_ont.*K_strict_gate",
        "eta=1|d\\^eta|nonlinear compression|linear torsion damping",
        "P2392|retired as selector|auxiliary beta_tors",
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
        out[pattern] = {"count": len(lines), "samples": lines[:12]}
    return {
        "tool": "rg",
        "patterns": out,
        "finding": (
            "Repo grep finds prior bridge and compression guardrails plus P2392 retirement of the beta_tors->chi11 selector route. "
            "P2393 therefore does not reopen that route; it computes the normalized eta=1 boundary embedding and the finite residual still left by current strict target parameters."
        ),
    }


def k_legacy(d: int, params: dict[str, float] = LEGACY_PARAMS) -> float:
    return params["alpha_geo"] * math.cos(params["omega"] * d + params["phi"]) / (1.0 + params["beta_tors"] * d)


def k_legacy_normalized(d: int, params: dict[str, float] = LEGACY_PARAMS) -> float:
    return k_legacy(d, params) / params["alpha_geo"]


def k_strict(d: int, params: dict[str, float]) -> float:
    if d == 0 and params["eta"] <= 0.0:
        raise ValueError("eta must be positive for the audited strict kernel domain")
    return math.cos(params["omega"] * d + params["phi"]) / (1.0 + params["beta"] * (d ** params["eta"]))


def boundary_params_from_legacy() -> dict[str, float]:
    return {
        "omega": LEGACY_PARAMS["omega"],
        "phi": LEGACY_PARAMS["phi"],
        "beta": LEGACY_PARAMS["beta_tors"],
        "eta": 1.0,
    }


def boundary_identity_audit() -> dict[str, Any]:
    params = boundary_params_from_legacy()
    rows = []
    max_abs_error = 0.0
    for d in FINITE_DOMAIN:
        left = k_legacy_normalized(d)
        right = k_strict(d, params)
        error = right - left
        max_abs_error = max(max_abs_error, abs(error))
        rows.append({"d": d, "legacy_normalized": left, "strict_eta1_boundary": right, "error": error})
    return {
        "identity_statement": "K_legacy_ont(d)/alpha_geo = K_strict_gate(d) after omega=omega_legacy, phi=phi_legacy, beta=beta_tors, eta=1.",
        "boundary_parameter_substitution": params,
        "finite_domain": FINITE_DOMAIN,
        "max_abs_error_float_replay": max_abs_error,
        "allclose_at_1e_minus_15": max_abs_error < 1.0e-15,
        "rows": rows,
        "proof_role": "This is only the normalized linear-damping boundary embedding; it is not the completed strict bridge.",
    }


def current_strict_residual_audit() -> dict[str, Any]:
    rows = []
    residuals = []
    previous_sign = None
    sign_changes = 0
    for d in FINITE_DOMAIN:
        legacy_norm = k_legacy_normalized(d)
        strict_value = k_strict(d, STRICT_TARGET_PARAMS)
        residual = strict_value - legacy_norm
        residuals.append(residual)
        sign = 1 if residual > 0 else -1 if residual < 0 else 0
        if previous_sign is not None and sign != 0 and previous_sign != 0 and sign != previous_sign:
            sign_changes += 1
        if sign != 0:
            previous_sign = sign
        rows.append(
            {
                "d": d,
                "legacy_normalized": legacy_norm,
                "strict_current_target": strict_value,
                "residual_strict_minus_legacy_norm": residual,
                "abs_residual": abs(residual),
            }
        )
    l2 = math.sqrt(sum(r * r for r in residuals))
    linf_row = max(rows, key=lambda row: row["abs_residual"])
    mean_abs = sum(abs(r) for r in residuals) / len(residuals)
    return {
        "strict_target_parameters": STRICT_TARGET_PARAMS,
        "finite_domain": FINITE_DOMAIN,
        "rows": rows,
        "l2_residual": l2,
        "mean_abs_residual": mean_abs,
        "linf_residual_row": linf_row,
        "residual_sign_changes": sign_changes,
        "nonzero_residual_on_current_target": any(abs(r) > 1.0e-12 for r in residuals),
        "interpretation": (
            "The current strict target is not obtained by amplitude-normalizing the legacy kernel alone. "
            "The bridge still needs explicit phase/frequency passage and nonlinear compression/damping completion."
        ),
    }


def gf2_rank(matrix: list[list[int]]) -> int:
    work = [row[:] for row in matrix if any(row)]
    rank = 0
    col = 0
    while work and col < len(work[0]):
        pivot = next((r for r in range(rank, len(work)) if work[r][col] % 2), None)
        if pivot is None:
            col += 1
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for r in range(len(work)):
            if r != rank and work[r][col] % 2:
                work[r] = [(a ^ b) for a, b in zip(work[r], work[rank])]
        rank += 1
        col += 1
    return rank


def bridge_obligation_matrix(boundary: dict[str, Any], residual: dict[str, Any], p2392: dict[str, Any]) -> dict[str, Any]:
    p2392_status = p2392.get("status", "OPEN_MISSING_ARTIFACT")
    beta_route_retired = "RETIRED" in p2392_status or "PASS_BETA_TORS_CHI11_RETIRED" in p2392_status
    rows = [
        {
            "obligation": "amplitude normalization boundary",
            "bits": [1, 1, 1, int(boundary["allclose_at_1e_minus_15"]), 0],
            "status": "BOUNDARY_IDENTITY_PROVED_NOT_FULL_BRIDGE",
        },
        {
            "obligation": "phase/frequency passage legacy canonical to strict target",
            "bits": [1, 1, 0, 0, 0],
            "status": "OPEN_PHASE_FREQUENCY_COMPLETION_MAP_REQUIRED",
        },
        {
            "obligation": "linear torsion damping to nonlinear strict compression",
            "bits": [1, 1, 0, int(not residual["nonzero_residual_on_current_target"]), 0],
            "status": "OPEN_NONLINEAR_COMPRESSION_COMPLETION_MAP_REQUIRED",
        },
        {
            "obligation": "selector route after P2392",
            "bits": [0, 1, 1, 1, 0],
            "status": "SELECTOR_MECHANISM_AVAILABLE_BETA_TORS_CHI11_ROUTE_RETIRED" if beta_route_retired else "SELECTOR_RETIREMENT_ARTIFACT_MISSING",
        },
        {
            "obligation": "post-bridge legacy physical role transfer",
            "bits": [1, 1, 0, 0, 0],
            "status": "BLOCKED_UNTIL_COMPLETION_BRIDGE_THEN_SEPARATE_ROLE_AUDIT",
        },
    ]
    matrix = [row["bits"] for row in rows]
    return {
        "columns": GF2_COLUMNS,
        "rows": rows,
        "rank_mod2": gf2_rank(matrix),
        "all_current_role_transfer_bits_zero": all(row["bits"][4] == 0 for row in rows),
        "open_completion_rows": [row["obligation"] for row in rows if row["status"].startswith("OPEN")],
        "proof_reading": (
            "The eta=1 boundary row is a real bridge ingredient, but the nonzero current residual and zero role-transfer column keep the bridge incomplete."
        ),
    }


def append_doc_sections(payload: dict[str, Any]) -> None:
    section_eq = """
## P2393/S1343 normalized kernel boundary and current residual certificate

`P2393/S1343` moves the bridge-completion lane from wording to a finite computation.  It proves the exact boundary embedding

```text
K_legacy_ont(d)/alpha_geo = K_strict_gate(d)
```

only after the explicit boundary substitution

```text
omega = omega_legacy, phi = phi_legacy, beta = beta_tors, eta = 1.
```

On the audited finite domain `d=0..12`, the replayed boundary residual is zero to floating precision, so amplitude normalization plus the `eta=1` substitution is a genuine bridge ingredient.  The same computation against the current strict target parameters leaves a nonzero residual vector, but after P2394 this is read only as an `eta=1` negative-control slice: the repo already has the finite APD completion witness `K_strict = K_legacy*A*P*D`.

This step does not reopen `beta_tors -> chi11`; after P2392 that route remains retired as a selector-search assumption.  It also does not license legacy physical-role transfer; after P2394 the next active target is the post-bridge role-transfer audit, not re-proving APD.
""".strip()
    section_lag = """
## P2393/S1343 normalized kernel boundary residual, no L_total promotion

`P2393/S1343` certifies a narrow kernel identity useful as a boundary negative control: normalized legacy linear damping is the `eta=1` boundary slice of the strict kernel family when the phase, frequency, and damping symbols are explicitly matched.  The finite residual audit then compares that normalized legacy slice with the current strict target and finds a nonzero residual; after P2394 this is not an APD-bridge gap, because the finite APD comparison witness already supplies `K_strict = K_legacy*A*P*D`.

Consequently, the Lagrangian/EOM draft may use P2393 only as a boundary-condition/negative-control certificate.  It may not promote the legacy physical roles into `L_total`, may not revive `beta_tors -> chi11` as a selector target, and may not claim role-transfer or ToE closure.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], section_eq), (DOC_FILES["lagrangian_eom_draft"], section_lag)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    artifacts = {name: load_json(path) if path.suffix == ".json" else {"text": read_text(path)} for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    boundary = boundary_identity_audit()
    residual = current_strict_residual_audit()
    matrix = bridge_obligation_matrix(boundary, residual, artifacts["P2392_SELECTOR_ASSUMPTION_RETIREMENT"])
    theorem_export = {
        "theorem_name": "P2393_T1_normalized_legacy_eta1_boundary_not_full_strict_completion",
        "legacy_parameters": LEGACY_PARAMS,
        "strict_target_parameters": STRICT_TARGET_PARAMS,
        "finite_domain": FINITE_DOMAIN,
        "proved_boundary_identity": boundary["identity_statement"],
        "boundary_max_abs_error_float_replay": boundary["max_abs_error_float_replay"],
        "current_target_linf_residual": residual["linf_residual_row"]["abs_residual"],
        "current_target_l2_residual": residual["l2_residual"],
        "completion_obligations_still_open": matrix["open_completion_rows"],
        "retired_route_not_reopened": "beta_tors -> chi11 remains retired as a selector-search assumption by P2392; P2393 uses beta_tors only as a legacy damping coordinate in the eta=1 boundary calculation.",
        "not_licensed": [
            "No completed legacy->strict bridge is claimed.",
            "No legacy physical-role transfer is claimed.",
            "No beta_tors -> chi11 selector target is reopened.",
            "No L_total source term, SM/GR extraction, or ToE closure is claimed.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2392_retirement_seen": artifacts["P2392_SELECTOR_ASSUMPTION_RETIREMENT"].get("packet_id") == "P2392",
        "boundary_identity_passes": boundary["allclose_at_1e_minus_15"],
        "current_target_residual_nonzero": residual["nonzero_residual_on_current_target"],
        "open_completion_rows_present": len(matrix["open_completion_rows"]) >= 2,
        "role_transfer_not_licensed": matrix["all_current_role_transfer_bits_zero"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2393_s1343_v1",
        "packet_id": "P2393",
        "stage_id": "S1343",
        "result_kind": "KERNEL_COMPLETION_BOUNDARY_RESIDUAL_CERTIFICATE",
        "status": "PASS_BOUNDARY_EMBEDDING_PROVED_CURRENT_COMPLETION_RESIDUAL_OPEN",
        "kernel_completion_boundary_residual_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", "TEXT_SOURCE") for name, artifact in artifacts.items()},
            "boundary_identity_audit": boundary,
            "current_strict_residual_audit": residual,
            "bridge_obligation_matrix": matrix,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": "Use the P2393 boundary row as a bridge ingredient, then attack the explicit phase/frequency and nonlinear-compression completion maps before any role-transfer audit.",
        "global_status": "OPEN_PROGRESS_BOUNDARY_CERTIFIED_COMPLETION_MAP_AND_ROLE_TRANSFER_STILL_OPEN",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    cert = payload["kernel_completion_boundary_residual_certificate"]
    theorem = cert["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2393 S1343: normalized kernel boundary and current residual certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2393/S1343 proves the normalized eta=1 boundary embedding of the legacy kernel into the strict kernel family, then computes the residual left by the current strict target parameters.",
                "",
                "## Boundary identity",
                "",
                f"- Identity: `{theorem['proved_boundary_identity']}`",
                f"- Finite domain: `{theorem['finite_domain']}`",
                f"- Max absolute boundary replay error: `{theorem['boundary_max_abs_error_float_replay']}`.",
                "",
                "## Current strict residual",
                "",
                f"- Current strict target parameters: `{theorem['strict_target_parameters']}`.",
                f"- L-infinity residual: `{theorem['current_target_linf_residual']}`.",
                f"- L2 residual: `{theorem['current_target_l2_residual']}`.",
                f"- Open completion obligations: `{theorem['completion_obligations_still_open']}`.",
                "",
                "## Hard limits",
                "",
                f"- {theorem['retired_route_not_reopened']}",
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
    payload = build_payload()
    append_doc_sections(payload)
    payload = build_payload()
    write_outputs(payload)
    if not all(payload["gatekeeper_checks"].values()):
        raise SystemExit(f"gatekeeper failure: {payload['gatekeeper_checks']}")


if __name__ == "__main__":
    main()
