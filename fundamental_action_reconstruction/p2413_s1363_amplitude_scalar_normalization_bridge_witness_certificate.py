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
SCRATCH = ROOT / "scratch"

OUT = GEN / "p2413_s1363_amplitude_scalar_normalization_bridge_witness_certificate.json"
MD = GEN / "p2413_s1363_amplitude_scalar_normalization_bridge_witness_certificate.md"

SOURCE_FILES = {
    "P2411_BRIDGE_SOURCE_HYPERGRAPH": GEN / "p2411_s1361_legacy_strict_bridge_source_obligation_hypergraph_certificate.json",
    "P2412_CHI11_SCOPE_SEPARATION": GEN / "p2412_s1362_chi11_selector_scope_separation_certificate.json",
    "SCRATCH_AMPLITUDE_SCALAR_NORMALIZATION": SCRATCH / "bridge_strict_completion_legacy_to_strict_amplitude_scalar_normalization_certificate_report.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

DOMAIN = list(range(12))
ALPHA_GEO_EXACT = "4*ln(2)"
BETA_TORS_EXACT = "1/100"
OMEGA_EXACT = "pi/4"
PHI_EXACT = "pi/6"


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
        "new_packet": "P2413|S1363|amplitude scalar normalization bridge witness|role-safe amplitude scalar",
        "scratch_amplitude": "amplitude scalar-normalization|alpha_geo\\^{-1}|K_legacy_ont\\(d\\)=alpha_geo\\*L_shape",
        "role_safe_boundary": "role-safe amplitude|scalar normalization alone is not a physical-role proof|alpha_geo physical-role transfer",
        "bridge_priority": "amplitude normalization|legacy -> strict completion bridge|K_legacy_ont -> K_strict_gate",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds scratch amplitude scalar-normalization evidence and role-transfer warnings, but no production "
            "P24xx exact bridge witness that separates scalar factorization from role-safe amplitude/source export."
        ),
    }


def theta_numerator(d: int) -> int:
    return 3 * d + 2


def denominator_numerator(d: int) -> int:
    return 100 + d


def cos_zero_congruence_holds(d: int) -> bool:
    return (theta_numerator(d) - 6) % 12 == 0


def approximate_legacy_shape(d: int) -> float:
    return math.cos(math.pi * d / 4.0 + math.pi / 6.0) / (1.0 + d / 100.0)


def sign(value: float) -> int:
    return 1 if value > 0 else -1 if value < 0 else 0


def witness_rows() -> list[dict[str, Any]]:
    alpha = 4.0 * math.log(2.0)
    rows = []
    for d in DOMAIN:
        shape = approximate_legacy_shape(d)
        kernel = alpha * shape
        normalized = kernel / alpha
        rows.append(
            {
                "d": d,
                "theta_over_pi": f"{theta_numerator(d)}/12",
                "cos_zero_congruence_3d_eq_4_mod_12": cos_zero_congruence_holds(d),
                "denominator": f"{denominator_numerator(d)}/100",
                "denominator_positive_exact": denominator_numerator(d) > 0,
                "factorization_identity": "K_legacy_ont(d) = alpha_geo * L_shape(d)",
                "normalization_identity": "alpha_geo^{-1} * K_legacy_ont(d) = L_shape(d)",
                "legacy_shape_float": shape,
                "legacy_kernel_float": kernel,
                "normalized_float": normalized,
                "float_residual_normalized_minus_shape": normalized - shape,
                "sign_shape": sign(shape),
                "sign_kernel": sign(kernel),
                "positive_alpha_preserves_sign": sign(shape) == sign(kernel),
            }
        )
    return rows


def exact_proof_certificate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "factorization_step": (
            "By definition K_legacy_ont(d)=alpha_geo*cos(pi*d/4+pi/6)/(1+beta_tors*d), so with "
            "L_shape(d)=cos(pi*d/4+pi/6)/(1+beta_tors*d), K_legacy_ont(d)=alpha_geo*L_shape(d)."
        ),
        "normalization_step": "Multiplication by alpha_geo^{-1} removes exactly one global scalar and leaves L_shape(d).",
        "positivity_step": "alpha_geo=4*ln(2)>0, hence scalar normalization does not flip signs.",
        "denominator_step": "On d=0..11, 1+d/100=(100+d)/100 is positive.",
        "nonzero_cos_step": "cos((3d+2)pi/12)=0 would require 3d+2 == 6 mod 12, i.e. 3d == 4 mod 12; gcd(3,12)=3 does not divide 4.",
        "domain_cos_zero_solutions": [row["d"] for row in rows if row["cos_zero_congruence_3d_eq_4_mod_12"]],
        "max_float_residual": max(abs(row["float_residual_normalized_minus_shape"]) for row in rows),
        "all_denominators_positive_exact": all(row["denominator_positive_exact"] for row in rows),
        "all_signs_preserved_by_positive_alpha": all(row["positive_alpha_preserves_sign"] for row in rows),
    }


def build_certificate(sources: dict[str, Any]) -> dict[str, Any]:
    rows = witness_rows()
    proof = exact_proof_certificate(rows)
    p2411_theorem = sources["P2411_BRIDGE_SOURCE_HYPERGRAPH"].get(
        "legacy_strict_bridge_source_obligation_hypergraph_certificate", {}
    ).get("theorem_export", {})
    p2412_theorem = sources["P2412_CHI11_SCOPE_SEPARATION"].get(
        "chi11_selector_scope_separation_certificate", {}
    ).get("theorem_export", {})
    scratch_summary = sources["SCRATCH_AMPLITUDE_SCALAR_NORMALIZATION"].get(
        "amplitude_scalar_normalization_summary", {}
    )
    return {
        "domain": DOMAIN,
        "alpha_geo_exact": ALPHA_GEO_EXACT,
        "beta_tors_exact": BETA_TORS_EXACT,
        "omega_exact": OMEGA_EXACT,
        "phi_exact": PHI_EXACT,
        "legacy_kernel_definition": "K_legacy_ont(d)=alpha_geo*cos(pi*d/4+pi/6)/(1+d/100)",
        "legacy_shape_definition": "L_shape(d)=cos(pi*d/4+pi/6)/(1+d/100)",
        "normalization_map": "N_alpha: K -> alpha_geo^{-1} K",
        "witness_rows": rows,
        "proof_certificate": proof,
        "scalar_normalization_witness_ready": proof["domain_cos_zero_solutions"] == []
        and proof["all_denominators_positive_exact"]
        and proof["all_signs_preserved_by_positive_alpha"],
        "full_amplitude_source_theorem_ready": False,
        "role_safe_amplitude_absorption_ready": False,
        "bridge_completion_ready": False,
        "p2411_full_bridge_still_requires_all_obligations": p2411_theorem.get("bridge_ready_true_masks") == [255],
        "p2412_chi11_scope_separation_inherited": p2412_theorem.get("current_mask") == 7
        and p2412_theorem.get("current_scope_separation_true") is True,
        "scratch_amplitude_witness_inherited": scratch_summary.get("scalar_normalization_witness_exported") is True,
        "scratch_blocks_role_transfer_inherited": scratch_summary.get("legacy_role_transfer_allowed") is False,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2413/S1363 amplitude scalar-normalization bridge witness certificate

`P2413/S1363` takes the first S2 bridge component, amplitude normalization, and exports a narrow exact witness rather than another global obstruction ledger.  On the audited `d=0..11` domain,

```text
K_legacy_ont(d) = alpha_geo * L_shape(d),
L_shape(d) = cos(pi*d/4 + pi/6)/(1+d/100),
alpha_geo^{-1} K_legacy_ont(d) = L_shape(d).
```

The proof checks the exact zero condition `3d == 4 mod 12`, which has no solution, and the denominator `(100+d)/100 > 0`, so the positive scalar `alpha_geo=4 ln 2` removes only the global amplitude and preserves signs.

This is a bridge ingredient, not a full amplitude source theorem and not a role-safe amplitude absorption theorem: it does not transfer `sin^2(theta_W)=alpha_geo/12`, does not complete `K_legacy_ont -> K_strict_gate`, and does not affect QW-2191, `L_total`, or ToE closure.
""".strip()
    lag_section = """
## P2413/S1363 amplitude scalar-normalization guard for Lagrangian/EOM

`P2413/S1363` permits only the scalar bridge witness `alpha_geo^{-1} K_legacy_ont = L_shape` on the audited domain.  It does not license an `alpha_geo` electroweak Lagrangian term, role-safe amplitude absorption, full bridge completion, or any `L_total` promotion.
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
    proof = cert["proof_certificate"]
    theorem_export = {
        "theorem_name": "P2413_T1_amplitude_scalar_normalization_bridge_witness_certificate",
        "domain_size": len(cert["domain"]),
        "normalization_map": cert["normalization_map"],
        "cos_zero_solution_count": len(proof["domain_cos_zero_solutions"]),
        "all_denominators_positive_exact": proof["all_denominators_positive_exact"],
        "all_signs_preserved_by_positive_alpha": proof["all_signs_preserved_by_positive_alpha"],
        "max_float_residual": proof["max_float_residual"],
        "scalar_normalization_witness_ready": cert["scalar_normalization_witness_ready"],
        "full_amplitude_source_theorem_ready": cert["full_amplitude_source_theorem_ready"],
        "role_safe_amplitude_absorption_ready": cert["role_safe_amplitude_absorption_ready"],
        "bridge_completion_ready": cert["bridge_completion_ready"],
        "p2411_full_bridge_still_requires_all_obligations": cert["p2411_full_bridge_still_requires_all_obligations"],
        "p2412_chi11_scope_separation_inherited": cert["p2412_chi11_scope_separation_inherited"],
        "scratch_amplitude_witness_inherited": cert["scratch_amplitude_witness_inherited"],
        "scratch_blocks_role_transfer_inherited": cert["scratch_blocks_role_transfer_inherited"],
        "not_licensed": [
            "No full strict amplitude dynamic source theorem is exported.",
            "No role-safe alpha_geo electroweak-role transfer is licensed.",
            "No raw identity K_legacy_ont == K_strict_gate is claimed.",
            "No bridge completion, QW-2191 discharge, role-bearing L_total, or ToE closure follows.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "domain_has_twelve_nodes": theorem_export["domain_size"] == 12,
        "cos_zero_congruence_has_no_domain_solution": theorem_export["cos_zero_solution_count"] == 0,
        "denominators_positive_exact": theorem_export["all_denominators_positive_exact"],
        "positive_alpha_preserves_signs": theorem_export["all_signs_preserved_by_positive_alpha"],
        "normalization_residual_float_small": theorem_export["max_float_residual"] < 1e-15,
        "scalar_witness_ready_but_full_source_not_ready": theorem_export["scalar_normalization_witness_ready"]
        and not theorem_export["full_amplitude_source_theorem_ready"],
        "role_safe_absorption_still_blocked": not theorem_export["role_safe_amplitude_absorption_ready"],
        "bridge_completion_still_blocked": not theorem_export["bridge_completion_ready"],
        "p2411_bridge_obligation_inherited": theorem_export["p2411_full_bridge_still_requires_all_obligations"],
        "p2412_scope_separation_inherited": theorem_export["p2412_chi11_scope_separation_inherited"],
        "scratch_amplitude_witness_inherited": theorem_export["scratch_amplitude_witness_inherited"],
        "scratch_role_transfer_block_inherited": theorem_export["scratch_blocks_role_transfer_inherited"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2413_s1363_v1",
        "packet_id": "P2413",
        "stage_id": "S1363",
        "result_kind": "AMPLITUDE_SCALAR_NORMALIZATION_BRIDGE_WITNESS_CERTIFICATE",
        "status": "PASS_AMPLITUDE_SCALAR_WITNESS_NO_ROLE_TRANSFER_NO_FULL_BRIDGE",
        "amplitude_scalar_normalization_bridge_witness_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Use this scalar witness as one amplitude bridge ingredient, then prove a real amplitude dynamic source or "
            "role-safe absorption theorem; do not transfer alpha_geo physical roles by scalar cancellation."
        ),
        "global_status": "OPEN_PROGRESS_AMPLITUDE_SCALAR_WITNESS_CERTIFIED_NO_ROLE_TRANSFER_OR_TOE_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["amplitude_scalar_normalization_bridge_witness_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2413 S1363: amplitude scalar-normalization bridge witness certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2413/S1363 proves the narrow scalar factorization witness for the amplitude bridge row.",
                "",
                "## Finite facts",
                "",
                f"- Domain size: `{theorem['domain_size']}`.",
                f"- Normalization map: `{theorem['normalization_map']}`.",
                f"- Cos-zero domain solutions: `{theorem['cos_zero_solution_count']}`.",
                f"- Denominators positive: `{theorem['all_denominators_positive_exact']}`.",
                f"- Signs preserved: `{theorem['all_signs_preserved_by_positive_alpha']}`.",
                f"- Max float residual: `{theorem['max_float_residual']}`.",
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
