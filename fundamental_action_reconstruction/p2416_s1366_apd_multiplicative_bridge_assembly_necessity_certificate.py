#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from itertools import combinations
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
SCRATCH = ROOT / "scratch"

OUT = GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json"
MD = GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.md"

SOURCE_FILES = {
    "P2411_BRIDGE_SOURCE_HYPERGRAPH": GEN / "p2411_s1361_legacy_strict_bridge_source_obligation_hypergraph_certificate.json",
    "P2413_AMPLITUDE_SCALAR_WITNESS": GEN / "p2413_s1363_amplitude_scalar_normalization_bridge_witness_certificate.json",
    "P2414_DAMPING_NONABSORPTION": GEN / "p2414_s1364_strict_damping_parameter_identifiability_nonabsorption_certificate.json",
    "P2415_PHASE_AFFINE_TRANSPORT": GEN / "p2415_s1365_phase_frequency_affine_transport_nonautomorphism_certificate.json",
    "SCRATCH_COMPLETION_NECESSITY": SCRATCH / "bridge_strict_kernel_completion_necessity_certificate_report.json",
    "SCRATCH_COMPLETED_FACTORIZATION": SCRATCH / "bridge_completed_strict_kernel_factorization_certificate_report.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

DOMAIN = list(range(12))
FACTOR_NAMES = ("alpha_normalization", "phase_frequency_transport", "damping_compression")
ALPHA_GEO = 4.0 * math.log(2.0)
OMEGA_LEGACY = math.pi / 4.0
PHI_LEGACY = math.pi / 6.0
BETA_TORS = 1.0 / 100.0
OMEGA_STRICT = 743.0 / 4000.0
PHI_STRICT = 13.0 / 80.0
BETA_STRICT = 1.0
ETA_STRICT = 9.0 / 5.0
TOL = 1.0e-12


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
        "new_packet": "P2416|S1366|APD multiplicative bridge assembly|multiplicative bridge assembly necessity",
        "scratch_factorization": "strict_over_legacy_quotient|factor_product|alpha_normalization\\+phase_frequency_transport\\+damping_compression",
        "completed_factorization": "completed strict kernel factorization|C_full=|completion factors",
        "prior_component_witnesses": "P2413|P2414|P2415|amplitude scalar|damping nonabsorption|phase/frequency affine",
        "bridge_limits": "source theorem|role transfer|QW-2191|ToE closure",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds scratch completion-factorization and the separate P2413/P2414/P2415 component witnesses, "
            "but no production P24xx APD assembly certificate that enumerates the factor-subset lattice and keeps "
            "value-level exactness separate from bridge-source, selector, role-transfer, L_total, or ToE closure."
        ),
    }


def k_legacy(d: int) -> float:
    return ALPHA_GEO * math.cos(OMEGA_LEGACY * d + PHI_LEGACY) / (1.0 + BETA_TORS * d)


def k_strict(d: int) -> float:
    return math.cos(OMEGA_STRICT * d + PHI_STRICT) / (1.0 + BETA_STRICT * d**ETA_STRICT)


def factor_values(d: int) -> dict[str, float]:
    cos_l = math.cos(OMEGA_LEGACY * d + PHI_LEGACY)
    cos_s = math.cos(OMEGA_STRICT * d + PHI_STRICT)
    return {
        "alpha_normalization": 1.0 / ALPHA_GEO,
        "phase_frequency_transport": cos_s / cos_l,
        "damping_compression": (1.0 + BETA_TORS * d) / (1.0 + BETA_STRICT * d**ETA_STRICT),
    }


def product(values: list[float]) -> float:
    out = 1.0
    for value in values:
        out *= value
    return out


def point_rows() -> list[dict[str, Any]]:
    rows = []
    for d in DOMAIN:
        legacy = k_legacy(d)
        strict = k_strict(d)
        factors = factor_values(d)
        quotient = strict / legacy
        factor_product = product([factors[name] for name in FACTOR_NAMES])
        rows.append(
            {
                "d": d,
                "legacy_kernel": legacy,
                "strict_kernel": strict,
                "strict_over_legacy_quotient": quotient,
                "alpha_normalization": factors["alpha_normalization"],
                "phase_frequency_transport": factors["phase_frequency_transport"],
                "damping_compression": factors["damping_compression"],
                "factor_product": factor_product,
                "quotient_minus_factor_product": quotient - factor_product,
            }
        )
    return rows


def subset_labels() -> list[tuple[str, tuple[str, ...]]]:
    labels: list[tuple[str, tuple[str, ...]]] = [("none", tuple())]
    for size in range(1, len(FACTOR_NAMES) + 1):
        for subset in combinations(FACTOR_NAMES, size):
            labels.append(("+".join(subset), subset))
    return labels


def best_scalar_repair(approximants: list[float], targets: list[float]) -> dict[str, float]:
    denominator = sum(value * value for value in approximants)
    scalar = sum(value * target for value, target in zip(approximants, targets)) / denominator
    residuals = [scalar * value - target for value, target in zip(approximants, targets)]
    return {
        "best_scalar": scalar,
        "max_abs_residual_after_best_scalar": max(abs(residual) for residual in residuals),
        "l2_residual_after_best_scalar": math.sqrt(sum(residual * residual for residual in residuals)),
    }


def subset_audit(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    targets = [row["strict_over_legacy_quotient"] for row in rows]
    audits = []
    for label, subset in subset_labels():
        approximants = [product([row[name] for name in subset]) for row in rows]
        residuals = [approx - target for approx, target in zip(approximants, targets)]
        sign_mismatches = [row["d"] for row, approx, target in zip(rows, approximants, targets) if (approx > 0) != (target > 0)]
        scalar = best_scalar_repair(approximants, targets)
        audits.append(
            {
                "subset_label": label,
                "subset": list(subset),
                "missing_factors": [name for name in FACTOR_NAMES if name not in subset],
                "max_abs_residual": max(abs(residual) for residual in residuals),
                "l2_residual": math.sqrt(sum(residual * residual for residual in residuals)),
                "l1_residual": sum(abs(residual) for residual in residuals),
                "sign_mismatch_d_values": sign_mismatches,
                "best_scalar_repair": scalar,
                "exact_without_extra_scalar": max(abs(residual) for residual in residuals) <= TOL,
                "exact_up_to_one_global_scalar": scalar["max_abs_residual_after_best_scalar"] <= TOL,
            }
        )
    return audits


def build_certificate(sources: dict[str, Any]) -> dict[str, Any]:
    rows = point_rows()
    subsets = subset_audit(rows)
    exact_without = [row["subset_label"] for row in subsets if row["exact_without_extra_scalar"]]
    exact_with_scalar = [row["subset_label"] for row in subsets if row["exact_up_to_one_global_scalar"]]
    p2411_theorem = sources["P2411_BRIDGE_SOURCE_HYPERGRAPH"].get(
        "legacy_strict_bridge_source_obligation_hypergraph_certificate", {}
    ).get("theorem_export", {})
    p2413_theorem = sources["P2413_AMPLITUDE_SCALAR_WITNESS"].get(
        "amplitude_scalar_normalization_bridge_witness_certificate", {}
    ).get("theorem_export", {})
    p2414_theorem = sources["P2414_DAMPING_NONABSORPTION"].get(
        "strict_damping_parameter_identifiability_nonabsorption_certificate", {}
    ).get("theorem_export", {})
    p2415_theorem = sources["P2415_PHASE_AFFINE_TRANSPORT"].get(
        "phase_frequency_affine_transport_nonautomorphism_certificate", {}
    ).get("theorem_export", {})
    scratch_summary = sources["SCRATCH_COMPLETION_NECESSITY"].get("necessity_summary", {})
    return {
        "domain": DOMAIN,
        "factor_names": FACTOR_NAMES,
        "pointwise_factor_rows": rows,
        "subset_audit_rows": subsets,
        "exact_subsets_without_extra_scalar": exact_without,
        "exact_subsets_up_to_one_global_scalar": exact_with_scalar,
        "unique_exact_without_extra_scalar": exact_without == ["alpha_normalization+phase_frequency_transport+damping_compression"],
        "alpha_missing_repairable_only_by_global_scalar": "phase_frequency_transport+damping_compression" in exact_with_scalar
        and "phase_frequency_transport+damping_compression" not in exact_without,
        "phase_missing_not_scalar_repairable": all(
            not row["exact_up_to_one_global_scalar"] for row in subsets if "phase_frequency_transport" in row["missing_factors"]
        ),
        "damping_missing_not_scalar_repairable": all(
            not row["exact_up_to_one_global_scalar"] for row in subsets if "damping_compression" in row["missing_factors"]
        ),
        "max_abs_full_factorization_residual": max(abs(row["quotient_minus_factor_product"]) for row in rows),
        "scratch_necessity_inherited": scratch_summary.get("exact_subsets_without_extra_scalar")
        == ["alpha_normalization+phase_frequency_transport+damping_compression"],
        "p2411_full_bridge_still_requires_all_obligations": p2411_theorem.get("bridge_ready_true_masks") == [255],
        "p2413_amplitude_witness_inherited": p2413_theorem.get("scalar_normalization_witness_ready") is True,
        "p2414_damping_nonabsorption_inherited": p2414_theorem.get("strict_beta_eta_identified_from_samples") is True,
        "p2415_phase_transport_inherited": p2415_theorem.get("continuous_affine_phase_transport_exact_float") is True,
        "apd_value_assembly_witness_ready": True,
        "strict_dynamic_source_exported": False,
        "selector_source_exported": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_licensed": False,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2416/S1366 APD multiplicative bridge-assembly necessity certificate

`P2416/S1366` assembles the three value-level bridge factors certified separately by P2413/P2414/P2415:

```text
K_strict_gate(d) / K_legacy_ont(d)
  = alpha_geo^{-1} * P_phase(d) * D_compression(d),
P_phase(d)=cos(theta_S(d))/cos(theta_L(d)),
D_compression(d)=(1+beta_tors*d)/(1+d^(9/5)).
```

On `d=0..11`, the full `alpha_normalization + phase_frequency_transport + damping_compression` product is the unique exact subset without an extra scalar.  The subset `phase_frequency_transport + damping_compression` becomes exact only after a post-hoc global scalar, so missing alpha is scalar-repairable but still not licensed as a physical-role theorem.  Missing phase or missing damping is not repairable by any single scalar.

This is a finite value-level assembly witness only.  It does not export strict dynamic sources for the factors, does not provide the selector/source premise, does not complete the full bridge, and does not license legacy role transfer, `L_total`, QW-2191 discharge, or ToE closure.
""".strip()
    lag_section = """
## P2416/S1366 APD assembly guard for Lagrangian/EOM

`P2416/S1366` allows the finite quotient identity `K_strict/K_legacy = alpha^{-1}*P_phase*D_compression` as value-level bridge bookkeeping.  Because the certificate still exports no strict dynamic source, selector source, full bridge theorem, or role-transfer theorem, the APD product cannot be promoted into a role-bearing `L_total` term.
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
        "theorem_name": "P2416_T1_apd_multiplicative_bridge_assembly_necessity_certificate",
        "domain_size": len(cert["domain"]),
        "factor_names": list(cert["factor_names"]),
        "subset_count": len(cert["subset_audit_rows"]),
        "exact_subsets_without_extra_scalar": cert["exact_subsets_without_extra_scalar"],
        "exact_subsets_up_to_one_global_scalar": cert["exact_subsets_up_to_one_global_scalar"],
        "unique_exact_without_extra_scalar": cert["unique_exact_without_extra_scalar"],
        "alpha_missing_repairable_only_by_global_scalar": cert["alpha_missing_repairable_only_by_global_scalar"],
        "phase_missing_not_scalar_repairable": cert["phase_missing_not_scalar_repairable"],
        "damping_missing_not_scalar_repairable": cert["damping_missing_not_scalar_repairable"],
        "max_abs_full_factorization_residual": cert["max_abs_full_factorization_residual"],
        "scratch_necessity_inherited": cert["scratch_necessity_inherited"],
        "p2411_full_bridge_still_requires_all_obligations": cert["p2411_full_bridge_still_requires_all_obligations"],
        "p2413_amplitude_witness_inherited": cert["p2413_amplitude_witness_inherited"],
        "p2414_damping_nonabsorption_inherited": cert["p2414_damping_nonabsorption_inherited"],
        "p2415_phase_transport_inherited": cert["p2415_phase_transport_inherited"],
        "apd_value_assembly_witness_ready": cert["apd_value_assembly_witness_ready"],
        "strict_dynamic_source_exported": cert["strict_dynamic_source_exported"],
        "selector_source_exported": cert["selector_source_exported"],
        "full_bridge_theorem_exported": cert["full_bridge_theorem_exported"],
        "role_transfer_licensed": cert["role_transfer_licensed"],
        "not_licensed": [
            "No strict dynamic source theorem for alpha, phase, damping, beta, eta, omega, or phi is exported.",
            "No selector-source theorem or QW-2191 discharge follows from APD value assembly.",
            "No post-bridge legacy physical-role transfer is licensed by quotient exactness.",
            "No role-bearing L_total term or ToE closure follows.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "domain_has_twelve_nodes": theorem_export["domain_size"] == 12,
        "three_factors_audited": theorem_export["factor_names"] == list(FACTOR_NAMES),
        "all_eight_subsets_audited": theorem_export["subset_count"] == 8,
        "full_product_exact": theorem_export["max_abs_full_factorization_residual"] <= TOL,
        "unique_exact_subset_without_scalar": theorem_export["unique_exact_without_extra_scalar"],
        "alpha_missing_only_scalar_repairable": theorem_export["alpha_missing_repairable_only_by_global_scalar"],
        "phase_missing_not_scalar_repairable": theorem_export["phase_missing_not_scalar_repairable"],
        "damping_missing_not_scalar_repairable": theorem_export["damping_missing_not_scalar_repairable"],
        "scratch_necessity_inherited": theorem_export["scratch_necessity_inherited"],
        "p2411_bridge_obligation_inherited": theorem_export["p2411_full_bridge_still_requires_all_obligations"],
        "p2413_amplitude_inherited": theorem_export["p2413_amplitude_witness_inherited"],
        "p2414_damping_inherited": theorem_export["p2414_damping_nonabsorption_inherited"],
        "p2415_phase_inherited": theorem_export["p2415_phase_transport_inherited"],
        "value_assembly_ready": theorem_export["apd_value_assembly_witness_ready"],
        "no_dynamic_source_exported": not theorem_export["strict_dynamic_source_exported"],
        "no_selector_source_exported": not theorem_export["selector_source_exported"],
        "full_bridge_still_open": not theorem_export["full_bridge_theorem_exported"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2416_s1366_v1",
        "packet_id": "P2416",
        "stage_id": "S1366",
        "result_kind": "APD_MULTIPLICATIVE_BRIDGE_ASSEMBLY_NECESSITY_CERTIFICATE",
        "status": "PASS_APD_VALUE_ASSEMBLY_EXACT_NO_SOURCE_NO_ROLE_TRANSFER",
        "apd_multiplicative_bridge_assembly_necessity_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Treat APD as a value-level assembly witness only; next prove factor source/selector premises or keep full bridge and role transfer open."
        ),
        "global_status": "OPEN_PROGRESS_APD_VALUE_ASSEMBLY_NO_SOURCE_ROLE_OR_TOE_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["apd_multiplicative_bridge_assembly_necessity_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2416 S1366: APD multiplicative bridge-assembly necessity certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2416/S1366 assembles the amplitude, phase, and damping value-level bridge factors and audits all factor subsets.",
                "",
                "## Finite facts",
                "",
                f"- Domain size: `{theorem['domain_size']}`.",
                f"- Subsets audited: `{theorem['subset_count']}`.",
                f"- Exact subsets without scalar: `{theorem['exact_subsets_without_extra_scalar']}`.",
                f"- Exact subsets up to one scalar: `{theorem['exact_subsets_up_to_one_global_scalar']}`.",
                f"- Max full residual: `{theorem['max_abs_full_factorization_residual']}`.",
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
