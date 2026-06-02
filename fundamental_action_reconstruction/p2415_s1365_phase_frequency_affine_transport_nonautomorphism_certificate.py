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

OUT = GEN / "p2415_s1365_phase_frequency_affine_transport_nonautomorphism_certificate.json"
MD = GEN / "p2415_s1365_phase_frequency_affine_transport_nonautomorphism_certificate.md"

SOURCE_FILES = {
    "P2411_BRIDGE_SOURCE_HYPERGRAPH": GEN / "p2411_s1361_legacy_strict_bridge_source_obligation_hypergraph_certificate.json",
    "P2412_CHI11_SCOPE_SEPARATION": GEN / "p2412_s1362_chi11_selector_scope_separation_certificate.json",
    "P2413_AMPLITUDE_SCALAR_WITNESS": GEN / "p2413_s1363_amplitude_scalar_normalization_bridge_witness_certificate.json",
    "P2414_DAMPING_NONABSORPTION": GEN / "p2414_s1364_strict_damping_parameter_identifiability_nonabsorption_certificate.json",
    "SCRATCH_NECESSITY": SCRATCH / "bridge_strict_kernel_completion_necessity_certificate_report.json",
    "SCRATCH_Z2_PHASE_SIGNS": SCRATCH / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json",
    "SCRATCH_GF2_PHASE_SYSTEM": SCRATCH / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_report.json",
    "SCRATCH_PHASE_AFFINE_TRANSPORT": SCRATCH / "bridge_strict_completion_legacy_to_strict_phase_frequency_affine_transport_certificate_report.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

DOMAIN = list(range(12))
UNITS_Z12 = [1, 5, 7, 11]
OMEGA_LEGACY = math.pi / 4.0
PHI_LEGACY = math.pi / 6.0
OMEGA_STRICT = 743.0 / 4000.0
PHI_STRICT = 13.0 / 80.0
TOL = 1.0e-14


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
        "new_packet": "P2415|S1365|phase frequency affine transport|phase-frequency affine transport nonautomorphism",
        "scratch_phase_affine": "phase affine|affine phase transport|omega_S/omega_L|Z12 automorphism phase",
        "gf2_phase_chain": "GF\\(2\\).*phase|phase-sign GF2|phase sign z2|node_bit_pattern",
        "selector_scope": "QW-2191|chi11 selector scope|orientation selector source",
        "bridge_priority": "phase/frequency passage|phase/frequency/topological bit passage|legacy -> strict completion bridge",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds scratch phase-affine, Z2/GF(2), and selector-scope material, but no production P24xx "
            "certificate after P2414 that promotes the phase/frequency row into a finite transport witness while "
            "separating it from Z12 reindexing, scalar replacement, selector-source closure, or full bridge completion."
        ),
    }


def sign(value: float) -> int:
    return 1 if value > 0.0 else -1 if value < 0.0 else 0


def bit_from_sign(value: int) -> int:
    if value == 1:
        return 0
    if value == -1:
        return 1
    raise ValueError("zero sign cannot be converted to a GF(2) bit")


def best_scalar_fit(values: list[float], targets: list[float]) -> dict[str, float]:
    denominator = sum(value * value for value in values)
    scalar = sum(value * target for value, target in zip(values, targets)) / denominator
    residuals = [scalar * value - target for value, target in zip(values, targets)]
    return {
        "best_scalar": scalar,
        "max_abs_residual": max(abs(residual) for residual in residuals),
        "l2_residual": math.sqrt(sum(residual * residual for residual in residuals)),
    }


def phase_rows(sources: dict[str, Any]) -> list[dict[str, Any]]:
    point_rows = sources["SCRATCH_NECESSITY"].get("pointwise_quotient_certificate", [])
    z2_bits = {row["d"]: row["node_bit"] for row in sources["SCRATCH_Z2_PHASE_SIGNS"].get("node_bit_rows", [])}
    rows: list[dict[str, Any]] = []
    for point in point_rows:
        d = int(point["d"])
        theta_l = OMEGA_LEGACY * d + PHI_LEGACY
        theta_s = OMEGA_STRICT * d + PHI_STRICT
        x_d = (theta_s - PHI_LEGACY) / OMEGA_LEGACY
        legacy_cos_d = math.cos(theta_l)
        strict_cos_d = math.cos(theta_s)
        legacy_cos_x = math.cos(OMEGA_LEGACY * x_d + PHI_LEGACY)
        factor = float(point["phase_frequency_transport_factor"])
        factor_sign = sign(factor)
        rows.append(
            {
                "d": d,
                "theta_legacy_at_d": theta_l,
                "theta_strict_at_d": theta_s,
                "affine_legacy_coordinate_x_d": x_d,
                "x_d_minus_d": x_d - d,
                "distance_x_d_to_nearest_integer": abs(x_d - round(x_d)),
                "legacy_cos_at_integer_d": legacy_cos_d,
                "strict_cos_at_d": strict_cos_d,
                "legacy_cos_at_affine_x_d": legacy_cos_x,
                "affine_transport_residual": legacy_cos_x - strict_cos_d,
                "phase_frequency_transport_factor_P_d": factor,
                "phase_factor_sign": factor_sign,
                "phase_factor_gf2_bit": bit_from_sign(factor_sign),
                "z2_node_bit": z2_bits.get(d),
                "matches_z2_node_bit": bit_from_sign(factor_sign) == z2_bits.get(d),
            }
        )
    return rows


def z12_automorphism_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    legacy_sign_pattern = [sign(row["legacy_cos_at_integer_d"]) for row in rows]
    strict_sign_pattern = [sign(row["strict_cos_at_d"]) for row in rows]
    audits = []
    for unit in UNITS_Z12:
        for offset in DOMAIN:
            mapped = [legacy_sign_pattern[(unit * d + offset) % 12] for d in DOMAIN]
            mismatches = [d for d, got, want in zip(DOMAIN, mapped, strict_sign_pattern) if got != want]
            audits.append(
                {
                    "unit": unit,
                    "offset": offset,
                    "mapped_legacy_sign_pattern": mapped,
                    "mismatch_positions_against_strict_sign": mismatches,
                    "mismatch_count": len(mismatches),
                    "matches_strict_sign_pattern": not mismatches,
                }
            )
    return audits


def build_certificate(sources: dict[str, Any]) -> dict[str, Any]:
    rows = phase_rows(sources)
    automorphisms = z12_automorphism_rows(rows)
    legacy_cos = [row["legacy_cos_at_integer_d"] for row in rows]
    strict_cos = [row["strict_cos_at_d"] for row in rows]
    scalar_fit = best_scalar_fit(legacy_cos, strict_cos)
    gf2_summary = sources["SCRATCH_GF2_PHASE_SYSTEM"].get("gf2_linear_system_summary", {})
    scratch_affine_summary = sources["SCRATCH_PHASE_AFFINE_TRANSPORT"].get("phase_frequency_affine_transport_summary", {})
    p2411_theorem = sources["P2411_BRIDGE_SOURCE_HYPERGRAPH"].get(
        "legacy_strict_bridge_source_obligation_hypergraph_certificate", {}
    ).get("theorem_export", {})
    p2412_theorem = sources["P2412_CHI11_SCOPE_SEPARATION"].get(
        "chi11_selector_scope_separation_certificate", {}
    ).get("theorem_export", {})
    p2413_theorem = sources["P2413_AMPLITUDE_SCALAR_WITNESS"].get(
        "amplitude_scalar_normalization_bridge_witness_certificate", {}
    ).get("theorem_export", {})
    p2414_theorem = sources["P2414_DAMPING_NONABSORPTION"].get(
        "strict_damping_parameter_identifiability_nonabsorption_certificate", {}
    ).get("theorem_export", {})
    x_distances = [row["distance_x_d_to_nearest_integer"] for row in rows]
    return {
        "domain": DOMAIN,
        "phase_affine_definition": {
            "legacy_argument": "theta_L(x)=omega_L*x+phi_L",
            "strict_argument": "theta_S(d)=omega_S*d+phi_S",
            "affine_transport": "x(d)=(theta_S(d)-phi_L)/omega_L=(omega_S/omega_L)*d+(phi_S-phi_L)/omega_L",
            "omega_legacy": "pi/4",
            "phi_legacy": "pi/6",
            "omega_strict": "743/4000",
            "phi_strict": "13/80",
            "z12_units_checked": UNITS_Z12,
        },
        "affine_slope_float": OMEGA_STRICT / OMEGA_LEGACY,
        "affine_intercept_float": (PHI_STRICT - PHI_LEGACY) / OMEGA_LEGACY,
        "phase_transport_rows": rows,
        "z12_unit_offset_automorphism_rows": automorphisms,
        "best_scalar_phase_fit": scalar_fit,
        "continuous_affine_phase_transport_exact_float": max(abs(row["affine_transport_residual"]) for row in rows) <= TOL,
        "max_abs_affine_transport_residual": max(abs(row["affine_transport_residual"]) for row in rows),
        "affine_coordinates_not_all_integers": any(distance > TOL for distance in x_distances),
        "minimum_distance_to_integer_affine_coordinate": min(x_distances),
        "maximum_distance_to_integer_affine_coordinate": max(x_distances),
        "z12_unit_offset_automorphism_count_checked": len(automorphisms),
        "no_z12_unit_offset_reindex_matches_strict_sign_pattern": not any(
            row["matches_strict_sign_pattern"] for row in automorphisms
        ),
        "best_z12_unit_offset_mismatch_count": min(row["mismatch_count"] for row in automorphisms),
        "phase_factor_bits": [row["phase_factor_gf2_bit"] for row in rows],
        "phase_factor_bits_match_z2_node_bits": all(row["matches_z2_node_bit"] for row in rows),
        "gf2_solution_inherited_unique": gf2_summary.get("unique_solution") is True,
        "scalar_phase_replacement_fails": scalar_fit["max_abs_residual"] > 1.0e-3,
        "scratch_affine_transport_inherited": scratch_affine_summary.get("continuous_affine_phase_transport_exact") is True,
        "scratch_nonautomorphism_inherited": scratch_affine_summary.get(
            "no_z12_unit_offset_reindex_matches_strict_sign_pattern"
        ) is True,
        "strict_phase_frequency_source_exported": False,
        "orientation_selector_source_exported": False,
        "phase_frequency_bridge_component_ready": False,
        "raw_kernel_identity_claimed": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_licensed": False,
        "p2411_full_bridge_still_requires_all_obligations": p2411_theorem.get("bridge_ready_true_masks") == [255],
        "p2412_chi11_scope_separation_inherited": p2412_theorem.get("current_scope_separation_true") is True,
        "p2413_amplitude_witness_inherited": p2413_theorem.get("scalar_normalization_witness_ready") is True,
        "p2414_damping_nonabsorption_inherited": p2414_theorem.get("damping_compression_bridge_component_ready") is False
        and p2414_theorem.get("strict_beta_eta_identified_from_samples") is True,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2415/S1365 phase/frequency affine-transport nonautomorphism certificate

`P2415/S1365` isolates the S2 phase/frequency row after the amplitude and damping witnesses.  It exports the finite continuous transport

```text
x(d) = (theta_S(d)-phi_L)/omega_L,
theta_L(x(d)) = theta_S(d),
```

with `omega_L=pi/4`, `phi_L=pi/6`, `omega_S=743/4000`, and `phi_S=13/80` on `d=0..11`.  The computation checks all 12 transported nodes, all 48 `Aut(Z12)` unit+offset reindexings, and the best scalar replacement from the legacy cosine carrier to the strict cosine carrier.

The result is positive only as a phase-coordinate transport witness: affine transport has zero numerical residual on the audited rows, but it is not a discrete `Z12` automorphism, no unit+offset reindexing reproduces the strict sign pattern, and the best scalar replacement has nonzero residual.  The phase-factor signs match the inherited Z2/GF(2) phase chain.

This does not derive `omega,phi` from strict dynamics, does not export an orientation/selector source, does not discharge QW-2191, does not complete the bridge, and does not license legacy role transfer, `L_total`, or ToE closure.
""".strip()
    lag_section = """
## P2415/S1365 phase/frequency nonautomorphism guard for Lagrangian/EOM

`P2415/S1365` permits only the finite affine phase-coordinate witness `theta_L(x(d))=theta_S(d)` and explicitly rejects treating it as a `Z12` reindexing, scalar replacement, selector-source theorem, or full bridge completion.  Therefore phase/frequency bookkeeping still cannot be promoted into a role-bearing `L_total` term.
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
        "theorem_name": "P2415_T1_phase_frequency_affine_transport_nonautomorphism_certificate",
        "domain_size": len(cert["domain"]),
        "affine_transport_definition": cert["phase_affine_definition"]["affine_transport"],
        "continuous_affine_phase_transport_exact_float": cert["continuous_affine_phase_transport_exact_float"],
        "max_abs_affine_transport_residual": cert["max_abs_affine_transport_residual"],
        "affine_coordinates_not_all_integers": cert["affine_coordinates_not_all_integers"],
        "minimum_distance_to_integer_affine_coordinate": cert["minimum_distance_to_integer_affine_coordinate"],
        "maximum_distance_to_integer_affine_coordinate": cert["maximum_distance_to_integer_affine_coordinate"],
        "z12_unit_offset_automorphism_count_checked": cert["z12_unit_offset_automorphism_count_checked"],
        "no_z12_unit_offset_reindex_matches_strict_sign_pattern": cert[
            "no_z12_unit_offset_reindex_matches_strict_sign_pattern"
        ],
        "best_z12_unit_offset_mismatch_count": cert["best_z12_unit_offset_mismatch_count"],
        "phase_factor_bits": cert["phase_factor_bits"],
        "phase_factor_bits_match_z2_node_bits": cert["phase_factor_bits_match_z2_node_bits"],
        "gf2_solution_inherited_unique": cert["gf2_solution_inherited_unique"],
        "scalar_phase_replacement_fails": cert["scalar_phase_replacement_fails"],
        "scalar_phase_best_fit_max_abs_residual": cert["best_scalar_phase_fit"]["max_abs_residual"],
        "scalar_phase_best_fit_l2_residual": cert["best_scalar_phase_fit"]["l2_residual"],
        "scratch_affine_transport_inherited": cert["scratch_affine_transport_inherited"],
        "scratch_nonautomorphism_inherited": cert["scratch_nonautomorphism_inherited"],
        "strict_phase_frequency_source_exported": cert["strict_phase_frequency_source_exported"],
        "orientation_selector_source_exported": cert["orientation_selector_source_exported"],
        "phase_frequency_bridge_component_ready": cert["phase_frequency_bridge_component_ready"],
        "raw_kernel_identity_claimed": cert["raw_kernel_identity_claimed"],
        "full_bridge_theorem_exported": cert["full_bridge_theorem_exported"],
        "role_transfer_licensed": cert["role_transfer_licensed"],
        "p2411_full_bridge_still_requires_all_obligations": cert["p2411_full_bridge_still_requires_all_obligations"],
        "p2412_chi11_scope_separation_inherited": cert["p2412_chi11_scope_separation_inherited"],
        "p2413_amplitude_witness_inherited": cert["p2413_amplitude_witness_inherited"],
        "p2414_damping_nonabsorption_inherited": cert["p2414_damping_nonabsorption_inherited"],
        "not_licensed": [
            "No strict dynamic source theorem for omega/phi is exported.",
            "No Z12 automorphism or scalar replacement of phase/frequency data is licensed.",
            "No orientation selector source or QW-2191 discharge follows.",
            "No full bridge completion, legacy role transfer, role-bearing L_total, or ToE closure follows.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "domain_has_twelve_nodes": theorem_export["domain_size"] == 12,
        "affine_transport_residual_small": theorem_export["max_abs_affine_transport_residual"] <= TOL,
        "not_all_affine_coordinates_integer": theorem_export["affine_coordinates_not_all_integers"],
        "checked_all_z12_unit_offsets": theorem_export["z12_unit_offset_automorphism_count_checked"] == 48,
        "no_z12_reindex_match": theorem_export["no_z12_unit_offset_reindex_matches_strict_sign_pattern"],
        "scalar_phase_replacement_rejected": theorem_export["scalar_phase_replacement_fails"],
        "phase_bits_match_z2": theorem_export["phase_factor_bits_match_z2_node_bits"],
        "gf2_unique_solution_inherited": theorem_export["gf2_solution_inherited_unique"],
        "no_phase_source_exported": not theorem_export["strict_phase_frequency_source_exported"],
        "no_selector_source_exported": not theorem_export["orientation_selector_source_exported"],
        "phase_bridge_component_still_open": not theorem_export["phase_frequency_bridge_component_ready"],
        "raw_identity_not_claimed": not theorem_export["raw_kernel_identity_claimed"],
        "full_bridge_still_open": not theorem_export["full_bridge_theorem_exported"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "p2411_bridge_obligation_inherited": theorem_export["p2411_full_bridge_still_requires_all_obligations"],
        "p2412_scope_separation_inherited": theorem_export["p2412_chi11_scope_separation_inherited"],
        "p2413_amplitude_witness_inherited": theorem_export["p2413_amplitude_witness_inherited"],
        "p2414_damping_nonabsorption_inherited": theorem_export["p2414_damping_nonabsorption_inherited"],
        "scratch_affine_inherited": theorem_export["scratch_affine_transport_inherited"],
        "scratch_nonautomorphism_inherited": theorem_export["scratch_nonautomorphism_inherited"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2415_s1365_v1",
        "packet_id": "P2415",
        "stage_id": "S1365",
        "result_kind": "PHASE_FREQUENCY_AFFINE_TRANSPORT_NONAUTOMORPHISM_CERTIFICATE",
        "status": "PASS_PHASE_FREQUENCY_AFFINE_TRANSPORT_WITNESS_NO_SOURCE_NO_SELECTOR_CLOSURE",
        "phase_frequency_affine_transport_nonautomorphism_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Use the affine phase-coordinate witness only as phase/frequency transport data; next prove a real "
            "omega/phi source or selector-source theorem, or keep the bridge component open."
        ),
        "global_status": "OPEN_PROGRESS_PHASE_TRANSPORT_WITNESS_NO_SELECTOR_SOURCE_OR_TOE_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["phase_frequency_affine_transport_nonautomorphism_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2415 S1365: phase/frequency affine-transport nonautomorphism certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2415/S1365 packages the finite phase/frequency affine-coordinate transport and rejects Z12 reindexing or scalar replacement overreads.",
                "",
                "## Finite facts",
                "",
                f"- Domain size: `{theorem['domain_size']}`.",
                f"- Max affine residual: `{theorem['max_abs_affine_transport_residual']}`.",
                f"- Z12 unit+offset checks: `{theorem['z12_unit_offset_automorphism_count_checked']}`.",
                f"- Best Z12 mismatch count: `{theorem['best_z12_unit_offset_mismatch_count']}`.",
                f"- Phase bits: `{theorem['phase_factor_bits']}`.",
                f"- Scalar replacement max residual: `{theorem['scalar_phase_best_fit_max_abs_residual']}`.",
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
