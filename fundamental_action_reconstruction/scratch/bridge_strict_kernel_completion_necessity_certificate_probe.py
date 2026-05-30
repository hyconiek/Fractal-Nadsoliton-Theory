#!/usr/bin/env python3
"""Scratch probe: necessity certificate for the legacy-carrier -> strict completion factors.

The completion ladder already shows how to multiply the legacy carrier into the
current strict kernel.  This probe asks a narrower proof/computation question:
which parts of that completion are actually forced on the finite Z12 audit
window?

On d=0..11 it enumerates all subsets of the three explicit completion factors

    A(d) = 1/alpha_geo
    P(d) = cos(omega_S*d+phi_S)/cos(omega_L*d+phi_L)
    D(d) = (1+beta_tors*d)/(1+beta*d^eta)

and checks the residual against K_strict_gate.  It also computes the best single
post-hoc scalar for every incomplete subset.  This distinguishes the global
normalization role of A from the genuinely d-dependent shape roles of P and D.

This is still a finite-domain completion certificate, not a derivation of the
factors from nadsoliton dynamics, not a beta_tors -> chi_11 theorem, and not a
role-transfer theorem.
"""
from __future__ import annotations

import itertools
import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_kernel_completion_necessity_certificate_report.json"
OUT_MD = HERE / "bridge_strict_kernel_completion_necessity_certificate_report.md"
LADDER_REPORT = HERE / "bridge_strict_kernel_completion_ladder_report.json"
BLOCKER_LATTICE = HERE / "bridge_completed_kernel_blocker_dependency_lattice_report.json"
INVERSE_BRANCH = HERE / "bridge_phase_normalized_z12_inverse_branch_interval_certificate_report.json"

ALPHA_GEO = 4.0 * math.log(2.0)
LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 9.0 / 5.0}
DOMAIN = list(range(12))
FACTOR_NAMES = ("alpha_normalization", "phase_frequency_transport", "damping_compression")
TOL = 1e-14


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def cos_legacy(d: float) -> float:
    return math.cos(LEGACY["omega"] * d + LEGACY["phi"])


def cos_strict(d: float) -> float:
    return math.cos(STRICT["omega"] * d + STRICT["phi"])


def legacy_full(d: float) -> float:
    return ALPHA_GEO * cos_legacy(d) / (1.0 + LEGACY["beta_tors"] * d)


def strict_kernel(d: float) -> float:
    return cos_strict(d) / (1.0 + STRICT["beta"] * d ** STRICT["eta"])


def factor_value(name: str, d: float) -> float:
    if name == "alpha_normalization":
        return 1.0 / ALPHA_GEO
    if name == "phase_frequency_transport":
        return cos_strict(d) / cos_legacy(d)
    if name == "damping_compression":
        return (1.0 + LEGACY["beta_tors"] * d) / (1.0 + STRICT["beta"] * d ** STRICT["eta"])
    raise ValueError(f"unknown factor: {name}")


def sign(value: float) -> int:
    if value > 0.0:
        return 1
    if value < 0.0:
        return -1
    return 0


def l2(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def subset_label(subset: tuple[str, ...]) -> str:
    return "+".join(subset) if subset else "none"


def apply_subset(d: int, subset: tuple[str, ...]) -> float:
    value = legacy_full(float(d))
    for name in subset:
        value *= factor_value(name, float(d))
    return value


def best_scalar_fit(values: list[float], targets: list[float]) -> dict[str, float]:
    denominator = sum(value * value for value in values)
    scalar = sum(value * target for value, target in zip(values, targets)) / denominator
    residuals = [scalar * value - target for value, target in zip(values, targets)]
    return {
        "best_scalar": scalar,
        "max_abs_residual_after_best_scalar": max(abs(value) for value in residuals),
        "l2_residual_after_best_scalar": l2(residuals),
    }


def subset_certificate(subset: tuple[str, ...], targets: list[float]) -> dict[str, Any]:
    values = [apply_subset(d, subset) for d in DOMAIN]
    residuals = [value - target for value, target in zip(values, targets)]
    sign_mismatches = [d for d, value, target in zip(DOMAIN, values, targets) if sign(value) != sign(target)]
    scalar_fit = best_scalar_fit(values, targets)
    return {
        "subset": list(subset),
        "subset_label": subset_label(subset),
        "missing_factors": [name for name in FACTOR_NAMES if name not in subset],
        "max_abs_residual": max(abs(value) for value in residuals),
        "l1_residual": sum(abs(value) for value in residuals),
        "l2_residual": l2(residuals),
        "sign_mismatch_d_values": sign_mismatches,
        "best_scalar_repair": scalar_fit,
        "exact_without_extra_scalar": max(abs(value) for value in residuals) <= TOL,
        "exact_up_to_one_global_scalar": scalar_fit["max_abs_residual_after_best_scalar"] <= TOL,
    }


def row_certificate(d: int) -> dict[str, Any]:
    df = float(d)
    legacy_value = legacy_full(df)
    strict_value = strict_kernel(df)
    alpha_factor = factor_value("alpha_normalization", df)
    phase_factor = factor_value("phase_frequency_transport", df)
    damping_factor = factor_value("damping_compression", df)
    quotient = strict_value / legacy_value
    product = alpha_factor * phase_factor * damping_factor
    return {
        "d": d,
        "legacy_full": legacy_value,
        "strict_kernel": strict_value,
        "strict_over_legacy_quotient": quotient,
        "alpha_normalization_factor": alpha_factor,
        "phase_frequency_transport_factor": phase_factor,
        "damping_compression_factor": damping_factor,
        "factor_product": product,
        "quotient_minus_factor_product": quotient - product,
        "legacy_sign": sign(legacy_value),
        "strict_sign": sign(strict_value),
        "phase_factor_sign": sign(phase_factor),
        "tail_amplification_if_damping_omitted": 1.0 / damping_factor,
    }


def build_payload() -> dict[str, Any]:
    ladder = load_json(LADDER_REPORT)
    blocker_lattice = load_json(BLOCKER_LATTICE)
    inverse_branch = load_json(INVERSE_BRANCH)

    targets = [strict_kernel(float(d)) for d in DOMAIN]
    subsets = [tuple(combo) for size in range(len(FACTOR_NAMES) + 1) for combo in itertools.combinations(FACTOR_NAMES, size)]
    subset_rows = [subset_certificate(subset, targets) for subset in subsets]
    point_rows = [row_certificate(d) for d in DOMAIN]

    exact_subsets = [row["subset_label"] for row in subset_rows if row["exact_without_extra_scalar"]]
    scalar_exact_subsets = [row["subset_label"] for row in subset_rows if row["exact_up_to_one_global_scalar"]]
    shape_exact_subsets = [label for label in scalar_exact_subsets if "phase_frequency_transport" in label and "damping_compression" in label]

    missing_phase = [row for row in subset_rows if "phase_frequency_transport" in row["missing_factors"]]
    missing_damping = [row for row in subset_rows if "damping_compression" in row["missing_factors"]]

    quotient_residuals = [row["quotient_minus_factor_product"] for row in point_rows]
    damping_values = [row["damping_compression_factor"] for row in point_rows]

    return {
        "result_kind": "SCRATCH_STRICT_KERNEL_COMPLETION_NECESSITY_CERTIFICATE__FINITE_Z12_NO_DYNAMICAL_DERIVATION",
        "status": "all-three-explicit-factors-are-necessary-for-exact-pointwise-completion; phase-and-damping-are-shape-critical",
        "source_reports": {
            "completion_ladder": str(LADDER_REPORT.relative_to(ROOT)),
            "blocker_lattice": str(BLOCKER_LATTICE.relative_to(ROOT)),
            "inverse_branch_interval_certificate": str(INVERSE_BRANCH.relative_to(ROOT)),
        },
        "constants": {
            "alpha_geo": ALPHA_GEO,
            "legacy": LEGACY,
            "strict": STRICT,
            "domain_d_values": DOMAIN,
            "tolerance": TOL,
        },
        "completion_factor_definitions": {
            "alpha_normalization": "A(d)=1/alpha_geo; global amplitude normalization only",
            "phase_frequency_transport": "P(d)=cos(omega_S*d+phi_S)/cos(omega_L*d+phi_L); sign/phase transport",
            "damping_compression": "D(d)=(1+beta_tors*d)/(1+beta*d^eta); denominator compression beta_tors*d -> beta*d^eta",
            "identity": "K_strict_gate(d)=K_legacy_ont(d)*A(d)*P(d)*D(d) on the audited finite domain",
        },
        "necessity_summary": {
            "exact_subsets_without_extra_scalar": exact_subsets,
            "exact_subsets_up_to_one_global_scalar": scalar_exact_subsets,
            "shape_exact_subsets_up_to_one_global_scalar": shape_exact_subsets,
            "phase_missing_sign_mismatch_union": sorted({d for row in missing_phase for d in row["sign_mismatch_d_values"]}),
            "minimum_best_scalar_l2_residual_when_phase_missing": min(row["best_scalar_repair"]["l2_residual_after_best_scalar"] for row in missing_phase),
            "minimum_best_scalar_l2_residual_when_damping_missing": min(row["best_scalar_repair"]["l2_residual_after_best_scalar"] for row in missing_damping),
            "alpha_is_global_scalar_not_shape_factor": "phase+damping without alpha becomes exact after multiplying by the global scalar 1/alpha_geo",
            "phase_is_sign_shape_critical": "without P(d), signs disagree at d=2,3,4,5,8,9 and no scalar removes the mismatch",
            "damping_is_envelope_shape_critical": "without D(d), the tail envelope is too large; at d=11 the omitted damping amplification is about 68.382268",
            "max_abs_quotient_identity_residual": max(abs(value) for value in quotient_residuals),
            "damping_factor_positive": all(value > 0.0 for value in damping_values),
            "damping_factor_strictly_decreasing_after_d0": all(left > right for left, right in zip(damping_values, damping_values[1:])),
        },
        "subset_enumeration": subset_rows,
        "pointwise_quotient_certificate": point_rows,
        "blocker_context": {
            "ladder_status": ladder["status"],
            "blocker_lattice_status": blocker_lattice["status"],
            "inverse_branch_status": inverse_branch["status"],
            "what_this_refines": "the completion ladder is not just an arbitrary display: the finite Z12 audit forces P and D as d-dependent shape corrections and A as the global normalization correction",
            "what_remains_open": [
                "strict_transport_derivation",
                "global_z12_map_derivation_as_nadsoliton_dynamics",
                "orientation_chi11_source",
                "chi11_uniqueness",
                "reynolds_obstruction_escape",
                "role_transfer_theorem",
            ],
        },
        "proof_certificate": {
            "quotient_step": "For every audited d, K_strict_gate(d)/K_legacy_ont(d) equals A(d)*P(d)*D(d) to machine precision.",
            "subset_exhaustion": "All 2^3 subsets of {A,P,D} were enumerated; only A+P+D is exact without an extra scalar.",
            "scalar_repair_step": "Allowing one post-hoc global scalar repairs only the missing-alpha case; it cannot repair missing phase or missing damping shape.",
            "sign_step": "The phase factor is necessary because it removes six legacy-vs-strict sign mismatches on d=0..11.",
            "envelope_step": "The damping factor is necessary because it supplies a positive strictly decreasing envelope compression from 1 at d=0 to the certified tail value at d=11.",
            "nonduplication": "This is a necessity/exhaustion certificate for the displayed completion factors, not another factorization-only or inverse-branch report.",
            "theoretical_limit": "The certificate proves finite-domain necessity relative to the chosen completion ansatz; it does not derive the ansatz from nadsoliton dynamics.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in solitonic state; this audit compares kernel carriers only.",
            "forbidden_reading": "No separate informational layer below the nadsoliton is introduced.",
        },
        "hard_limits": [
            "K_strict_gate remains the current live/full operational kernel.",
            "No unqualified identity K_legacy_ont == K_strict_gate is claimed.",
            "No proof derives A(d), P(d), or D(d) from strict nadsoliton dynamics.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer to K_strict_gate is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    summary = payload["necessity_summary"]
    lines = [
        "# Strict Kernel Completion Necessity Certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This finite-domain audit exhausts all subsets of the explicit completion factors",
        "`A(d)`, `P(d)`, and `D(d)` that complete the historical legacy carrier into",
        "the current strict working kernel on `d=0..11`.",
        "",
        "## Factor definitions",
        "",
    ]
    for name, definition in payload["completion_factor_definitions"].items():
        lines.append(f"- `{name}`: {definition}")
    lines.extend([
        "",
        "## Necessity summary",
        "",
        f"- Exact subsets without extra scalar: `{summary['exact_subsets_without_extra_scalar']}`",
        f"- Exact subsets up to one global scalar: `{summary['exact_subsets_up_to_one_global_scalar']}`",
        f"- Missing-phase sign mismatch union: `{summary['phase_missing_sign_mismatch_union']}`",
        f"- Minimum best-scalar L2 residual when phase is missing: `{summary['minimum_best_scalar_l2_residual_when_phase_missing']:.12e}`",
        f"- Minimum best-scalar L2 residual when damping is missing: `{summary['minimum_best_scalar_l2_residual_when_damping_missing']:.12e}`",
        f"- Max quotient identity residual: `{summary['max_abs_quotient_identity_residual']:.12e}`",
        f"- Damping factor positive: `{summary['damping_factor_positive']}`",
        f"- Damping factor strictly decreasing after d=0: `{summary['damping_factor_strictly_decreasing_after_d0']}`",
        "",
        "## Subset enumeration",
        "",
        "| subset | missing | max residual | best scalar | max residual after scalar | sign mismatches |",
        "|---|---|---:|---:|---:|---|",
    ])
    for row in payload["subset_enumeration"]:
        lines.append(
            "| {subset} | {missing} | {max_res:.12e} | {scalar:.12e} | {scalar_res:.12e} | {signs} |".format(
                subset=row["subset_label"],
                missing=",".join(row["missing_factors"]) or "none",
                max_res=row["max_abs_residual"],
                scalar=row["best_scalar_repair"]["best_scalar"],
                scalar_res=row["best_scalar_repair"]["max_abs_residual_after_best_scalar"],
                signs=row["sign_mismatch_d_values"],
            )
        )
    lines.extend([
        "",
        "## Guarded interpretation",
        "",
        "The audit shows which explicit completion factors are necessary for the finite",
        "Z12 completion identity.  It does not prove that strict nadsoliton dynamics must",
        "generate those factors, does not prove `beta_tors -> chi_11`, and does not",
        "transfer legacy physical roles onto `K_strict_gate`.",
        "",
        "## Hard limits",
        "",
    ])
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
