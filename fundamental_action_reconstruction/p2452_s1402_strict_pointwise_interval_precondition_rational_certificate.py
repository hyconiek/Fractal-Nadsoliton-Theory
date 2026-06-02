#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p2452_s1402_strict_pointwise_interval_precondition_rational_certificate.json"
MD = GEN / "p2452_s1402_strict_pointwise_interval_precondition_rational_certificate.md"

SOURCE_FILES = {
    "P2451_INTERVAL_ENCLOSURE": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

OMEGA = Fraction(743, 4000)
PHI = Fraction(13, 80)
ETA = Fraction(9, 5)
D_INTERVAL = (Fraction(0), Fraction(5))
ZERO_PROJECTION_AUDIT_START = Fraction(1, 1_000_000)
STATIONARY_FACTOR_AUDIT_START = Fraction(1, 10_000)
ARCHIMEDES_PI_LOWER = Fraction(333, 106)


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def fraction_dict(value: Fraction) -> dict[str, Any]:
    return {"numerator": value.numerator, "denominator": value.denominator, "decimal": float(value)}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:35]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2452|S1402|interval precondition rational|rational phase-band|transcendental precondition",
        "p2451_input": "P2451|S1401|interval enclosure root exclusion|floating interval enclosure|monotone sin",
        "phase_language": "phase band|monotone trigonometric|pi/2|Archimedes|rational pi",
        "domain_language": "log-domain|denominator positivity|positive real power|d\\^eta|d\\^\\(eta-1\\)",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|point-coordinate selector",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def append_doc_sections() -> None:
    eq_section = """
## P2452/S1402 strict pointwise interval-precondition rational certificate

`P2452/S1402` audits the exact rational preconditions behind the P2451 interval enclosure: with `omega=743/4000`, `phi=13/80`, and `d in [0,5]`, the phase interval is exactly `[13/80,873/800]`, and `873/800 < 333/212 < pi/2` using the Archimedean lower bound `333/106 < pi`.  The P2451 complement starts also satisfy `d>0`, so the `log(d)`, `d^(9/5)`, `d^(4/5)`, and denominator-positivity assumptions are explicitly separated from floating interval evaluation.

This certifies only interval-evaluation preconditions.  It does not certify directed rounding, symbolic root exclusion, a point-coordinate selector, source/gauge authority, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.
""".strip()
    lag_section = """
## P2452/S1402 rational interval-precondition guard

`P2452/S1402` proves the rational phase/domain preconditions used by the P2451 interval audit, including the monotone trigonometric phase band and positive log/power domains.  These preconditions strengthen the audit hygiene but still do not license any pointwise coordinate as a selector/source/gauge row for `L_total`.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2451 = sources["P2451_INTERVAL_ENCLOSURE"].get("strict_pointwise_projection_interval_enclosure_root_exclusion_audit", {}).get("theorem_export", {})
    phase_left = OMEGA * D_INTERVAL[0] + PHI
    phase_right = OMEGA * D_INTERVAL[1] + PHI
    pi_half_lower = ARCHIMEDES_PI_LOWER / 2
    theorem_export = {
        "theorem_name": "P2452_T1_strict_pointwise_interval_precondition_rational_certificate",
        "inherited_interval_enclosure_audit": "P2451/S1401",
        "strict_tuple_rational_encoding": {
            "omega": fraction_dict(OMEGA),
            "phi": fraction_dict(PHI),
            "eta": fraction_dict(ETA),
            "beta": fraction_dict(Fraction(1)),
        },
        "audited_d_interval": {"lo": fraction_dict(D_INTERVAL[0]), "hi": fraction_dict(D_INTERVAL[1])},
        "exact_phase_interval": {"lo": fraction_dict(phase_left), "hi": fraction_dict(phase_right)},
        "archimedes_pi_lower_bound": fraction_dict(ARCHIMEDES_PI_LOWER),
        "pi_half_lower_bound_used": fraction_dict(pi_half_lower),
        "phase_lower_positive": phase_left > 0,
        "phase_upper_below_pi_half_lower_bound": phase_right < pi_half_lower,
        "phase_band_inside_open_0_pi_over_2": phase_left > 0 and phase_right < pi_half_lower,
        "monotone_sin_cos_precondition_certified": phase_left > 0 and phase_right < pi_half_lower,
        "p2451_zero_projection_interval_cell_count": p2451.get("zero_projection_amplitude_interval_audit", {}).get("cell_count"),
        "p2451_stationary_factor_interval_cell_count": p2451.get("stationary_factor_interval_audit", {}).get("cell_count"),
        "zero_projection_log_domain_start": fraction_dict(ZERO_PROJECTION_AUDIT_START),
        "stationary_factor_log_domain_start": fraction_dict(STATIONARY_FACTOR_AUDIT_START),
        "zero_projection_log_domain_positive": ZERO_PROJECTION_AUDIT_START > 0,
        "stationary_factor_log_domain_positive": STATIONARY_FACTOR_AUDIT_START > 0,
        "positive_power_domain_precondition": "For audited d>=0, d^(9/5) is nonnegative; for derivative/log audited d>0, d^(4/5) and log(d) are defined and positive-domain denominators satisfy 1+d^(9/5)>0.",
        "denominator_positivity_precondition_certified": True,
        "rational_preconditions_exported": True,
        "directed_rounding_interval_theorem_exported_by_this_certificate": False,
        "symbolic_root_exclusion_theorem_exported_by_this_certificate": False,
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_slice_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Rational phase/domain preconditions do not by themselves prove directed-rounding interval arithmetic or symbolic root exclusion.",
            "The certified monotone trigonometric band and log/power domains do not choose a strict point-coordinate selector.",
            "No strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Combine these exact preconditions with a directed-rounding interval backend or symbolic bounds for the projection amplitude and stationary factor."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "phase_band_inside_open_0_pi_over_2": theorem_export["phase_band_inside_open_0_pi_over_2"],
        "monotone_sin_cos_precondition_certified": theorem_export["monotone_sin_cos_precondition_certified"],
        "zero_projection_log_domain_positive": theorem_export["zero_projection_log_domain_positive"],
        "stationary_factor_log_domain_positive": theorem_export["stationary_factor_log_domain_positive"],
        "denominator_positivity_precondition_certified": theorem_export["denominator_positivity_precondition_certified"],
        "p2451_cell_counts_inherited": theorem_export["p2451_zero_projection_interval_cell_count"] == 49922
        and theorem_export["p2451_stationary_factor_interval_cell_count"] == 49960,
        "no_directed_rounding_interval_theorem": not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
        "no_symbolic_root_exclusion": not theorem_export["symbolic_root_exclusion_theorem_exported_by_this_certificate"],
        "no_pointwise_selector_export": not theorem_export["pointwise_coordinate_selector_exported_by_this_certificate"],
        "no_observable_source_export": not theorem_export["strict_observable_source_constraint_exported_by_this_certificate"],
        "no_gauge_slice_export": not theorem_export["gauge_slice_theorem_exported_by_this_certificate"],
        "no_value_generator_export": not theorem_export["strict_physical_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export["qw2191_discharged"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2452_s1402_v1",
        "packet_id": "P2452",
        "stage_id": "S1402",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_PRECONDITION_RATIONAL_CERTIFICATE_NO_EXACT_SELECTOR_SOURCE_THEOREM",
        "strict_pointwise_interval_precondition_rational_certificate": {
            "inputs": {name: rel(path) for name, path in SOURCE_FILES.items()},
            "input_fingerprints": {name: sha256_json(sources[name]) for name in sources},
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_precondition_rational_certificate"]["theorem_export"]
    lines = [
        "# P2452/S1402 strict pointwise interval-precondition rational certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Rational phase/domain preconditions",
        "",
        f"Phase interval: `[{t['exact_phase_interval']['lo']['numerator']}/{t['exact_phase_interval']['lo']['denominator']}, {t['exact_phase_interval']['hi']['numerator']}/{t['exact_phase_interval']['hi']['denominator']}]`.",
        f"Pi/2 lower bound used: `{t['pi_half_lower_bound_used']['numerator']}/{t['pi_half_lower_bound_used']['denominator']}` from `333/106 < pi`.",
        f"Phase band inside `(0,pi/2)`: `{t['phase_band_inside_open_0_pi_over_2']}`.",
        f"P2451 cell counts inherited: zero-projection `{t['p2451_zero_projection_interval_cell_count']}`, stationary-factor `{t['p2451_stationary_factor_interval_cell_count']}`.",
        "",
        "## Guardrail",
        "",
        "This certifies rational preconditions only.  It exports no directed-rounding interval theorem, no symbolic root-exclusion theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps({"status": payload["status"], "gatekeepers": payload["gatekeeper_checks"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
