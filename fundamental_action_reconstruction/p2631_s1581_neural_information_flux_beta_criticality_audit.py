#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any, Callable

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2631_s1581_neural_information_flux_beta_criticality_audit.json"
MD = GEN / "p2631_s1581_neural_information_flux_beta_criticality_audit.md"

REPO_ROOT = REPO
P2629_JSON = GEN / "p2629_s1579_zbeta_normalization_gauge_obstruction.json"
P2630_JSON = GEN / "p2630_s1580_strict_beta_source_vs_bridge_zbeta_separation.json"

DF = 9.0 / 5.0
ETA = 9.0 / 5.0
OMEGA = math.pi / 4.0
PHI = math.pi / 6.0
DOMAIN_DMAX = 10.0
DOMAIN_DMIN = 1.0e-6
BETA_STRICT = 1.0

NEGATIVE_EXPORT_FLAGS = [
    "information_conservation_beta_identity_exported",
    "edge_of_chaos_beta_identity_exported",
    "uv_gauge_obstruction_closed",
    "positive_beta_renormalization_source_exported",
    "nonlinear_damping_completion_source_exported",
    "full_legacy_to_strict_bridge_revalidated",
    "role_transfer_revalidated",
    "role_bearing_ltotal_reenabled",
    "qw2191_discharged",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.lean", "-g", "*.json",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO_ROOT,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def semantic_rg_audit() -> dict[str, Any]:
    # Content-first anti-duplication: neural/information-critical content, not ticket IDs.
    patterns = {
        "neural_kernel_prism_content": (
            "positional encoding|heavy-tailed attention|attention bias|ALiBi|Energy-Based|"
            "Boltzmann|Gaussian Process|fractal graph|network temperature"
        ),
        "information_flux_norm_content": (
            "information flux|zachowanie strumienia|flux conservation|unitarity|norm conservation|"
            "fractal measure|dmu_f|measure.*fractal|pure information"
        ),
        "criticality_edge_chaos_content": (
            "edge of chaos|krawedz chaosu|krawędź chaosu|critical point|phase transition|"
            "frozen network|overheated|mutual information|przepustowosc|przepustowość"
        ),
        "uv_normalization_invariance_content": (
            "UV normalization|beta_uv|normalization gauge|normalization-invariant|scale convention|"
            "target-independent.*normalization|beta_micro/beta_strict"
        ),
        "bridge_closure_guard_content": (
            "positive_beta_renormalization_source|nonlinear_damping_completion_source|"
            "role-bearing L_total|role transfer|QW-2191|ToE closure|K_legacy_ont|K_strict_gate"
        ),
    }
    return {"tool": "rg", "mode": "content-first semantic audit for neural information-flux beta criticality", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def kernel_sq_weight(d: float, beta: float, include_phase: bool = True) -> float:
    phase = math.cos(OMEGA * d + PHI) ** 2 if include_phase else 0.5
    return phase * (d ** (DF - 1.0)) / ((1.0 + beta * (d ** ETA)) ** 2)


def simpson_integral(func: Callable[[float], float], a: float, b: float, n: int = 20000) -> float:
    if n % 2:
        n += 1
    h = (b - a) / n
    total = func(a) + func(b)
    odd = 0.0
    even = 0.0
    for i in range(1, n):
        val = func(a + i * h)
        if i % 2:
            odd += val
        else:
            even += val
    return h * (total + 4.0 * odd + 2.0 * even) / 3.0


def flux(beta: float, include_phase: bool = True, dmax: float = DOMAIN_DMAX) -> float:
    return simpson_integral(lambda d: kernel_sq_weight(d, beta, include_phase=include_phase), DOMAIN_DMIN, dmax)


def flux_derivative_integral(beta: float, include_phase: bool = True, dmax: float = DOMAIN_DMAX) -> float:
    def integrand(d: float) -> float:
        phase = math.cos(OMEGA * d + PHI) ** 2 if include_phase else 0.5
        return -2.0 * phase * (d ** (DF - 1.0 + ETA)) / ((1.0 + beta * (d ** ETA)) ** 3)
    return simpson_integral(integrand, DOMAIN_DMIN, dmax)


def entropy_for_beta(beta: float, bins: int = 4000, dmax: float = DOMAIN_DMAX) -> dict[str, float]:
    # Differential entropy and effective support of the normalized squared-kernel attention density.
    # This is a standard throughput proxy; it is not asserted to be the unique FIN mutual information.
    dx = (dmax - DOMAIN_DMIN) / bins
    weights = []
    points = []
    for i in range(bins):
        d = DOMAIN_DMIN + (i + 0.5) * dx
        w = kernel_sq_weight(d, beta, include_phase=True)
        weights.append(w)
        points.append(d)
    z = sum(weights) * dx
    probs = [w / z for w in weights]
    entropy = -sum(p * math.log(max(p, 1e-300)) for p in probs) * dx
    mean_d = sum(p * d for p, d in zip(probs, points)) * dx
    second_d = sum(p * d * d for p, d in zip(probs, points)) * dx
    return {
        "beta": beta,
        "flux": z,
        "differential_entropy_proxy": entropy,
        "effective_length_exp_entropy": math.exp(entropy),
        "mean_distance": mean_d,
        "std_distance": math.sqrt(max(second_d - mean_d * mean_d, 0.0)),
    }


def flux_conservation_certificate() -> dict[str, Any]:
    beta_rows = []
    for beta in [0.25, 0.5, 1.0, 1.1473958, 2.0, 4.0]:
        f = flux(beta, include_phase=True)
        beta_rows.append({
            "beta": beta,
            "finite_window_flux": f,
            "relative_to_beta_1_flux": f / flux(BETA_STRICT, include_phase=True),
            "flux_derivative": flux_derivative_integral(beta, include_phase=True),
        })
    phase_averaged_beta1 = flux(BETA_STRICT, include_phase=False)
    return {
        "theorem": "For beta>0 on a fixed finite window, F(beta)=int |K_beta(d)|^2 dmu_f(d) is strictly decreasing because dF/dbeta is a negative integral wherever cos^2 is not identically zero.",
        "measure": "dmu_f(d)=d^(D_f-1) dd with D_f=eta=9/5 on 0<d<=10",
        "consequence": (
            "A conservation equation F(beta)=C selects beta only after an external value C is supplied.  "
            "For every beta0>0, choosing C=F(beta0) makes beta0 the conserved-flux solution; choosing C=F(1) makes beta=1 a tautological calibration, not a target-independent derivation."
        ),
        "phase_averaged_closed_form": (
            "If cos^2 is replaced by its resonance average 1/2 and D_f=eta, then on 0<d<infinity, "
            "F_avg(beta)=1/(2*eta*beta).  Thus F_avg(beta)=C gives beta=1/(2*eta*C), so beta=1 follows only from the extra normalization C=1/(2*eta)."
        ),
        "beta_rows": beta_rows,
        "phase_averaged_beta_1_finite_window": phase_averaged_beta1,
        "strict_beta_is_solution_only_for_const": flux(BETA_STRICT, include_phase=True),
        "exports_beta_equals_1_identity": False,
    }


def edge_of_chaos_scan() -> dict[str, Any]:
    betas = [0.1, 0.15, 0.2, 0.25, 0.35, 0.5, 0.75, 1.0, 1.1473958, 1.5, 2.0, 3.0, 5.0, 8.0]
    rows = [entropy_for_beta(beta) for beta in betas]
    max_entropy = max(rows, key=lambda row: row["differential_entropy_proxy"])
    max_effective_length = max(rows, key=lambda row: row["effective_length_exp_entropy"])
    return {
        "proxy_definition": "Normalize |K_beta(d)|^2 dmu_f on 0<d<=10 and use differential entropy/effective support as a finite-window information-throughput proxy.",
        "rows": rows,
        "max_entropy_row": max_entropy,
        "max_effective_length_row": max_effective_length,
        "beta_1_row": next(row for row in rows if row["beta"] == 1.0),
        "verdict": (
            "Under this explicit target-blind proxy the broadest information distribution is reached at the smallest audited beta, not uniquely at beta=1.  "
            "Therefore the neural edge-of-chaos analogy is useful language but does not by itself export a theorem-level critical beta identity."
        ),
        "exports_edge_of_chaos_beta_identity": False,
    }


def uv_invariance_certificate() -> dict[str, Any]:
    rows = []
    invariant = BETA_STRICT * (1.0 ** ETA)
    for scale in [0.5, 1.0, 2.0, 10.0]:
        # d' = scale*d; beta' must be beta/scale^eta to keep beta*d^eta invariant.
        beta_prime = BETA_STRICT / (scale ** ETA)
        rows.append({
            "length_unit_scale_a": scale,
            "beta_prime_preserving_denominator": beta_prime,
            "invariant_beta_times_length_unit_to_eta": beta_prime * (scale ** ETA),
            "bare_beta_equals_1_after_rescaling": beta_prime == 1.0,
        })
    return {
        "group_action": "d -> a*d requires beta -> beta/a^eta to preserve 1+beta*d^eta",
        "invariant": "beta * L^eta, not the bare numeral beta",
        "rows": rows,
        "consequence": (
            "Information conservation can be formulated invariantly, but the bare statement beta=1 is still a convention unless a separate canonical length/UV normalization fixes L=1.  "
            "This does not close the P2629 normalization-gauge obstruction or the P2630 source-type separation."
        ),
        "exports_uv_independent_bare_beta_identity": False,
    }


def admissible_beta_classification() -> dict[str, Any]:
    return {
        "finite_window_norm_stability": {
            "beta_positive": "admissible as a stable damping denominator; flux finite and strictly decreasing in beta",
            "beta_zero": "finite on a finite window but not a heavy-tailed damping source and fails the intended strict denominator semantics",
            "beta_negative": "rejected because 1+beta*d^eta can vanish on positive d, creating poles/sign instabilities",
        },
        "infinite_phase_averaged_tail": {
            "beta_positive": "finite with F_avg(beta)=1/(2*eta*beta) when D_f=eta",
            "beta_zero": "diverges at infinity for the undamped phase-averaged measure",
            "beta_negative": "rejected by a positive-distance pole",
        },
        "critical_beta_status": "No unique beta=1 critical point follows without an extra target-independent conservation constant, canonical length/UV unit, or phase-stability functional whose extremum is proven at beta=1 before comparison to K_strict_gate.",
    }


def source_acceptance(flux_cert: dict[str, Any], chaos_cert: dict[str, Any], uv_cert: dict[str, Any]) -> dict[str, Any]:
    gates = {
        "flux_conservation_selects_beta_1_without_calibrating_const": False,
        "edge_of_chaos_proxy_has_unique_beta_1_maximum": bool(chaos_cert["max_entropy_row"]["beta"] == 1.0),
        "bare_beta_1_invariant_under_uv_length_rescaling": False,
        "normalization_invariant_match_beta_micro_over_beta_strict_equals_1": False,
        "no_role_transfer_or_toe_overclaim": True,
    }
    return {
        "gates": gates,
        "accepts_information_conservation_beta_identity": all(gates.values()),
        "failed_gates": [name for name, value in gates.items() if not value],
        "status": "OBSTRUCTION_CLASSIFICATION_NO_BETA_1_IDENTITY_EXPORT",
        "reason": (
            "P2631 turns the neural-network prism into finite mathematical tests.  The flux equation is monotone but calibratable to any beta0>0; "
            "the explicit entropy-throughput proxy does not peak at beta=1; and bare beta=1 is not invariant under length/UV rescaling.  "
            "Thus this packet classifies the obstruction instead of promoting a new positive_beta_renormalization_source."
        ),
    }


def recommendation() -> str:
    return (
        "Next honest step: stop treating the neural edge-of-chaos analogy as a proof of beta=1.  Either (i) derive an independent, "
        "dimensionless conservation constant and canonical length/UV unit from nadsoliton dynamics that makes F(beta)=F(1) noncircular, plus a "
        "phase-stability functional whose unique extremum is beta=1 before any strict-kernel comparison; or (ii) explicitly downgrade the damping bridge "
        "to an approximate finite-domain effective-kernel theorem with declared beta range, epsilon, and domain."
    )


def write_markdown(payload: dict[str, Any]) -> None:
    flux_cert = payload["flux_conservation_certificate"]
    chaos = payload["edge_of_chaos_scan"]
    uv = payload["uv_invariance_certificate"]
    acceptance = payload["source_acceptance"]
    lines = [
        "# P2631/S1581 neural information-flux beta criticality audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps neural/information-critical research content rather than only ticket names or numbers before adding the audit.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Flux conservation result",
        "",
        flux_cert["theorem"],
        "",
        flux_cert["consequence"],
        "",
        flux_cert["phase_averaged_closed_form"],
        "",
        "| beta | finite-window flux | relative to beta=1 | dF/dbeta |",
        "| ---: | ---: | ---: | ---: |",
    ])
    for row in flux_cert["beta_rows"]:
        lines.append(f"| {row['beta']:.7g} | {row['finite_window_flux']:.10f} | {row['relative_to_beta_1_flux']:.10f} | {row['flux_derivative']:.10f} |")
    lines.extend([
        "",
        "## Edge-of-chaos proxy scan",
        "",
        chaos["proxy_definition"],
        "",
        f"Maximum audited entropy proxy occurs at `beta={chaos['max_entropy_row']['beta']}`, not uniquely at `beta=1`.",
        "",
        chaos["verdict"],
        "",
        "## UV/length normalization invariance",
        "",
        uv["consequence"],
        "",
        "## Verdict",
        "",
        acceptance["reason"],
        "",
        "P2631 therefore does not repair P2625/P2629/P2630 and does not reopen bridge completion, role-transfer, role-bearing `L_total`, `QW-2191`, or ToE closure.",
        "",
        "## Recommended next honest step",
        "",
        payload["recommended_next_honest_step"],
        "",
    ])
    MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2629 = load_json(P2629_JSON)
    p2630 = load_json(P2630_JSON)
    flux_cert = flux_conservation_certificate()
    chaos_cert = edge_of_chaos_scan()
    uv_cert = uv_invariance_certificate()
    acceptance = source_acceptance(flux_cert, chaos_cert, uv_cert)
    payload: dict[str, Any] = {
        "status": "P2631_NEURAL_INFORMATION_FLUX_BETA_CRITICALITY_AUDIT_OBSTRUCTION_NO_SOURCE_EXPORT_NO_BRIDGE_NO_LTOTAL_NO_QW2191_NO_TOE",
        "sources": {
            "P2629_JSON": rel(P2629_JSON),
            "P2630_JSON": rel(P2630_JSON),
        },
        "upstream_state": {
            "p2629_positive_beta_source_accepts": bool(p2629.get("exact_source_gate", {}).get("accepts_positive_beta_renormalization_source", False)),
            "p2630_current_bridge_accepts": bool(p2630.get("bridge_source_truth_table", {}).get("current_accepts_bridge_positive_zbeta_source", False)),
        },
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "neural_prism_mapping": {
            "cos_omega_d_phi": "sinusoidal positional/resonance encoding",
            "one_over_1_plus_beta_d_eta": "heavy-tailed attention bias on a fractal distance measure",
            "beta": "length/temperature parameter whose bare numeric value is not invariant until a length/UV unit is fixed",
            "eta_equals_df": "tail exponent equals fractal spectral measure exponent in the audited strict layer",
            "guard": "analogy is heuristic unless converted into invariant conservation and stability theorems",
        },
        "flux_conservation_certificate": flux_cert,
        "edge_of_chaos_scan": chaos_cert,
        "uv_invariance_certificate": uv_cert,
        "admissible_beta_classification": admissible_beta_classification(),
        "source_acceptance": acceptance,
        "recommended_next_honest_step": recommendation(),
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["fingerprint_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "fingerprint_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)
    append_once(
        ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
        "## P2631/S1581 neural information-flux beta criticality audit guard",
        "\n## P2631/S1581 neural information-flux beta criticality audit guard\n\n"
        "`P2631/S1581` audits the strict kernel through the neural-network prism: the cosine factor can be read as sinusoidal positional/resonance encoding and `1/(1+beta*d^eta)` as a heavy-tailed attention bias on the `D_f=eta=9/5` fractal measure.  The finite-window information-flux functional `F(beta)=int |K_beta(d)|^2 dmu_f(d)` is strictly decreasing for `beta>0`, so `F(beta)=C` selects `beta=1` only if the conservation constant is already calibrated to `F(1)`; for every `beta0>0`, `C=F(beta0)` selects `beta0`.  The audited entropy-throughput proxy does not uniquely peak at `beta=1`, and bare `beta=1` is not invariant under `d -> a*d`, `beta -> beta/a^eta`.  Therefore the neural criticality route does not export `positive_beta_renormalization_source`, does not close the P2629/P2630 obstruction, and does not reopen role-transfer, role-bearing `L_total`, `QW-2191`, or ToE closure.\n",
    )
    append_once(
        ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
        "## P2631/S1581 neural beta criticality Ltotal guard",
        "\n## P2631/S1581 neural beta criticality Ltotal guard\n\n"
        "`P2631/S1581` does not license a role-bearing `L_total` damping term from neural edge-of-chaos language.  Information-flux conservation needs an independent dimensionless constant and canonical UV/length normalization before the bare value `beta=1` can be read as a theorem rather than calibration.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(json.dumps({
        "status": result["status"],
        "flux_beta_1_const": result["flux_conservation_certificate"]["strict_beta_is_solution_only_for_const"],
        "entropy_max_beta": result["edge_of_chaos_scan"]["max_entropy_row"]["beta"],
        "accepts_beta_identity": result["source_acceptance"]["accepts_information_conservation_beta_identity"],
        "recommended_next": result["recommended_next_honest_step"],
        "out": rel(OUT),
        "md": rel(MD),
    }, indent=2, sort_keys=True, ensure_ascii=False))
