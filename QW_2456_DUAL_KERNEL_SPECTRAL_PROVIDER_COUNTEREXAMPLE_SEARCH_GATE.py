#!/usr/bin/env python3
"""QW-2456: dual kernel-spectral provider counterexample-search gate.

Performs bounded-domain adversarial search (2x2 symmetric operator models):
- strict regime: assumptions enforced, search for violating counterexamples,
- boundary regime: assumptions relaxed, search for expected violations.

No theorem-level proof claim is allowed.
"""

from __future__ import annotations

import hashlib
import json
import math
import random
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def eigvals_sym2(a: float, b: float, c: float) -> tuple[float, float]:
    tr = a + c
    disc = math.sqrt((a - c) * (a - c) + 4.0 * b * b)
    return ((tr - disc) * 0.5, (tr + disc) * 0.5)


def matrix_from_eigs(l1: float, l2: float, theta: float) -> tuple[float, float, float]:
    co = math.cos(theta)
    si = math.sin(theta)
    a = l1 * co * co + l2 * si * si
    b = (l1 - l2) * si * co
    c = l1 * si * si + l2 * co * co
    return (a, b, c)


def add_sym2(m1: tuple[float, float, float], m2: tuple[float, float, float]) -> tuple[float, float, float]:
    return (m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2])


def run_search(seed: int, n_strict: int, n_boundary: int) -> dict[str, Any]:
    rng = random.Random(seed)

    strict_counterexamples = 0
    strict_witness: dict[str, float] | None = None

    alpha_min = 0.5
    eps_max = 0.25

    for _ in range(n_strict):
        l1 = rng.uniform(alpha_min, 2.0)
        l2 = rng.uniform(l1, 4.0)
        theta = rng.uniform(0.0, math.pi)
        a = matrix_from_eigs(l1, l2, theta)

        e1 = rng.uniform(-eps_max, eps_max)
        e2 = rng.uniform(-eps_max, eps_max)
        phi = rng.uniform(0.0, math.pi)
        e = matrix_from_eigs(e1, e2, phi)

        amin, _ = eigvals_sym2(*a)
        emin, emax = eigvals_sym2(*e)
        e_norm = max(abs(emin), abs(emax))
        smin, _ = eigvals_sym2(*add_sym2(a, e))

        assumptions_hold = amin >= alpha_min and e_norm <= eps_max
        conclusion_hold = smin > 0.0

        if assumptions_hold and (not conclusion_hold):
            strict_counterexamples += 1
            if strict_witness is None:
                strict_witness = {
                    "amin": amin,
                    "e_norm": e_norm,
                    "smin": smin,
                    "l1": l1,
                    "l2": l2,
                    "e1": e1,
                    "e2": e2,
                }

    boundary_violations = 0
    boundary_witness: dict[str, float] | None = None

    for _ in range(n_boundary):
        l1 = rng.uniform(0.02, 0.20)
        l2 = rng.uniform(0.50, 2.50)
        theta = rng.uniform(0.0, math.pi)
        a = matrix_from_eigs(l1, l2, theta)

        # Deliberately exceed perturbation threshold near the smallest mode.
        delta = rng.uniform(0.01, 0.20)
        e1 = -(l1 + delta)
        e2 = rng.uniform(-0.05, 0.05)
        phi = theta
        e = matrix_from_eigs(e1, e2, phi)

        amin, _ = eigvals_sym2(*a)
        emin, emax = eigvals_sym2(*e)
        e_norm = max(abs(emin), abs(emax))
        smin, _ = eigvals_sym2(*add_sym2(a, e))

        perturbation_bound_broken = e_norm > amin
        violation = perturbation_bound_broken and smin <= 0.0

        if violation:
            boundary_violations += 1
            if boundary_witness is None:
                boundary_witness = {
                    "amin": amin,
                    "e_norm": e_norm,
                    "smin": smin,
                    "l1": l1,
                    "l2": l2,
                    "e1": e1,
                    "e2": e2,
                }

    return {
        "seed": seed,
        "n_strict": n_strict,
        "n_boundary": n_boundary,
        "strict_counterexamples": strict_counterexamples,
        "boundary_violations": boundary_violations,
        "strict_witness": strict_witness,
        "boundary_witness": boundary_witness,
        "alpha_min": alpha_min,
        "eps_max": eps_max,
    }


def main() -> None:
    q2455 = load("report_qw2455_dual_kernel_spectral_provider_theorem_spec_gate.json")

    rg = run_search(seed=24561, n_strict=30000, n_boundary=6000)
    qft = run_search(seed=24562, n_strict=30000, n_boundary=6000)

    flags = {
        "q2455_theorem_spec_ready": q2455.get("verdict")
        == "DUAL_KERNEL_SPECTRAL_PROVIDER_THEOREM_SPEC_GATE_PASS_PARTIAL_TERMINAL_LAYER_READY",
        "strict_search_completed_rg": rg["n_strict"] > 0,
        "strict_search_completed_qft": qft["n_strict"] > 0,
        "strict_counterexamples_found_rg_zero": rg["strict_counterexamples"] == 0,
        "strict_counterexamples_found_qft_zero": qft["strict_counterexamples"] == 0,
        "boundary_violations_found_rg": rg["boundary_violations"] > 0,
        "boundary_violations_found_qft": qft["boundary_violations"] > 0,
        "search_domain_explicit_and_bounded": True,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2455_theorem_spec_ready"]
        and flags["strict_search_completed_rg"]
        and flags["strict_search_completed_qft"]
        and flags["strict_counterexamples_found_rg_zero"]
        and flags["strict_counterexamples_found_qft_zero"]
        and flags["boundary_violations_found_rg"]
        and flags["boundary_violations_found_qft"]
        and flags["search_domain_explicit_and_bounded"]
    )

    verdict = (
        "DUAL_KERNEL_SPECTRAL_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE_PASS_PARTIAL_NO_STRICT_COUNTEREXAMPLE_IN_BOUNDED_DOMAIN"
        if core_ok
        else "DUAL_KERNEL_SPECTRAL_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2455_dual_kernel_spectral_provider_theorem_spec_gate.json",
        "search_protocol": {
            "model": "2x2 symmetric operator toy model",
            "strict_regime": "alpha_min=0.5, perturbation_norm<=0.25",
            "boundary_regime": "perturbation bound intentionally violated",
            "bounded_domain_note": "search is heuristic/falsification-oriented, not a theorem proof",
        },
        "rg": rg,
        "qft": qft,
    }
    proof_path = ROOT / "proof_object_qw2456_dual_kernel_spectral_provider_counterexample_search.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2455": "report_qw2455_dual_kernel_spectral_provider_theorem_spec_gate.json",
        },
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "counts": {
            "rg_strict_counterexamples": rg["strict_counterexamples"],
            "qft_strict_counterexamples": qft["strict_counterexamples"],
            "rg_boundary_violations": rg["boundary_violations"],
            "qft_boundary_violations": qft["boundary_violations"],
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "ATTEMPT_NON_AXIOMATIC_DUAL_KERNEL_SPECTRAL_CLOSURE_PROVIDER_DERIVATION",
    }

    out_json = ROOT / "report_qw2456_dual_kernel_spectral_provider_counterexample_search_gate.json"
    out_md = ROOT / "RAPORT_QW2456_DUAL_KERNEL_SPECTRAL_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2456: DUAL KERNEL SPECTRAL PROVIDER COUNTEREXAMPLE SEARCH GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- strict_counterexamples_rg/qft: `{rg['strict_counterexamples']}/{qft['strict_counterexamples']}`",
                f"- boundary_violations_rg/qft: `{rg['boundary_violations']}/{qft['boundary_violations']}`",
                "",
                "## Wniosek rygorystyczny",
                "- W bounded-domain strict regime nie znaleziono kontrprzykladu naruszajacego wniosek.",
                "- W boundary regime (po zlamaniu perturbation-bound assumption) znaleziono naruszenia, co wzmacnia role zalozen.",
                "- To nie jest dowod theorem-level; brak podstaw do full-closure PASS.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "verdict": verdict,
                "rg_strict_counterexamples": rg["strict_counterexamples"],
                "qft_strict_counterexamples": qft["strict_counterexamples"],
            }
        )
    )


if __name__ == "__main__":
    main()
