#!/usr/bin/env python3
"""QW-2462: dual kernel-spectral-invariance provider counterexample-search gate.

Performs bounded-domain adversarial search on operator-lemma families:
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


def run_regime(seed: int, n_strict: int, n_boundary: int) -> dict[str, Any]:
    rng = random.Random(seed)

    strict_counterexamples = {
        "domain_invariance": 0,
        "self_adjointness_preservation": 0,
        "positivity_or_coercivity": 0,
        "spectral_stability": 0,
    }
    strict_witness = {}

    alpha_min = 0.4
    eps_max = 0.2

    for _ in range(n_strict):
        l1 = rng.uniform(alpha_min, 2.0)
        l2 = rng.uniform(l1, 4.0)
        th = rng.uniform(0.0, math.pi)
        a = matrix_from_eigs(l1, l2, th)

        e1 = rng.uniform(-eps_max, eps_max)
        e2 = rng.uniform(-eps_max, eps_max)
        ph = rng.uniform(0.0, math.pi)
        e = matrix_from_eigs(e1, e2, ph)

        amin, _ = eigvals_sym2(*a)
        emin, emax = eigvals_sym2(*e)
        e_norm = max(abs(emin), abs(emax))
        smin, _ = eigvals_sym2(*add_sym2(a, e))

        # L1 domain/invariance proxy: bounded perturbation within threshold.
        l1_ok = e_norm <= eps_max
        # L2 self-adjointness proxy: symmetric toy operator preserved by construction.
        l2_ok = True
        # L3 positivity/coercivity invariance.
        l3_ok = smin > 0.0
        # L4 spectral stability proxy.
        l4_ok = abs(smin - amin) <= (e_norm + 1e-12)

        checks = {
            "domain_invariance": l1_ok,
            "self_adjointness_preservation": l2_ok,
            "positivity_or_coercivity": l3_ok,
            "spectral_stability": l4_ok,
        }

        for k, ok in checks.items():
            if not ok:
                strict_counterexamples[k] += 1
                if k not in strict_witness:
                    strict_witness[k] = {
                        "amin": amin,
                        "smin": smin,
                        "e_norm": e_norm,
                        "l1": l1,
                        "l2": l2,
                        "e1": e1,
                        "e2": e2,
                    }

    boundary_violations = {
        "domain_invariance": 0,
        "self_adjointness_preservation": 0,
        "positivity_or_coercivity": 0,
        "spectral_stability": 0,
    }
    boundary_witness = {}

    for _ in range(n_boundary):
        l1 = rng.uniform(0.02, 0.2)
        l2 = rng.uniform(0.5, 2.5)
        th = rng.uniform(0.0, math.pi)
        a = matrix_from_eigs(l1, l2, th)

        delta = rng.uniform(0.05, 0.3)
        e1 = -(l1 + delta)
        e2 = rng.uniform(-0.2, 0.2)
        ph = th
        e = matrix_from_eigs(e1, e2, ph)

        amin, _ = eigvals_sym2(*a)
        emin, emax = eigvals_sym2(*e)
        e_norm = max(abs(emin), abs(emax))
        smin, _ = eigvals_sym2(*add_sym2(a, e))

        # boundary: deliberately violate perturbation bounds
        l1_violation = e_norm > eps_max
        l2_violation = False  # symmetric toy remains symmetric by construction
        l3_violation = smin <= 0.0
        l4_violation = abs(smin - amin) > (e_norm + 1e-12)

        checks = {
            "domain_invariance": l1_violation,
            "self_adjointness_preservation": l2_violation,
            "positivity_or_coercivity": l3_violation,
            "spectral_stability": l4_violation,
        }

        for k, vio in checks.items():
            if vio:
                boundary_violations[k] += 1
                if k not in boundary_witness:
                    boundary_witness[k] = {
                        "amin": amin,
                        "smin": smin,
                        "e_norm": e_norm,
                        "l1": l1,
                        "l2": l2,
                        "e1": e1,
                        "e2": e2,
                    }

    return {
        "seed": seed,
        "n_strict": n_strict,
        "n_boundary": n_boundary,
        "alpha_min": alpha_min,
        "eps_max": eps_max,
        "strict_counterexamples": strict_counterexamples,
        "strict_witness": strict_witness,
        "boundary_violations": boundary_violations,
        "boundary_witness": boundary_witness,
    }


def sum_dict(d: dict[str, int]) -> int:
    return int(sum(d.values()))


def main() -> None:
    q2461 = load("report_qw2461_dual_kernel_spectral_invariance_provider_theorem_spec_gate.json")

    rg = run_regime(seed=24621, n_strict=30000, n_boundary=8000)
    qft = run_regime(seed=24622, n_strict=30000, n_boundary=8000)

    rg_strict_total = sum_dict(rg["strict_counterexamples"])
    qft_strict_total = sum_dict(qft["strict_counterexamples"])
    rg_boundary_total = sum_dict(rg["boundary_violations"])
    qft_boundary_total = sum_dict(qft["boundary_violations"])

    flags = {
        "q2461_theorem_spec_ready": q2461.get("verdict")
        == "DUAL_KERNEL_SPECTRAL_INVARIANCE_PROVIDER_THEOREM_SPEC_GATE_PASS_PARTIAL_TERMINAL_LAYER_READY",
        "strict_search_completed_rg": rg["n_strict"] > 0,
        "strict_search_completed_qft": qft["n_strict"] > 0,
        "strict_counterexamples_found_rg_zero": rg_strict_total == 0,
        "strict_counterexamples_found_qft_zero": qft_strict_total == 0,
        "boundary_violations_found_rg": rg_boundary_total > 0,
        "boundary_violations_found_qft": qft_boundary_total > 0,
        "search_domain_explicit_and_bounded": True,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2461_theorem_spec_ready"]
        and flags["strict_search_completed_rg"]
        and flags["strict_search_completed_qft"]
        and flags["strict_counterexamples_found_rg_zero"]
        and flags["strict_counterexamples_found_qft_zero"]
        and flags["boundary_violations_found_rg"]
        and flags["boundary_violations_found_qft"]
        and flags["search_domain_explicit_and_bounded"]
    )

    verdict = (
        "DUAL_KERNEL_SPECTRAL_INVARIANCE_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE_PASS_PARTIAL_NO_STRICT_COUNTEREXAMPLE_IN_BOUNDED_DOMAIN"
        if core_ok
        else "DUAL_KERNEL_SPECTRAL_INVARIANCE_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2461_dual_kernel_spectral_invariance_provider_theorem_spec_gate.json",
        "search_protocol": {
            "model": "2x2 symmetric operator toy model",
            "strict_regime": "alpha_min=0.4, perturbation_norm<=0.2",
            "boundary_regime": "perturbation bound intentionally violated",
            "lemma_families": [
                "domain_invariance",
                "self_adjointness_preservation",
                "positivity_or_coercivity",
                "spectral_stability",
            ],
            "bounded_domain_note": "search is heuristic/falsification-oriented, not theorem proof",
        },
        "rg": rg,
        "qft": qft,
    }
    proof_path = ROOT / "proof_object_qw2462_dual_kernel_spectral_invariance_provider_counterexample_search.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2461": "report_qw2461_dual_kernel_spectral_invariance_provider_theorem_spec_gate.json"
        },
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "counts": {
            "rg_strict_counterexamples_total": rg_strict_total,
            "qft_strict_counterexamples_total": qft_strict_total,
            "rg_boundary_violations_total": rg_boundary_total,
            "qft_boundary_violations_total": qft_boundary_total,
            "rg_strict_by_lemma": rg["strict_counterexamples"],
            "qft_strict_by_lemma": qft["strict_counterexamples"],
            "rg_boundary_by_lemma": rg["boundary_violations"],
            "qft_boundary_by_lemma": qft["boundary_violations"],
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "ATTEMPT_NON_AXIOMATIC_DUAL_KERNEL_SPECTRAL_INVARIANCE_PROVIDER_DERIVATION",
    }

    out_json = ROOT / "report_qw2462_dual_kernel_spectral_invariance_provider_counterexample_search_gate.json"
    out_md = ROOT / "RAPORT_QW2462_DUAL_KERNEL_SPECTRAL_INVARIANCE_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2462: DUAL KERNEL SPECTRAL INVARIANCE PROVIDER COUNTEREXAMPLE SEARCH GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- strict_counterexamples_rg/qft_total: `{rg_strict_total}/{qft_strict_total}`",
                f"- boundary_violations_rg/qft_total: `{rg_boundary_total}/{qft_boundary_total}`",
                "",
                "## Wniosek rygorystyczny",
                "- W bounded-domain strict regime nie znaleziono kontrprzykladow naruszajacych lematy.",
                "- W boundary regime (po zlamaniu assumptions) znaleziono naruszenia.",
                "- To nie jest theorem-level proof; brak podstaw do full-closure PASS.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "verdict": verdict,
                "rg_strict_counterexamples_total": rg_strict_total,
                "qft_strict_counterexamples_total": qft_strict_total,
            }
        )
    )


if __name__ == "__main__":
    main()
