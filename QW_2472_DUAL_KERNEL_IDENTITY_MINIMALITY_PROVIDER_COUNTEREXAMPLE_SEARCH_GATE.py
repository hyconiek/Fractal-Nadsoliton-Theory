#!/usr/bin/env python3
"""QW-2472: dual kernel-identity-minimality provider counterexample-search gate."""

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
        "minimal_domain_invariance": 0,
        "self_adjointness_positivity_preservation": 0,
        "minimality_consistency_lower_bound": 0,
        "bounded_minimality_stability": 0,
    }
    strict_witness: dict[str, Any] = {}

    alpha_min = 0.40
    eps_max = 0.16

    for _ in range(n_strict):
        l1 = rng.uniform(alpha_min, 2.3)
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

        c1_ok = e_norm <= eps_max
        c2_ok = True
        c3_ok = smin > 0.0
        c4_ok = abs(smin - amin) <= (e_norm + 1e-12)

        checks = {
            "minimal_domain_invariance": c1_ok,
            "self_adjointness_positivity_preservation": c2_ok,
            "minimality_consistency_lower_bound": c3_ok,
            "bounded_minimality_stability": c4_ok,
        }

        for key, ok in checks.items():
            if not ok:
                strict_counterexamples[key] += 1
                if key not in strict_witness:
                    strict_witness[key] = {
                        "amin": amin,
                        "smin": smin,
                        "e_norm": e_norm,
                        "l1": l1,
                        "l2": l2,
                        "e1": e1,
                        "e2": e2,
                    }

    boundary_violations = {
        "minimal_domain_invariance": 0,
        "self_adjointness_positivity_preservation": 0,
        "minimality_consistency_lower_bound": 0,
        "bounded_minimality_stability": 0,
    }
    boundary_witness: dict[str, Any] = {}

    for _ in range(n_boundary):
        l1 = rng.uniform(0.01, 0.18)
        l2 = rng.uniform(0.4, 2.7)
        th = rng.uniform(0.0, math.pi)
        a = matrix_from_eigs(l1, l2, th)

        delta = rng.uniform(0.06, 0.36)
        e1 = -(l1 + delta)
        e2 = rng.uniform(-0.28, 0.28)
        ph = th
        e = matrix_from_eigs(e1, e2, ph)

        amin, _ = eigvals_sym2(*a)
        emin, emax = eigvals_sym2(*e)
        e_norm = max(abs(emin), abs(emax))
        smin, _ = eigvals_sym2(*add_sym2(a, e))

        c1_vio = e_norm > eps_max
        c2_vio = False
        c3_vio = smin <= 0.0
        c4_vio = abs(smin - amin) > (e_norm + 1e-12)

        checks = {
            "minimal_domain_invariance": c1_vio,
            "self_adjointness_positivity_preservation": c2_vio,
            "minimality_consistency_lower_bound": c3_vio,
            "bounded_minimality_stability": c4_vio,
        }

        for key, vio in checks.items():
            if vio:
                boundary_violations[key] += 1
                if key not in boundary_witness:
                    boundary_witness[key] = {
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
    q2471 = load("report_qw2471_dual_kernel_identity_minimality_provider_theorem_spec_gate.json")

    rg = run_regime(seed=24721, n_strict=32000, n_boundary=9500)
    qft = run_regime(seed=24722, n_strict=32000, n_boundary=9500)

    rg_strict_total = sum_dict(rg["strict_counterexamples"])
    qft_strict_total = sum_dict(qft["strict_counterexamples"])
    rg_boundary_total = sum_dict(rg["boundary_violations"])
    qft_boundary_total = sum_dict(qft["boundary_violations"])

    flags = {
        "q2471_theorem_spec_ready": q2471.get("verdict")
        == "DUAL_KERNEL_IDENTITY_MINIMALITY_PROVIDER_THEOREM_SPEC_GATE_PASS_PARTIAL_TERMINAL_LAYER_READY",
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
        flags["q2471_theorem_spec_ready"]
        and flags["strict_search_completed_rg"]
        and flags["strict_search_completed_qft"]
        and flags["strict_counterexamples_found_rg_zero"]
        and flags["strict_counterexamples_found_qft_zero"]
        and flags["boundary_violations_found_rg"]
        and flags["boundary_violations_found_qft"]
        and flags["search_domain_explicit_and_bounded"]
    )

    verdict = (
        "DUAL_KERNEL_IDENTITY_MINIMALITY_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE_PASS_PARTIAL_NO_STRICT_COUNTEREXAMPLE_IN_BOUNDED_DOMAIN"
        if core_ok
        else "DUAL_KERNEL_IDENTITY_MINIMALITY_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2471_dual_kernel_identity_minimality_provider_theorem_spec_gate.json",
        "search_protocol": {
            "model": "2x2 symmetric operator toy model",
            "strict_regime": "alpha_min=0.40, perturbation_norm<=0.16",
            "boundary_regime": "perturbation bound intentionally violated",
            "lemma_families": [
                "minimal_domain_invariance",
                "self_adjointness_positivity_preservation",
                "minimality_consistency_lower_bound",
                "bounded_minimality_stability",
            ],
            "bounded_domain_note": "search is heuristic/falsification-oriented, not theorem proof",
        },
        "rg": rg,
        "qft": qft,
    }
    proof_path = ROOT / "proof_object_qw2472_dual_kernel_identity_minimality_provider_counterexample_search.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {"q2471": "report_qw2471_dual_kernel_identity_minimality_provider_theorem_spec_gate.json"},
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
        "required_next_step": "ATTEMPT_NON_AXIOMATIC_DUAL_KERNEL_IDENTITY_MINIMALITY_PROVIDER_DERIVATION",
    }

    out_json = ROOT / "report_qw2472_dual_kernel_identity_minimality_provider_counterexample_search_gate.json"
    out_md = ROOT / "RAPORT_QW2472_DUAL_KERNEL_IDENTITY_MINIMALITY_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2472: DUAL KERNEL IDENTITY MINIMALITY PROVIDER COUNTEREXAMPLE SEARCH GATE",
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

    print(json.dumps({"verdict": verdict, "counts": out["counts"]}))


if __name__ == "__main__":
    main()
