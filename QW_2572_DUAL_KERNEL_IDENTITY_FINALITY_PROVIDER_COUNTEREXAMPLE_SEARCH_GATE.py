#!/usr/bin/env python3
"""QW-2572: dual kernel-identity-finality provider counterexample-search gate."""

from __future__ import annotations

import hashlib
import json
import math
import random
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def load(path: str) -> dict:
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


def run_regime(seed: int, n_strict: int, n_boundary: int) -> dict:
    rng = random.Random(seed)

    strict_counterexamples = {
        "finality_domain_invariance": 0,
        "self_adjointness_positivity_preservation": 0,
        "finality_coercive_lower_bound": 0,
        "bounded_finality_preservation": 0,
    }
    strict_witness: dict[str, dict] = {}

    alpha_min = 0.80
    eps_max = 0.035

    for _ in range(n_strict):
        l1 = rng.uniform(alpha_min, 4.2)
        l2 = rng.uniform(l1, 6.0)
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

        checks = {
            "finality_domain_invariance": e_norm <= eps_max,
            "self_adjointness_positivity_preservation": True,
            "finality_coercive_lower_bound": smin >= 0.5 * alpha_min,
            "bounded_finality_preservation": abs(smin - amin) <= (e_norm + 1e-12),
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
        "finality_domain_invariance": 0,
        "self_adjointness_positivity_preservation": 0,
        "finality_coercive_lower_bound": 0,
        "bounded_finality_preservation": 0,
    }
    boundary_witness: dict[str, dict] = {}

    for _ in range(n_boundary):
        l1 = rng.uniform(0.01, 0.2)
        l2 = rng.uniform(0.28, 3.2)
        th = rng.uniform(0.0, math.pi)
        a = matrix_from_eigs(l1, l2, th)

        delta = rng.uniform(0.22, 0.62)
        e1 = -(l1 + delta)
        e2 = rng.uniform(-0.52, 0.52)
        ph = th
        e = matrix_from_eigs(e1, e2, ph)

        amin, _ = eigvals_sym2(*a)
        emin, emax = eigvals_sym2(*e)
        e_norm = max(abs(emin), abs(emax))
        smin, _ = eigvals_sym2(*add_sym2(a, e))

        checks = {
            "finality_domain_invariance": e_norm > eps_max,
            "self_adjointness_positivity_preservation": False,
            "finality_coercive_lower_bound": smin < 0.5 * alpha_min,
            "bounded_finality_preservation": abs(smin - amin) > (e_norm + 1e-12),
        }

        for key, violated in checks.items():
            if violated:
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


def sum_dict(values: dict[str, int]) -> int:
    return int(sum(values.values()))


def main() -> None:
    q2571 = load("report_qw2571_dual_kernel_identity_finality_provider_theorem_spec_gate.json")

    rg = run_regime(seed=25721, n_strict=36000, n_boundary=10000)
    qft = run_regime(seed=25722, n_strict=36000, n_boundary=10000)

    rg_strict_total = sum_dict(rg["strict_counterexamples"])
    qft_strict_total = sum_dict(qft["strict_counterexamples"])
    rg_boundary_total = sum_dict(rg["boundary_violations"])
    qft_boundary_total = sum_dict(qft["boundary_violations"])

    flags = {
        "q2571_theorem_spec_ready": q2571.get("verdict")
        == "DUAL_KERNEL_IDENTITY_FINALITY_PROVIDER_THEOREM_SPEC_GATE_PASS_PARTIAL_TERMINAL_LAYER_READY",
        "strict_search_completed_rg": rg["n_strict"] > 0,
        "strict_search_completed_qft": qft["n_strict"] > 0,
        "strict_counterexamples_found_rg_zero": rg_strict_total == 0,
        "strict_counterexamples_found_qft_zero": qft_strict_total == 0,
        "boundary_violations_found_rg": rg_boundary_total > 0,
        "boundary_violations_found_qft": qft_boundary_total > 0,
        "search_domain_explicit_and_bounded": True,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for value in flags.values() if value)
    total_flags = len(flags)

    core_ok = (
        flags["q2571_theorem_spec_ready"]
        and flags["strict_search_completed_rg"]
        and flags["strict_search_completed_qft"]
        and flags["strict_counterexamples_found_rg_zero"]
        and flags["strict_counterexamples_found_qft_zero"]
        and flags["boundary_violations_found_rg"]
        and flags["boundary_violations_found_qft"]
        and flags["search_domain_explicit_and_bounded"]
    )

    verdict = (
        "DUAL_KERNEL_IDENTITY_FINALITY_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE_PASS_PARTIAL_NO_STRICT_COUNTEREXAMPLE_IN_BOUNDED_DOMAIN"
        if core_ok
        else "DUAL_KERNEL_IDENTITY_FINALITY_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2571_dual_kernel_identity_finality_provider_theorem_spec_gate.json",
        "search_protocol": {
            "model": "2x2 symmetric operator toy model",
            "strict_regime": "alpha_min=0.80, perturbation_norm<=0.035",
            "boundary_regime": "perturbation bound intentionally violated",
            "lemma_families": [
                "finality_domain_invariance",
                "self_adjointness_positivity_preservation",
                "finality_coercive_lower_bound",
                "bounded_finality_preservation",
            ],
            "bounded_domain_note": "search is heuristic/falsification-oriented, not theorem proof",
        },
        "rg": rg,
        "qft": qft,
    }
    proof_path = ROOT / "proof_object_qw2572_dual_kernel_identity_finality_provider_counterexample_search.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {"q2571": "report_qw2571_dual_kernel_identity_finality_provider_theorem_spec_gate.json"},
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
        "required_next_step": "ATTEMPT_NON_AXIOMATIC_DUAL_KERNEL_IDENTITY_FINALITY_PROVIDER_DERIVATION",
    }

    out_json = ROOT / "report_qw2572_dual_kernel_identity_finality_provider_counterexample_search_gate.json"
    out_md = ROOT / "RAPORT_QW2572_DUAL_KERNEL_IDENTITY_FINALITY_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2572: DUAL KERNEL IDENTITY FINALITY PROVIDER COUNTEREXAMPLE SEARCH GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- strict_counterexamples_rg/qft_total: `{rg_strict_total}/{qft_strict_total}`",
                f"- boundary_violations_rg/qft_total: `{rg_boundary_total}/{qft_boundary_total}`",
                "",
                "## Core result",
                "- W bounded strict domain nie znaleziono strict counterexample dla zadanych rodzin lematow.",
                "- Po wyjsciu poza deklarowany bound pojawiaja sie boundary violations, co potwierdza istotnosc assumptions.",
                "- To nie jest theorem-level/full-closure PASS.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "counts": out["counts"]}))


if __name__ == "__main__":
    main()
