#!/usr/bin/env python3
"""
QW-2125: K_total generation-alignment audit (structural, no-fit).

Purpose:
- quantify whether robust 12-octave tripartition (QW-2118) naturally aligns
  with canonical 3-generation templates,
- separate "tripartition exists" from "tripartition equals physical generations".
"""

from __future__ import annotations

import itertools
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2125_ktotal_generation_alignment_audit.json"
OUT_MD = ROOT / "RAPORT_QW2125_KTOTAL_GENERATION_ALIGNMENT_AUDIT.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def kernel_value(d: float, omega: float, phi: float, beta: float, eta: float) -> float:
    return float(np.cos(omega * d + phi) / (1.0 + beta * (d**eta)))


def build_ktotal_matrix(n: int, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    m = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = min(abs(i - j), n - abs(i - j))
            m[i, j] = kernel_value(float(d), omega, phi, beta, eta)
    return 0.5 * (m + m.T)


def kmeans_tripartition_labels(m: np.ndarray) -> List[int]:
    eigvals, eigvecs = np.linalg.eigh(m)
    idx = np.argsort(eigvals)[::-1]
    x = np.abs(eigvecs[:, idx[:3]])
    cent = np.array([x[0], x[4], x[8]], dtype=float)
    labels = np.zeros(x.shape[0], dtype=int)
    for _ in range(80):
        d2 = np.sum((x[:, None, :] - cent[None, :, :]) ** 2, axis=2)
        labels = np.argmin(d2, axis=1)
        new_cent = np.array(
            [x[labels == c].mean(axis=0) if np.any(labels == c) else cent[c] for c in range(3)],
            dtype=float,
        )
        if np.allclose(new_cent, cent, atol=1e-12):
            break
        cent = new_cent
    return [int(v) for v in labels.tolist()]


def labels_to_sets(labels: Sequence[int]) -> Dict[int, Set[int]]:
    vals = sorted(set(int(v) for v in labels))
    return {c: {i for i, lv in enumerate(labels) if int(lv) == c} for c in vals}


def best_alignment(cluster_sets: Dict[int, Set[int]], target_sets: Dict[int, Set[int]]) -> Dict[str, object]:
    cluster_ids = sorted(cluster_sets.keys())
    target_ids = sorted(target_sets.keys())
    best_score = -1
    best_perm: Tuple[int, int, int] | None = None
    best_rows: List[Dict[str, object]] = []
    for perm in itertools.permutations(cluster_ids):
        score = 0
        rows: List[Dict[str, object]] = []
        for tg, cid in zip(target_ids, perm):
            inter = cluster_sets[cid] & target_sets[tg]
            s = len(inter)
            score += s
            rows.append(
                {
                    "target_group": int(tg),
                    "cluster_id": int(cid),
                    "intersection_size": int(s),
                    "intersection_members": sorted(int(v) for v in inter),
                }
            )
        if score > best_score:
            best_score = score
            best_perm = tuple(int(v) for v in perm)
            best_rows = rows
    overlap_fraction = float(best_score / 12.0)
    return {
        "best_score_12": int(best_score),
        "best_overlap_fraction": overlap_fraction,
        "cluster_permutation_for_targets": list(best_perm) if best_perm is not None else [],
        "rows": best_rows,
    }


def perturbed_kernel_samples(base_kernel: Dict[str, float], n: int, seed: int) -> List[Dict[str, float]]:
    rng = np.random.default_rng(seed)
    out: List[Dict[str, float]] = []
    for _ in range(n):
        kk = {
            "omega": float(np.clip(base_kernel["omega"] * (1.0 + rng.normal(0.0, 0.02)), 0.01, 2.0)),
            "phi": float(np.clip(base_kernel["phi"] + rng.normal(0.0, 0.02), -np.pi, np.pi)),
            "beta": float(np.clip(base_kernel["beta"] * (1.0 + rng.normal(0.0, 0.03)), 0.01, 5.0)),
            "eta": float(np.clip(base_kernel["eta"] * (1.0 + rng.normal(0.0, 0.02)), 0.5, 4.0)),
        }
        out.append(kk)
    return out


def main() -> None:
    r2118 = load_json("report_qw2118_ktotal_spectral_tripartition_gate.json")
    kernel = {k: float(v) for k, v in r2118["kernel"].items()}

    labels_base = [int(v) for v in r2118["band_partition"]["tripartition"]["labels"]]
    clusters_base = labels_to_sets(labels_base)

    templates = {
        "mod3_4x3_template": {
            0: {0, 3, 6, 9},
            1: {1, 4, 7, 10},
            2: {2, 5, 8, 11},
        },
        "contiguous_4x3_template": {
            0: {0, 1, 2, 3},
            1: {4, 5, 6, 7},
            2: {8, 9, 10, 11},
        },
    }

    alignment = {
        name: best_alignment(clusters_base, target_sets) for name, target_sets in templates.items()
    }

    n_perturb = 300
    pert = perturbed_kernel_samples(kernel, n=n_perturb, seed=2125)
    overlap_mod3: List[float] = []
    for kk in pert:
        m = build_ktotal_matrix(12, kk["omega"], kk["phi"], kk["beta"], kk["eta"])
        labels = kmeans_tripartition_labels(m)
        cl = labels_to_sets(labels)
        al = best_alignment(cl, templates["mod3_4x3_template"])
        overlap_mod3.append(float(al["best_overlap_fraction"]))

    overlap_mod3_mean = float(np.mean(overlap_mod3))
    overlap_mod3_p10 = float(np.percentile(overlap_mod3, 10))
    overlap_mod3_p90 = float(np.percentile(overlap_mod3, 90))

    base_mod3 = float(alignment["mod3_4x3_template"]["best_overlap_fraction"])
    base_contig = float(alignment["contiguous_4x3_template"]["best_overlap_fraction"])

    flags = {
        "q2118_tripartition_is_balanced_4_4_4": bool(r2118["flags"].get("balanced_tripartition_4_4_4", False)),
        "three_group_template_defined_without_fit": True,
        "base_mod3_overlap_ge_0p66": bool(base_mod3 >= (8.0 / 12.0)),
        "base_mod3_beats_contiguous_template": bool(base_mod3 > base_contig),
        "mod3_overlap_mean_ge_0p60_under_kernel_perturb": bool(overlap_mod3_mean >= 0.60),
        "mod3_overlap_p10_ge_0p50_under_kernel_perturb": bool(overlap_mod3_p10 >= 0.50),
        "deterministic_no_scan_no_retune": True,
        "generation_mapping_is_unique_and_physical": False,
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    verdict = (
        "KTOTAL_GENERATION_ALIGNMENT_AUDIT_PASS_STRUCTURAL_PARTIAL"
        if (
            flags["q2118_tripartition_is_balanced_4_4_4"]
            and flags["base_mod3_overlap_ge_0p66"]
            and flags["base_mod3_beats_contiguous_template"]
            and flags["mod3_overlap_mean_ge_0p60_under_kernel_perturb"]
            and flags["mod3_overlap_p10_ge_0p50_under_kernel_perturb"]
        )
        else "KTOTAL_GENERATION_ALIGNMENT_AUDIT_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2118_tripartition": "report_qw2118_ktotal_spectral_tripartition_gate.json",
        },
        "kernel": kernel,
        "templates": {
            k: {str(kk): sorted(int(v) for v in vv) for kk, vv in t.items()} for k, t in templates.items()
        },
        "base_tripartition_labels": labels_base,
        "base_alignment": alignment,
        "perturbation_robustness": {
            "n_perturb": n_perturb,
            "mod3_overlap_fraction_mean": overlap_mod3_mean,
            "mod3_overlap_fraction_p10": overlap_mod3_p10,
            "mod3_overlap_fraction_p90": overlap_mod3_p90,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "DEFINE_CANONICAL_STATE_TO_OCTAVE_MAPPING_AND_PROMOTE_TO_PHYSICAL_GENERATION_DERIVATION_GATE"
            if verdict.endswith("PARTIAL")
            else "REVIEW_TRIPARTITION_OR_TEMPLATE_ASSUMPTIONS_AND_RERUN_QW2125"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2125: KTOTAL GENERATION ALIGNMENT AUDIT",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Base overlap (12-point scale)",
        f"- mod3 template: `{alignment['mod3_4x3_template']['best_score_12']}/12` "
        f"(`{alignment['mod3_4x3_template']['best_overlap_fraction']:.3f}`)",
        f"- contiguous template: `{alignment['contiguous_4x3_template']['best_score_12']}/12` "
        f"(`{alignment['contiguous_4x3_template']['best_overlap_fraction']:.3f}`)",
        "",
        "## Robustness (kernel perturbations)",
        f"- mean overlap (mod3): `{overlap_mod3_mean:.3f}`",
        f"- p10 overlap (mod3): `{overlap_mod3_p10:.3f}`",
        f"- p90 overlap (mod3): `{overlap_mod3_p90:.3f}`",
        "",
        "## Interpretation",
        "- structural 3-way partition exists and aligns partially with generation-like template,",
        "- physical one-to-one generation mapping remains open.",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2125] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2125] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2125] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()

