#!/usr/bin/env python3
"""Build an exhaustive, explicitly heuristic FIN "puzzle atlas".

The tag scores are discovery aids, not evidence of physical equivalence.
The numerical block-reduction section is independent linear algebra and is
reported separately with machine-checkable residuals.
"""

from __future__ import annotations

import csv
import itertools
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm


HERE = Path(__file__).resolve().parent
FIG = HERE / "FIN_Nadsoliton_Puzzle_Atlas_Figures"
FIG.mkdir(exist_ok=True)

N = 12
K = np.array(
    [
        0.46998567264502017,
        0.19204355169010282,
        0.09142861427792497,
        0.04702916874565040,
        0.02413122336363006,
        0.011070817321442113,
    ]
)


CHARACTERISTICS = [
    ("C01", "strict oscillatory-damped kernel", {"spectral", "graph", "nonlocal", "oscillation", "damping", "finite"}),
    ("C02", "legacy intermediate bridge profile", {"bridge", "information", "amplitude", "phase", "damping", "green"}),
    ("C03", "positive graph Laplacian A", {"graph", "spectral", "positivity", "variational", "conservation", "semigroup"}),
    ("C04", "unitary/heat dual functional calculus", {"spectral", "unitary", "diffusion", "semigroup", "analytic", "dynamics"}),
    ("C05", "Green/resolvent/action reconstruction", {"green", "resolvent", "variational", "inverse", "field", "spectral"}),
    ("C06", "orientation cocycle and chiral receivers", {"topology", "cohomology", "chirality", "symmetry", "resource", "phase"}),
    ("C07", "fractal Schur compression", {"coarse_graining", "schur", "renormalization", "multiscale", "fractal", "graph"}),
    ("C08", "adaptive operator and memory law", {"learning", "memory", "feedback", "nonlinear", "adaptive", "operator"}),
    ("C09", "information/entropy geometry", {"information", "entropy", "geometry", "probability", "metric", "inference"}),
    ("C10", "operational process tuple", {"process", "measurement", "environment", "instrument", "observer", "record"}),
    ("C11", "dimensional-scale obstruction", {"dimensions", "units", "scale", "calibration", "no_go", "provenance"}),
    ("C12", "selector/source obstruction", {"selector", "symmetry", "chirality", "resource", "no_go", "provenance"}),
]


DOMAINS = [
    ("D01", "spectral graphs and quantum walks", {"spectral", "graph", "unitary", "diffusion", "dynamics", "finite"}, "Exact at the finite operator level; not a continuum or laboratory identification."),
    ("D02", "lattice Green functions and inverse field reconstruction", {"green", "resolvent", "field", "inverse", "variational", "spectral"}, "Exact for chosen quadratic operators; interaction and physical units are extra."),
    ("D03", "graph signal processing and wavelets", {"graph", "spectral", "filter", "multiscale", "oscillation", "nonlocal"}, "Strong filter analogy; it does not supply ontology or dynamics by itself."),
    ("D04", "Kron reduction, Schur complements and multigrid", {"graph", "schur", "coarse_graining", "multiscale", "renormalization", "resolvent"}, "Exact block algebra; closure inside the strict kernel family is separately false."),
    ("D05", "Mori-Zwanzig reduction and process tensors", {"memory", "environment", "process", "coarse_graining", "dynamics", "resolvent"}, "Exact reduced-memory mechanism after a declared observed/hidden split; the split is operational input."),
    ("D06", "resource theory of asymmetry", {"resource", "symmetry", "chirality", "selector", "measurement", "no_go"}, "Exact no-creation logic for symmetric free operations; no nonpremise source is obtained."),
    ("D07", "stochastic thermodynamics and feedback", {"entropy", "information", "feedback", "measurement", "environment", "units"}, "Requires bath, temperature, energy scale and reset ledger."),
    ("D08", "information geometry and diffusion maps", {"information", "geometry", "metric", "probability", "diffusion", "inference"}, "Produces intrinsic dimensionless geometry, not SI length or spacetime."),
    ("D09", "noncommutative geometry and spectral action", {"spectral", "geometry", "operator", "variational", "dimensions", "algebra"}, "Spectral analogy is real; algebra, grading, real structure and scale remain additional data."),
    ("D10", "renormalization and stable fractional limits", {"renormalization", "scale", "multiscale", "coarse_graining", "nonlocal", "damping"}, "Possible universality language; the audited strict family is not Schur-closed."),
    ("D11", "control, observability and process tomography", {"measurement", "instrument", "record", "calibration", "inverse", "process"}, "Operationally exact when preparations, clock and POVM are supplied."),
    ("D12", "adaptive neural operators and reservoir memory", {"learning", "adaptive", "operator", "memory", "feedback", "nonlinear"}, "Mechanistic analogy; training objective and data provenance are not derived."),
    ("D13", "topological phases, cohomology and holonomy", {"topology", "cohomology", "phase", "chirality", "symmetry", "geometry"}, "Correct language for orientation carriers; no automatic signed section."),
    ("D14", "causal-operational quantum theory", {"process", "observer", "instrument", "measurement", "environment", "dynamics"}, "Clarifies interventions and records; does not choose a physical realization."),
    ("D15", "statistical field theory and random operators", {"field", "probability", "spectral", "entropy", "nonlinear", "scale"}, "Useful null ensembles and effective-action comparisons; not an equivalence theorem."),
]


def structural_score(tags: set[str], domain_tags: set[str]) -> tuple[int, list[str]]:
    shared = sorted(tags & domain_tags)
    cosine = len(shared) / math.sqrt(len(tags) * len(domain_tags))
    score = min(4, int(round(4.0 * cosine)))
    return score, shared


def write_atlas_tables() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    singles: list[dict[str, object]] = []
    for cid, cname, ctags in CHARACTERISTICS:
        for did, dname, dtags, caveat in DOMAINS:
            score, shared = structural_score(ctags, dtags)
            singles.append(
                {
                    "characteristic_id": cid,
                    "characteristic": cname,
                    "domain_id": did,
                    "domain": dname,
                    "structural_score_0_4": score,
                    "shared_tags": ";".join(shared),
                    "mandatory_caveat": caveat,
                }
            )

    pairs: list[dict[str, object]] = []
    for left, right in itertools.combinations(CHARACTERISTICS, 2):
        lid, lname, ltags = left
        rid, rname, rtags = right
        union = ltags | rtags
        for did, dname, dtags, caveat in DOMAINS:
            base, shared = structural_score(union, dtags)
            touches_both = bool(ltags & dtags) and bool(rtags & dtags)
            score = min(4, base + int(touches_both and len(shared) >= 3))
            pairs.append(
                {
                    "pair_id": f"{lid}+{rid}",
                    "pair": f"{lname} + {rname}",
                    "domain_id": did,
                    "domain": dname,
                    "structural_score_0_4": score,
                    "touches_both_puzzles": touches_both,
                    "shared_tags": ";".join(shared),
                    "mandatory_caveat": caveat,
                }
            )

    for name, rows in [
        ("FIN_Nadsoliton_Puzzle_Atlas_All_Singles.csv", singles),
        ("FIN_Nadsoliton_Puzzle_Atlas_All_Pairs.csv", pairs),
    ]:
        with (HERE / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    return singles, pairs


def strict_matrix() -> np.ndarray:
    w = np.zeros((N, N))
    for x in range(N):
        for y in range(N):
            delta = abs(x - y) % N
            d = min(delta, N - delta)
            if d:
                w[x, y] = K[d - 1]
    return w[0].sum() * np.eye(N) - w


def dynamic_schur_study() -> dict[str, object]:
    """Test the new hidden-node memory/self-energy puzzle."""
    a = strict_matrix()
    retained = np.arange(0, N, 2)
    hidden = np.arange(1, N, 2)
    a_ee = a[np.ix_(retained, retained)]
    a_eo = a[np.ix_(retained, hidden)]
    a_oe = a[np.ix_(hidden, retained)]
    a_oo = a[np.ix_(hidden, hidden)]
    size = len(retained)

    static = a_ee - a_eo @ np.linalg.solve(a_oo, a_oe)
    z_values = [0.05, 0.2, 1.0]
    self_energy = {}
    resolvent_residual = {}
    derivative_norm = {}
    for z in z_values:
        sigma = a_eo @ np.linalg.solve(z * np.eye(size) + a_oo, a_oe)
        reduced = np.linalg.inv(z * np.eye(size) + a_ee - sigma)
        exact = np.linalg.inv(z * np.eye(N) + a)[np.ix_(retained, retained)]
        sigma_prime = -a_eo @ np.linalg.solve(
            (z * np.eye(size) + a_oo) @ (z * np.eye(size) + a_oo), a_oe
        )
        self_energy[str(z)] = float(np.linalg.norm(sigma, 2))
        resolvent_residual[str(z)] = float(np.linalg.norm(reduced - exact, 2))
        derivative_norm[str(z)] = float(np.linalg.norm(sigma_prime, 2))

    times = [0.1, 0.25, 0.5, 1.0]
    heat_defect, wave_defect = {}, {}
    static_heat_error, static_wave_error = {}, {}
    factorization_residual_heat, factorization_residual_wave = {}, {}
    selector = np.zeros((size, N))
    injector = np.zeros((N, size))
    for j, x in enumerate(retained):
        selector[j, x] = 1.0
        injector[x, j] = 1.0
    q_hidden = np.eye(N) - injector @ selector
    for time in times:
        heat = expm(-time * a)
        heat2 = expm(-2 * time * a)
        h_e = selector @ heat @ injector
        h_e2 = selector @ heat2 @ injector
        h_def = h_e2 - h_e @ h_e
        h_fact = selector @ heat @ q_hidden @ heat @ injector
        heat_defect[str(time)] = float(np.linalg.norm(h_def, 2))
        factorization_residual_heat[str(time)] = float(np.linalg.norm(h_def - h_fact, 2))
        static_heat_error[str(time)] = float(
            np.linalg.norm(h_e - expm(-time * static), 2)
        )

        wave = expm(-1j * time * a)
        wave2 = expm(-2j * time * a)
        u_e = selector @ wave @ injector
        u_e2 = selector @ wave2 @ injector
        u_def = u_e2 - u_e @ u_e
        u_fact = selector @ wave @ q_hidden @ wave @ injector
        wave_defect[str(time)] = float(np.linalg.norm(u_def, 2))
        factorization_residual_wave[str(time)] = float(np.linalg.norm(u_def - u_fact, 2))
        static_wave_error[str(time)] = float(
            np.linalg.norm(u_e - expm(-1j * time * static), 2)
        )

    return {
        "retained_nodes": retained.tolist(),
        "hidden_nodes": hidden.tolist(),
        "static_schur_row_sum_residual": float(np.max(np.abs(static.sum(axis=1)))),
        "static_schur_eigenvalues": np.linalg.eigvalsh(static).tolist(),
        "self_energy_operator_norm": self_energy,
        "self_energy_derivative_operator_norm": derivative_norm,
        "block_resolvent_identity_residual": resolvent_residual,
        "heat_semigroup_compression_defect": heat_defect,
        "wave_group_compression_defect": wave_defect,
        "static_schur_heat_propagator_error": static_heat_error,
        "static_schur_wave_propagator_error": static_wave_error,
        "exact_hidden_excursion_identity_residual_heat": factorization_residual_heat,
        "exact_hidden_excursion_identity_residual_wave": factorization_residual_wave,
        "interpretation": (
            "Eliminating odd nodes generates a nonconstant frequency-dependent "
            "self-energy and therefore exact reduced memory. This is a typed "
            "mathematical object, not a physical environment claim."
        ),
    }


def make_figures(
    singles: list[dict[str, object]], numerical: dict[str, object]
) -> None:
    matrix = np.zeros((len(CHARACTERISTICS), len(DOMAINS)))
    for row in singles:
        i = int(str(row["characteristic_id"])[1:]) - 1
        j = int(str(row["domain_id"])[1:]) - 1
        matrix[i, j] = int(row["structural_score_0_4"])
    fig, ax = plt.subplots(figsize=(13.5, 7.2))
    image = ax.imshow(matrix, vmin=0, vmax=4, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(DOMAINS)), [d[0] for d in DOMAINS], rotation=45, ha="right")
    ax.set_yticks(range(len(CHARACTERISTICS)), [c[0] for c in CHARACTERISTICS])
    ax.set_xlabel("comparison domain")
    ax.set_ylabel("FIN puzzle")
    ax.set_title("Heuristic structural-fit atlas (discovery aid, not equivalence evidence)")
    fig.colorbar(image, ax=ax, label="structural score 0–4")
    fig.tight_layout()
    fig.savefig(FIG / "single_puzzle_domain_atlas.png", dpi=190)
    plt.close(fig)

    times = np.array([0.1, 0.25, 0.5, 1.0])
    heat = [numerical["heat_semigroup_compression_defect"][str(v)] for v in times]
    wave = [numerical["wave_group_compression_defect"][str(v)] for v in times]
    sheat = [numerical["static_schur_heat_propagator_error"][str(v)] for v in times]
    swave = [numerical["static_schur_wave_propagator_error"][str(v)] for v in times]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.7))
    axes[0].plot(times, heat, "o-", label="heat")
    axes[0].plot(times, wave, "s-", label="wave")
    axes[0].set_title("Projected composition defect")
    axes[0].set_xlabel("dimensionless time")
    axes[0].set_ylabel("operator norm")
    axes[0].grid(alpha=0.25)
    axes[0].legend()
    axes[1].plot(times, sheat, "o-", label="heat vs static Schur")
    axes[1].plot(times, swave, "s-", label="wave vs static Schur")
    axes[1].set_title("Static effective-generator error")
    axes[1].set_xlabel("dimensionless time")
    axes[1].set_ylabel("operator norm")
    axes[1].grid(alpha=0.25)
    axes[1].legend()
    fig.suptitle("Hidden-node elimination creates exact reduced memory on strict Z12")
    fig.tight_layout()
    fig.savefig(FIG / "dynamic_schur_memory.png", dpi=190)
    plt.close(fig)


def main() -> None:
    singles, pairs = write_atlas_tables()
    numerical = dynamic_schur_study()
    ranked_pairs = sorted(
        (
            row
            for row in pairs
            if row["touches_both_puzzles"] and int(row["structural_score_0_4"]) >= 3
        ),
        key=lambda row: (
            -int(row["structural_score_0_4"]),
            str(row["pair_id"]),
            str(row["domain_id"]),
        ),
    )
    payload = {
        "warning": "Tag scores are heuristic search indices, not scientific evidence.",
        "coverage": {
            "single_puzzles": len(CHARACTERISTICS),
            "domains": len(DOMAINS),
            "single_domain_comparisons": len(singles),
            "unordered_puzzle_pairs": math.comb(len(CHARACTERISTICS), 2),
            "pair_domain_comparisons": len(pairs),
        },
        "dynamic_schur_study": numerical,
        "top_pair_domain_candidates": ranked_pairs[:30],
    }
    (HERE / "FIN_Nadsoliton_Puzzle_Atlas_Results.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    make_figures(singles, numerical)
    print(json.dumps(payload["coverage"], indent=2))
    print(
        "max block-resolvent residual:",
        max(numerical["block_resolvent_identity_residual"].values()),
    )
    print(
        "min nonzero self-energy derivative norm:",
        min(numerical["self_energy_derivative_operator_norm"].values()),
    )


if __name__ == "__main__":
    main()
