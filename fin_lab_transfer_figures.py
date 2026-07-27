#!/usr/bin/env python3
"""Generate code-native figures for the FIN laboratory transfer package."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np


ROOT = Path(__file__).resolve().parent
FIG = ROOT / "FIN_Lab_P240_242_Transfer_Figures"

BLUE = "#1F5A99"
GREEN = "#19733A"
ORANGE = "#D55E00"
RED = "#A61B1B"
VIOLET = "#6A3D9A"
GREY = "#EDF1F5"
DARK = "#1B2631"


def box(ax, xy, width, height, text, color, fontsize=9, text_color="white"):
    text = text.replace("\\n", "\n")
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.025,rounding_size=0.03",
        facecolor=color,
        edgecolor="white",
        linewidth=1.5,
    )
    ax.add_patch(patch)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        text,
        ha="center",
        va="center",
        color=text_color,
        fontsize=fontsize,
        wrap=True,
    )
    return patch


def arrow(ax, start, end, color=DARK, label=None, yoffset=0.0):
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=13,
        linewidth=1.5,
        color=color,
    )
    ax.add_patch(patch)
    if label:
        ax.text(
            (start[0] + end[0]) / 2,
            (start[1] + end[1]) / 2 + yoffset,
            label,
            ha="center",
            va="center",
            fontsize=8,
            color=color,
            bbox=dict(facecolor="white", edgecolor="none", pad=1.5),
        )


def duality_map():
    fig, ax = plt.subplots(figsize=(12.0, 6.0), constrained_layout=True)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6)
    ax.axis("off")
    box(
        ax,
        (4.25, 4.65),
        3.5,
        0.8,
        r"strict finite target\n$W=K_s,\ s=1.660307278766099$\n$A=sI-W\geq0$",
        BLUE,
        fontsize=10,
    )
    box(ax, (0.65, 2.75), 3.2, 0.9, r"coherent ray\n$U_t=e^{-itA}$\nBorn populations $|U_t|^2$", VIOLET)
    box(ax, (8.15, 2.75), 3.2, 0.9, r"diffusive ray\n$P_t=e^{-tA}$\nreversible Markov channel", GREEN)
    arrow(ax, (4.7, 4.65), (3.1, 3.65), label="choose quantum state cone", yoffset=0.18)
    arrow(ax, (7.3, 4.65), (8.9, 3.65), label="choose probability simplex", yoffset=0.18)
    box(
        ax,
        (3.05, 0.65),
        5.9,
        1.05,
        "operational completion\nstate + preparation + calibrated clock + channel + instrument/POVM\n"
        "+ environment/apparatus + event record + custody/holdout",
        ORANGE,
        fontsize=9.5,
    )
    arrow(ax, (2.25, 2.75), (4.35, 1.7), label="records", yoffset=0.15)
    arrow(ax, (9.75, 2.75), (7.65, 1.7), label="records", yoffset=0.15)
    ax.text(
        6,
        0.17,
        "The operator fixes spatial spectral structure; it does not select the physical temporal category.",
        ha="center",
        fontsize=10,
        color=RED,
        weight="bold",
    )
    ax.set_title("FIN duality and the missing operational object", fontsize=15, weight="bold")
    fig.savefig(FIG / "duality_operational_map.png", dpi=220)
    plt.close(fig)


def platform_matrix():
    labels = [
        "1. complete\npreparations",
        "2. calibrated\nclock",
        "3. resolving\nmeasurement",
        "4. finite raw\ncounts",
        "5. independent\ncustody",
    ]
    platforms = [
        "P-A 12-mode photonic / stochastic process",
        "P-B qubit Ramsey / M-Z auxiliary",
        "P-C event-level double slit",
    ]
    values = np.array(
        [
            [1.0, 1.0, 1.0, 1.0, 0.25],
            [0.2, 1.0, 0.2, 1.0, 0.25],
            [0.25, 1.0, 0.35, 1.0, 0.55],
        ]
    )
    fig, ax = plt.subplots(figsize=(11.6, 4.8), constrained_layout=True)
    image = ax.imshow(values, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(labels)), labels)
    ax.set_yticks(range(len(platforms)), platforms)
    ax.tick_params(axis="x", labelsize=9)
    ax.tick_params(axis="y", labelsize=9)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            if j == 4:
                text = "human gate"
            elif values[i, j] >= 0.9:
                text = "direct"
            elif values[i, j] >= 0.5:
                text = "conditional"
            else:
                text = "not P227"
            ax.text(
                j,
                i,
                text,
                ha="center",
                va="center",
                fontsize=8,
                color="black",
                weight="bold" if values[i, j] >= 0.9 else "normal",
            )
    ax.set_title(
        "Platform fit to the five measurement obligations (feasibility hypotheses)",
        fontsize=13,
        weight="bold",
    )
    colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
    colorbar.set_label("direct match to P227/P241 duty")
    fig.savefig(FIG / "platform_obligation_matrix.png", dpi=220)
    plt.close(fig)


def photonic_architecture():
    fig, ax = plt.subplots(figsize=(12.4, 6.2), constrained_layout=True)
    ax.set_xlim(0, 12.4)
    ax.set_ylim(0, 6.2)
    ax.axis("off")
    box(ax, (0.35, 4.65), 2.0, 0.85, "heralded SPDC\nor characterized\nattenuated laser", BLUE)
    box(ax, (2.8, 4.65), 2.0, 0.85, r"12 selectable" "\n" r"basis inputs" "\n" r"$j=0,\ldots,11$", BLUE)
    arrow(ax, (2.35, 5.075), (2.8, 5.075))

    box(ax, (5.25, 4.65), 2.65, 0.85, "shared calibrated target\ncoupling matrix $A$\n(not just a visual fit)", ORANGE)
    arrow(ax, (4.8, 5.075), (5.25, 5.075))

    box(ax, (4.15, 2.65), 3.0, 0.9, "coherent implementation\nlaser-written guides or\nprogrammable M-Z mesh", VIOLET)
    box(ax, (8.2, 2.65), 3.25, 0.9, "dissipative implementation\nverified jumps $L_{xy}$ or\nprogrammed stochastic routing", GREEN)
    arrow(ax, (6.0, 4.65), (5.65, 3.55), label=r"$U_\tau$", yoffset=0.1)
    arrow(ax, (7.1, 4.65), (9.65, 3.55), label=r"$P_\tau$", yoffset=0.1)

    box(ax, (4.15, 0.75), 3.0, 0.85, "12-channel SPAD/SNSPD\n+ TCSPC; phase-sensitive\ninterferometric controls", VIOLET)
    box(ax, (8.2, 0.75), 3.25, 0.85, "12-channel event counts\nat $\tau$ and sealed $2\tau$\n+ calibration/controls", GREEN)
    arrow(ax, (5.65, 2.65), (5.65, 1.6))
    arrow(ax, (9.8, 2.65), (9.8, 1.6))
    ax.text(
        1.4,
        2.0,
        "Critical falsification:\nordinary loss/dephasing\nneed not equal\n$e^{-\\tau A}$.",
        ha="center",
        va="center",
        fontsize=11,
        color=RED,
        weight="bold",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#FFF4F1", edgecolor=RED),
    )
    arrow(ax, (2.35, 2.0), (4.15, 2.95), color=RED)
    ax.set_title(
        "Recommended direct platform P-A: paired coherent and heat-process realization",
        fontsize=14,
        weight="bold",
    )
    fig.savefig(FIG / "platform_A_apparatus.png", dpi=220)
    plt.close(fig)


def double_slit_architecture():
    fig, ax = plt.subplots(figsize=(12.0, 5.2), constrained_layout=True)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 5.2)
    ax.axis("off")
    box(ax, (0.35, 2.9), 2.0, 1.0, "single-photon source\nSPDC herald or\ncharacterized weak source", BLUE)
    box(ax, (3.0, 2.9), 1.85, 1.0, "two slits\nmotorized shutters\n4 configurations", ORANGE)
    box(ax, (5.55, 2.9), 2.0, 1.0, "free propagation\nfixed geometry\nstability monitors", GREY, text_color=DARK)
    box(ax, (8.25, 2.9), 2.45, 1.0, "position-sensitive\nSPAD array + TCSPC\n(or gated ICCD/EMCCD)", GREEN)
    arrow(ax, (2.35, 3.4), (3.0, 3.4))
    arrow(ax, (4.85, 3.4), (5.55, 3.4))
    arrow(ax, (7.55, 3.4), (8.25, 3.4))
    box(
        ax,
        (3.0, 0.65),
        7.7,
        1.05,
        "RAW event row: event_id, UTC timestamp, detector x/y, configuration,\n"
        "run_id, subset, shutter intervention, calibration reference\n"
        "Required configurations: both_open, left_only, right_only, both_closed",
        VIOLET,
        fontsize=9.5,
    )
    arrow(ax, (9.45, 2.9), (8.8, 1.7), label="event stream", yoffset=0.1)
    ax.text(
        6,
        0.18,
        "A rendered interference image is not admissible evidence. This platform does not by itself reconstruct the 12x12 FIN generator.",
        ha="center",
        fontsize=9.5,
        color=RED,
        weight="bold",
    )
    ax.set_title(
        "Platform P-C: event-level double-slit custody pilot", fontsize=14, weight="bold"
    )
    fig.savefig(FIG / "platform_C_double_slit.png", dpi=220)
    plt.close(fig)


def custody_pipeline():
    fig, ax = plt.subplots(figsize=(12.2, 5.4), constrained_layout=True)
    ax.set_xlim(0, 12.2)
    ax.set_ylim(0, 5.4)
    ax.axis("off")
    box(ax, (0.3, 3.05), 2.2, 1.05, "PROVIDER\nacquires raw events\n+ calibration/controls", BLUE)
    box(ax, (3.2, 3.05), 2.2, 1.05, "REGISTRAR\nverifies SHA-256\nseals holdout, signs", ORANGE)
    box(ax, (6.1, 3.05), 2.2, 1.05, "ANALYST\nlocks P240/P242\nbefore unblinding", VIOLET)
    box(ax, (9.0, 3.05), 2.7, 1.05, "P242 ONE SHOT\natomically consumes token\nreports pass/fail/controls", GREEN)
    arrow(ax, (2.5, 3.58), (3.2, 3.58), label="raw bundle")
    arrow(ax, (5.4, 3.58), (6.1, 3.58), label="calibration only")
    arrow(ax, (8.3, 3.58), (9.0, 3.58), label="signed unblinding")
    box(
        ax,
        (1.2, 0.7),
        3.2,
        1.0,
        "P241 validator\n11/11 schema + chronology\n+ trusted registrar signature",
        GREY,
        text_color=DARK,
    )
    box(
        ax,
        (7.2, 0.7),
        3.2,
        1.0,
        "immutable ledger\nbundle digest + analysis hash\n+ execution receipt",
        GREY,
        text_color=DARK,
    )
    arrow(ax, (4.3, 3.05), (3.0, 1.7))
    arrow(ax, (10.35, 3.05), (8.8, 1.7))
    ax.text(
        6.1,
        0.15,
        "provider != registrar != analyst; no threshold repair and no second run after seeing the holdout",
        ha="center",
        fontsize=10,
        color=RED,
        weight="bold",
    )
    ax.set_title("Blind-custody transfer and one-execution rule", fontsize=14, weight="bold")
    fig.savefig(FIG / "custody_pipeline.png", dpi=220)
    plt.close(fig)


def main():
    FIG.mkdir(exist_ok=True)
    duality_map()
    platform_matrix()
    photonic_architecture()
    double_slit_architecture()
    custody_pipeline()
    print(FIG)


if __name__ == "__main__":
    main()
