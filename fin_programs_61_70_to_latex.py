#!/usr/bin/env python3
"""Build archival LaTeX for FIN Programs 61--70."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Programs_61_70_Continuum_Operational_Physics_Monograph.md")
TARGET = Path("FIN_Programs_61_70_Continuum_Operational_Physics_Monograph.tex")

PREAMBLE = r"""\documentclass[11pt,a4paper,openany]{report}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
\usepackage{microtype}
\usepackage{amsmath,amssymb,amsthm,mathtools,bm}
\usepackage{booktabs,array,longtable,tabularx}
\usepackage{xcolor}
\usepackage[margin=24mm,headheight=23pt]{geometry}
\usepackage{enumitem}
\usepackage{fancyhdr}
\usepackage{fancyvrb}
\usepackage{graphicx}
\usepackage{url}
\usepackage{hyperref}

\definecolor{finblue}{HTML}{1F5A99}
\definecolor{fingreen}{HTML}{19733A}
\definecolor{finorange}{HTML}{D55E00}
\definecolor{finred}{HTML}{A61B1B}
\definecolor{finviolet}{HTML}{6A3D9A}

\hypersetup{
  unicode=true,
  colorlinks=true,
  linkcolor=finblue,
  citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Programs 61--70: Continuum and Operational Physics},
  pdfauthor={Krzysztof \.Zuchowski},
  pdfsubject={Continuum tests and falsifiable information-to-experiment interfaces},
  pdfkeywords={FIN, Schur compression, continuum, tomography, causality, Landauer}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs 61--70}
\fancyhead[R]{\small Release 10.7}
\fancyfoot[C]{\thepage}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.55em}
\setlist{nosep,leftmargin=2em}
\setcounter{tocdepth}{2}
\setcounter{secnumdepth}{3}
\emergencystretch=2em

\newcommand{\one}{\bm 1}
\newcommand{\C}{\mathbb C}
\newcommand{\R}{\mathbb R}
\newcommand{\Z}{\mathbb Z}
\newcommand{\statusProven}{\textcolor{fingreen}{[Proven]}}
\newcommand{\statusStrong}{\textcolor{finblue}{[Strong evidence]}}
\newcommand{\statusModerate}{\textcolor{finviolet}{[Moderate evidence]}}
\newcommand{\statusConditional}{\textcolor{finorange}{[Conditional]}}
\newcommand{\statusRefuted}{\textcolor{finred}{[Refuted]}}

\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{9mm}
{\Large\bfseries FIN Research Monograph --- Release 10.7\par}
\vspace{11mm}
{\Huge\bfseries FIN Programs 61--70\par}
\vspace{6mm}
{\Large Continuum Tests, Operational Physics, and\\Falsifiable Information-to-Experiment Interfaces\par}
\vspace{17mm}
{\Large Krzysztof \.Zuchowski\par}
\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
\vspace{2mm}
{\normalsize ORCID: \href{https://orcid.org/0009-0002-0909-3613}{0009-0002-0909-3613}\par}
\vspace{10mm}
{\large 27 July 2026\par}
{\normalsize Publication --- Preprint; record version 1.0.0\par}
\vfill
\begin{minipage}{0.88\textwidth}
\small
\textbf{Scope.} Ten executed programs on projective continuum compatibility,
regularizer independence, compression composition, signed/Krein structure,
operational tomography, causal tails, calibration, chirality, Landauer
thermodynamics, and blinded model scoring.

\medskip
\textbf{Guardrail.} No silent legacy--strict identity, strict selector,
physical-unit source, causal/Lorentz closure, role transfer, or
Theory-of-Everything closure is claimed.

\medskip
\textbf{Reproducibility.}
\path{fin_programs_61_70_continuum_operational_physics.py},
\path{test_fin_programs_61_70_continuum_operational_physics.py}, JSON, and figures.

\medskip
\textbf{License.} Creative Commons Attribution 4.0 International (CC BY 4.0).
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}

\pagenumbering{roman}
\tableofcontents
\clearpage
\pagenumbering{arabic}
"""

POSTAMBLE = r"""
\clearpage
\chapter*{Archival note}
\addcontentsline{toc}{chapter}{Archival note}
Machine-readable results are stored in
\path{FIN_Programs_61_70_Continuum_Operational_Physics_Results.json}. The
executable and ten regression tests reproduce the finite results and figures.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Programs 61--70: Continuum Tests,
Operational Physics, and Falsifiable Information-to-Experiment Interfaces}
(FIN Research Monograph, Release 10.7; Version 1.0.0) [Preprint].
\end{document}
"""

FIGURES = {
    "FINFIG61": (
        "FIN_Programs_61_70_Continuum_Operational_Physics_Figures/program61_projective_continuum.png",
        "Projective naturality defects under fixed mean-diagonal normalization.",
    ),
    "FINFIG62": (
        "FIN_Programs_61_70_Continuum_Operational_Physics_Figures/program62_regularizer_limit.png",
        r"Shifted-resolvent resistance converges as \(\delta\to0^+\) while the full Green norm diverges.",
    ),
    "FINFIG65": (
        "FIN_Programs_61_70_Continuum_Operational_Physics_Figures/program65_operational_tomography.png",
        "Single-preparation distinguishability of unitary and diffusive record laws.",
    ),
    "FINFIG66": (
        "FIN_Programs_61_70_Continuum_Operational_Physics_Figures/program66_causal_tails.png",
        "Opposite-site tails for the full strict kernel and a nearest-neighbour control.",
    ),
    "FINFIG69": (
        "FIN_Programs_61_70_Continuum_Operational_Physics_Figures/program69_landauer_protocol.png",
        "Reversible work, bath heat, and reset error in the explicit erasure protocol.",
    ),
    "FINFIG70": (
        "FIN_Programs_61_70_Continuum_Operational_Physics_Figures/program70_blinded_challenge.png",
        "Held-out synthetic model ranking with the exact hidden generator excluded.",
    ),
}


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    replacements = {
        "![Projective continuum](FIN_Programs_61_70_Continuum_Operational_Physics_Figures/program61_projective_continuum.png)": "FINFIG61",
        "![Regularizer limit](FIN_Programs_61_70_Continuum_Operational_Physics_Figures/program62_regularizer_limit.png)": "FINFIG62",
        "![Operational tomography](FIN_Programs_61_70_Continuum_Operational_Physics_Figures/program65_operational_tomography.png)": "FINFIG65",
        "![Causal tails](FIN_Programs_61_70_Continuum_Operational_Physics_Figures/program66_causal_tails.png)": "FINFIG66",
        "![Landauer protocol](FIN_Programs_61_70_Continuum_Operational_Physics_Figures/program69_landauer_protocol.png)": "FINFIG69",
        "![Blinded challenge](FIN_Programs_61_70_Continuum_Operational_Physics_Figures/program70_blinded_challenge.png)": "FINFIG70",
    }
    for markdown, token in replacements.items():
        text = text.replace(markdown, f"\n\n{token}\n\n")
    text = re.sub(r"```[a-zA-Z]*\n(.*?)```", r"\1", text, flags=re.DOTALL)
    text = re.sub(r"`([^`]+)`", r"\1", text)
    text = re.sub(r"^(#{1,4})\s+\d+(?:\.\d+)*\.\s+", r"\1 ", text, flags=re.MULTILINE)

    body = converter.body_from_markdown(text)
    for token, (path, caption) in FIGURES.items():
        figure = (
            "\\begin{figure}[htbp]\n"
            "\\centering\n"
            f"\\includegraphics[width=0.92\\textwidth]{{{path}}}\n"
            f"\\caption{{{caption}}}\n"
            "\\end{figure}\n"
        )
        body = body.replace(token, figure)
    for plain, macro in [
        ("[Proven]", r"\statusProven{}"),
        ("[Strong evidence]", r"\statusStrong{}"),
        ("[Moderate evidence]", r"\statusModerate{}"),
        ("[Conditional]", r"\statusConditional{}"),
        ("[Refuted]", r"\statusRefuted{}"),
    ]:
        body = body.replace(plain, macro)
    for filename in [
        "fin_programs_61_70_continuum_operational_physics.py",
        "test_fin_programs_61_70_continuum_operational_physics.py",
        "FIN_Programs_61_70_Continuum_Operational_Physics_Results.json",
        "FIN_Programs_61_70_Continuum_Operational_Physics_Figures/",
    ]:
        body = body.replace(converter.escape_normal(filename), f"\\path{{{filename}}}")
    digest = "94b3bef894ce659a1c42e2783f7a478846c0114ea45d754e33f3fe7712f641bf"
    digest_tex = "\\texttt{" + "}\\allowbreak\\texttt{".join(
        digest[i : i + 16] for i in range(0, len(digest), 16)
    ) + "}"
    body = body.replace(digest, digest_tex)
    body = body.replace(
        f"Hidden generator digest: {digest_tex}.",
        "Hidden generator digest:\\par\\smallskip\\noindent"
        f"\\texttt{{{digest[:32]}}}\\par\\noindent"
        f"\\texttt{{{digest[32:]}}}.",
    )
    TARGET.write_text(PREAMBLE + "\n" + body + "\n" + POSTAMBLE, encoding="utf-8")
    print(f"Wrote {TARGET} ({TARGET.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
