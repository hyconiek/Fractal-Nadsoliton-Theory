#!/usr/bin/env python3
"""Build archival LaTeX for FIN Programs 125--137."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Programs_125_137_Trace_Localizer_Physics_Monograph.md")
TARGET = Path("FIN_Programs_125_137_Trace_Localizer_Physics_Monograph.tex")

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
\usepackage{xurl}
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
  pdftitle={FIN Programs 125--137: Trace Selection, Natural Localizers, Fractional Physics, and Operational Sources},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Invariant trace obstruction, natural C12 localizer, fractional dynamics, dimensional calibration, apparatus tomography, and signed operational sources},
  pdfkeywords={FIN, trace selection, cyclic groups, fractional Laplacian, stable process, wave dynamics, operational physics, apparatus tomography}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs 125--137}
\fancyhead[R]{\small Release 10.13}
\fancyfoot[C]{\thepage}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.55em}
\setlist{nosep,leftmargin=2em}
\setcounter{tocdepth}{2}
\setcounter{secnumdepth}{3}
\emergencystretch=3em

\newcommand{\one}{\bm 1}
\newcommand{\C}{\mathbb C}
\newcommand{\R}{\mathbb R}
\newcommand{\Z}{\mathbb Z}
\newcommand{\statusProven}{\textcolor{fingreen}{[Proven]}}
\newcommand{\statusComputer}{\textcolor{fingreen}{[Proven, computer-assisted]}}
\newcommand{\statusStrong}{\textcolor{finblue}{[Strong evidence]}}
\newcommand{\statusConditional}{\textcolor{finorange}{[Conditional]}}
\newcommand{\statusRefuted}{\textcolor{finred}{[Refuted in scope]}}
\newcommand{\statusOpen}{\textcolor{finviolet}{[Open]}}

\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{9mm}
{\Large\bfseries FIN Research Monograph --- Release 10.13\par}
\vspace{11mm}
{\Huge\bfseries FIN Programs 125--137\par}
\vspace{6mm}
{\Large Trace Selection, Natural Localizers,\\
Fractional Physics, and Operational Sources\par}
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
\textbf{Scope.} Thirteen executed studies of invariant trace selection,
the natural \(C_{12}\) localizer, quantitative fractional limits,
wave ultraviolet resolution, physical calibration, apparatus tomography,
parameter-source falsification, conditional damping completion, and signed
operational states.

\medskip
\textbf{Central result.} The localized fibre is natural, but its invariant
probability traces form a simplex. The uniform trace gives \(9/5\), while the
normalized Hilbert trace gives \(17/9\). Present symmetry does not select the
strict damping exponent.

\medskip
\textbf{Guardrail.} No strict selector, internal physical units, completed
legacy--strict bridge, role transfer, role-bearing \(L_{\mathrm{total}}\),
Theory-of-Everything closure, or external validation is claimed.

\medskip
\textbf{License.} Creative Commons Attribution 4.0 International (CC BY 4.0).
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}

\pagenumbering{roman}
\chapter*{Publication statement}
\addcontentsline{toc}{chapter}{Publication statement}
This release reports Programs 125--137 and recommends Programs 138--150.
Proofs, computer-assisted calculations, strong numerical evidence,
conditional constructions, falsified candidates, and open obligations are
kept separate.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""

POSTAMBLE = r"""
\clearpage
\chapter*{Archival note}
\addcontentsline{toc}{chapter}{Archival note}
Machine-readable results are stored in
\path{FIN_Programs_125_137_Trace_Localizer_Physics_Results.json}.
The executable, thirty-five tests, and twelve generated figures reproduce the
release.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Programs 125--137: Trace Selection,
Natural Localizers, Fractional Physics, and Operational Sources}
(FIN Research Monograph, Release 10.13; Version 1.0.0) [Preprint]. Zenodo.
\end{document}
"""

FIGURES = {
    "program125_invariant_trace_simplex.png": "Invariant trace simplex. Symmetry leaves two independent state parameters and does not select the uniform sector trace.",
    "program126_natural_fibre_localizer.png": r"Natural decomposition of the \(C_{12}\) homological--character fibres.",
    "program127_continuous_fractional_enclosure.png": r"Continuous finite-window enclosure around the \(C|q|^{4/5}\) symbol law.",
    "program128_quantitative_stable_window.png": r"Finite-step characteristic-function distance bounds in the certified frequency window.",
    "program129_fractional_wave_uv_obstruction.png": r"The standard dyadic wave bound grows as \(\Lambda^{3/5}\) and is not ultraviolet summable.",
    "program130_dimensional_calibration.png": "The minimal physical calibration map has rank three.",
    "program131_apparatus_tomography.png": "Held-out log loss for nonparametric apparatus-memory models.",
    "program132_crossover_rg_flow.png": "Exact projective flow from the local ultraviolet endpoint to the fractional infrared endpoint.",
    "program133_phase_frequency_source_test.png": "Exact arithmetic parameter encodings fail robustness and source-selection tests.",
    "program134_amplitude_projectivization.png": r"Amplitude projectivization removes \(\alpha_{\rm geo}\) but leaves a large shape mismatch.",
    "program135_conditional_damping_bridge.png": "The trace-dependent family of conditional damping completion factors.",
    "program136_signed_state_receiver.png": "The signed receiver distinguishes prepared opposite branches but does not choose a branch.",
}


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    text = re.sub(
        r"\$\$\s*(.*?)\s*\$\$",
        lambda match: "\\[\n" + match.group(1).strip() + "\n\\]",
        text,
        flags=re.DOTALL,
    )
    text = re.sub(
        r"(?<!\\)\$([^\n$]+?)(?<!\\)\$",
        lambda match: "\\(" + match.group(1) + "\\)",
        text,
    )
    text = (
        text.replace("–", "--")
        .replace("—", "---")
        .replace("’", "'")
        .replace("“", '"')
        .replace("”", '"')
    )
    for filename in FIGURES:
        pattern = re.compile(
            r"!\[[^\]]*\]\(FIN_Programs_125_137_Trace_Localizer_Physics_Figures/"
            + re.escape(filename)
            + r"\)"
        )
        text = pattern.sub(f"\n\nFINFIG_{filename}\n\n", text)
    text = re.sub(r"```[a-zA-Z]*\n(.*?)```", r"\1", text, flags=re.DOTALL)
    text = re.sub(r"`([^`]+)`", r"\1", text)
    text = re.sub(
        r"^(#{1,4})\s+\d+(?:\.\d+)*(?:\.)?\s+",
        r"\1 ",
        text,
        flags=re.MULTILINE,
    )
    body = converter.body_from_markdown(text)
    body = body.replace("Ż", r"\.Z").replace("ż", r"\.z")
    for filename, caption in FIGURES.items():
        token = f"FINFIG\\_{filename.replace('_', r'\_')}"
        path = (
            "FIN_Programs_125_137_Trace_Localizer_Physics_Figures/"
            + filename
        )
        body = body.replace(
            token,
            "\\begin{figure}[htbp]\n"
            "\\centering\n"
            f"\\includegraphics[width=0.94\\textwidth]{{{path}}}\n"
            f"\\caption{{{caption}}}\n"
            "\\end{figure}",
        )
    for plain, macro in [
        ("[Proven, computer-assisted]", r"\statusComputer{}"),
        ("[Proven]", r"\statusProven{}"),
        ("[Strong evidence]", r"\statusStrong{}"),
        ("[Conditional]", r"\statusConditional{}"),
        ("[Refuted in scope]", r"\statusRefuted{}"),
        ("[Open]", r"\statusOpen{}"),
    ]:
        body = body.replace(plain, macro)
    TARGET.write_text(PREAMBLE + "\n" + body + "\n" + POSTAMBLE, encoding="utf-8")
    print(TARGET)


if __name__ == "__main__":
    main()
