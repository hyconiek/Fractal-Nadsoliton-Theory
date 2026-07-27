#!/usr/bin/env python3
"""Build archival LaTeX for FIN Programs 91--100."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Programs_91_100_Critical_Nonlocal_Operational_Monograph.md")
TARGET = Path("FIN_Programs_91_100_Critical_Nonlocal_Operational_Monograph.tex")

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
  pdftitle={FIN Programs 91--100: Critical Scaling, Nonlocal Continuum Structure, and Operational Completion Tests},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Schur nonclosure, critical scaling, nonlocal continuum limits, process tensors, and damping completion},
  pdfkeywords={FIN, Schur complement, renormalization, graphon, nonlocal operator, process tensor, information thermodynamics}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs 91--100}
\fancyhead[R]{\small Release 10.10}
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
\newcommand{\statusModerate}{\textcolor{finviolet}{[Moderate evidence]}}
\newcommand{\statusConditional}{\textcolor{finorange}{[Conditional]}}
\newcommand{\statusRefuted}{\textcolor{finred}{[Refuted]}}
\newcommand{\statusOpen}{\textcolor{finviolet}{[Open]}}

\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{9mm}
{\Large\bfseries FIN Research Monograph --- Release 10.10\par}
\vspace{11mm}
{\Huge\bfseries FIN Programs 91--100\par}
\vspace{6mm}
{\Large Critical Scaling, Nonlocal Continuum Structure,\\
and Operational Completion Tests\par}
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
\textbf{Scope.} Ten executed studies of strict-kernel Schur nonclosure,
critical scaling, graphon continuum limits, long-range propagation,
adaptive quotient structure, selector-source provenance, process-tensor
design, feedback thermodynamics, data acquisition, and the
legacy--strict damping-completion atom.

\medskip
\textbf{Central theorem.} A tail-controlled scale-invariant spectral ratio
proves that the fixed-mass strict lattice family is not closed, even up to
scalar normalization, under alternating-site Schur reduction.

\medskip
\textbf{Guardrail.} No strict selector, physical-unit source, Lorentz closure,
legacy role transfer, completed legacy--strict bridge, role-bearing
\(L_{\mathrm{total}}\), or external physical validation is claimed.

\medskip
\textbf{License.} Creative Commons Attribution 4.0 International (CC BY 4.0).
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}

\pagenumbering{roman}
\chapter*{Publication statement}
\addcontentsline{toc}{chapter}{Publication statement}
This release reports Programs 91--100 and recommends Programs 101--112.
Theorem-level statements, computer-assisted certificates, numerical evidence,
conditioned operational models, failed inferences, and open obligations are
labeled separately.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""

POSTAMBLE = r"""
\clearpage
\chapter*{Archival note}
\addcontentsline{toc}{chapter}{Archival note}
Machine-readable results are stored in
\path{FIN_Programs_91_100_Critical_Nonlocal_Operational_Results.json}.
The executable, seventeen tests, seven figures, and external-data intake
template reproduce the release.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Programs 91--100: Critical Scaling,
Nonlocal Continuum Structure, and Operational Completion Tests}
(FIN Research Monograph, Release 10.10; Version 1.0.0) [Preprint].
\end{document}
"""

FIGURES = {
    "FINFIG91": (
        "FIN_Programs_91_100_Critical_Nonlocal_Operational_Figures/program91_strict_nonclosure_certificate.png",
        "Tail-certified separation of the native and Schur scale-invariant spectral ratios.",
    ),
    "FINFIG92": (
        "FIN_Programs_91_100_Critical_Nonlocal_Operational_Figures/program92_critical_tuning.png",
        r"Exact nearest-neighbour critical-tuning defect and its \(N^{-2}\) relative rate.",
    ),
    "FINFIG93": (
        "FIN_Programs_91_100_Critical_Nonlocal_Operational_Figures/program93_graphon_continuum.png",
        "First-order discretization convergence and bounded high-mode continuum spectrum.",
    ),
    "FINFIG94": (
        "FIN_Programs_91_100_Critical_Nonlocal_Operational_Figures/program94_long_range_propagation.png",
        "Certified coupling-tail bounds for the strict long-range lattice profile.",
    ),
    "FINFIG97": (
        "FIN_Programs_91_100_Critical_Nonlocal_Operational_Figures/program97_process_tensor_design.png",
        "Maximin two-time instrument design under symmetric detector confusion.",
    ),
    "FINFIG98": (
        "FIN_Programs_91_100_Critical_Nonlocal_Operational_Figures/program98_feedback_ledger.png",
        "System feedback work and apparatus-inclusive work after memory reset.",
    ),
    "FINFIG100": (
        "FIN_Programs_91_100_Critical_Nonlocal_Operational_Figures/program100_damping_completion_atom.png",
        r"Legacy and strict damping envelopes and their necessary \(d^{-4/5}\) completion multiplier.",
    ),
}


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    text = re.sub(
        r"\$\$\s*(.*?)\s*\$\$",
        lambda match: "\\[\n" + match.group(1).strip() + "\n\\]",
        text,
        flags=re.DOTALL,
    )
    text = (
        text.replace("–", "--")
        .replace("—", "---")
        .replace("’", "'")
        .replace("“", '"')
        .replace("”", '"')
    )
    replacements = {
        "![Strict nonclosure certificate](FIN_Programs_91_100_Critical_Nonlocal_Operational_Figures/program91_strict_nonclosure_certificate.png)": "FINFIG91",
        "![Critical tuning](FIN_Programs_91_100_Critical_Nonlocal_Operational_Figures/program92_critical_tuning.png)": "FINFIG92",
        "![Graphon continuum](FIN_Programs_91_100_Critical_Nonlocal_Operational_Figures/program93_graphon_continuum.png)": "FINFIG93",
        "![Long-range propagation](FIN_Programs_91_100_Critical_Nonlocal_Operational_Figures/program94_long_range_propagation.png)": "FINFIG94",
        "![Process-tensor design](FIN_Programs_91_100_Critical_Nonlocal_Operational_Figures/program97_process_tensor_design.png)": "FINFIG97",
        "![Feedback ledger](FIN_Programs_91_100_Critical_Nonlocal_Operational_Figures/program98_feedback_ledger.png)": "FINFIG98",
        "![Damping completion](FIN_Programs_91_100_Critical_Nonlocal_Operational_Figures/program100_damping_completion_atom.png)": "FINFIG100",
    }
    for markdown, token in replacements.items():
        text = text.replace(markdown, f"\n\n{token}\n\n")
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
        ("[Proven, computer-assisted]", r"\statusComputer{}"),
        ("[Proven]", r"\statusProven{}"),
        ("[Strong evidence]", r"\statusStrong{}"),
        ("[Moderate evidence]", r"\statusModerate{}"),
        ("[Conditional]", r"\statusConditional{}"),
        ("[Refuted]", r"\statusRefuted{}"),
        ("[Open]", r"\statusOpen{}"),
    ]:
        body = body.replace(plain, macro)

    for filename in [
        "fin_programs_91_100_critical_nonlocal_operational.py",
        "test_fin_programs_91_100_critical_nonlocal_operational.py",
        "FIN_Programs_91_100_Critical_Nonlocal_Operational_Results.json",
        "FIN_Programs_91_100_External_Data_Intake_Template.json",
        "FIN_Programs_91_100_Critical_Nonlocal_Operational_Figures/",
    ]:
        body = body.replace(
            converter.escape_normal(filename), f"\\path{{{filename}}}"
        )

    TARGET.write_text(
        PREAMBLE + "\n" + body + "\n" + POSTAMBLE, encoding="utf-8"
    )
    print(f"Wrote {TARGET} ({TARGET.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
