#!/usr/bin/env python3
"""Build archival LaTeX for FIN Programs 81--90."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Programs_81_90_Asymptotic_Operational_Completion_Monograph.md")
TARGET = Path("FIN_Programs_81_90_Asymptotic_Operational_Completion_Monograph.tex")

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
  pdftitle={FIN Programs 81--90: Asymptotic and Operational Completion},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Asymptotic spectral obstructions, quotient actions, feedback thermodynamics, and process tensors},
  pdfkeywords={FIN, Schur complement, harmonic aliasing, graph limits, process tensor, feedback thermodynamics}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs 81--90}
\fancyhead[R]{\small Release 10.9}
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
{\Large\bfseries FIN Research Monograph --- Release 10.9\par}
\vspace{11mm}
{\Huge\bfseries FIN Programs 81--90\par}
\vspace{6mm}
{\Large Asymptotic Spectral Obstructions\\
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
\textbf{Scope.} Ten executed studies of harmonic-alias Schur reduction,
large-\(N\) spectral obstructions, dense-kernel norm artefacts, fixed symbols,
approximate locality, quotient actions, robust calibration, chiral-state
dynamics, feedback thermodynamics, process tensors, and external-data
admission.

\medskip
\textbf{Central correction.} Frobenius convergence of dense normalized rows
does not imply uniform symbol convergence or a continuum fixed point.

\medskip
\textbf{Guardrail.} No strict selector, physical-unit source, Lorentz closure,
legacy role transfer, completed legacy--strict bridge, or external physical
validation is claimed.

\medskip
\textbf{License.} Creative Commons Attribution 4.0 International (CC BY 4.0).
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}

\pagenumbering{roman}
\chapter*{Publication statement}
\addcontentsline{toc}{chapter}{Publication statement}
This release reports Programs 81--90 and recommends Programs 91--102.
Theorem-level statements, numerical evidence, conditioned operational models,
failed inferences, and open obligations are labeled separately.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""

POSTAMBLE = r"""
\clearpage
\chapter*{Archival note}
\addcontentsline{toc}{chapter}{Archival note}
Machine-readable results are stored in
\path{FIN_Programs_81_90_Asymptotic_Operational_Completion_Results.json}.
The executable, fifteen tests, six figures, and external-data intake template
reproduce the release.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Programs 81--90: Asymptotic Spectral
Obstructions and Operational Completion Tests} (FIN Research Monograph,
Release 10.9; Version 1.0.0) [Preprint].
\end{document}
"""

FIGURES = {
    "FINFIG81": (
        "FIN_Programs_81_90_Asymptotic_Operational_Completion_Figures/program81_asymptotic_symbol_audit.png",
        "Frobenius and uniform-symbol naturality metrics through retained size 49152.",
    ),
    "FINFIG82": (
        "FIN_Programs_81_90_Asymptotic_Operational_Completion_Figures/program82_fixed_symbol_classification.png",
        "Normalized harmonic-alias flow toward the constant-symbol boundary for positive-mass families.",
    ),
    "FINFIG83": (
        "FIN_Programs_81_90_Asymptotic_Operational_Completion_Figures/program83_dense_kernel_universality.png",
        "Generic dense-row Frobenius dilution alongside persistent uniform spectral defects.",
    ),
    "FINFIG84": (
        "FIN_Programs_81_90_Asymptotic_Operational_Completion_Figures/program84_locality_bounds.png",
        "Finite-volume locality bounds combining a factorial series tail with Duhamel truncation error.",
    ),
    "FINFIG87": (
        "FIN_Programs_81_90_Asymptotic_Operational_Completion_Figures/program87_chiral_state_stability.png",
        "Paired reflection-symmetric chiral sectors and an explicitly biased state law.",
    ),
    "FINFIG88": (
        "FIN_Programs_81_90_Asymptotic_Operational_Completion_Figures/program88_feedback_thermodynamics.png",
        "Feedback work, free-energy change, and mutual information versus binary measurement error.",
    ),
}


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    text = (
        text.replace("–", "--")
        .replace("—", "---")
        .replace("’", "'")
        .replace("“", '"')
        .replace("”", '"')
    )
    replacements = {
        "![Asymptotic symbol audit](FIN_Programs_81_90_Asymptotic_Operational_Completion_Figures/program81_asymptotic_symbol_audit.png)": "FINFIG81",
        "![Fixed-symbol classification](FIN_Programs_81_90_Asymptotic_Operational_Completion_Figures/program82_fixed_symbol_classification.png)": "FINFIG82",
        "![Dense-kernel universality](FIN_Programs_81_90_Asymptotic_Operational_Completion_Figures/program83_dense_kernel_universality.png)": "FINFIG83",
        "![Approximate locality bounds](FIN_Programs_81_90_Asymptotic_Operational_Completion_Figures/program84_locality_bounds.png)": "FINFIG84",
        "![Chiral state stability](FIN_Programs_81_90_Asymptotic_Operational_Completion_Figures/program87_chiral_state_stability.png)": "FINFIG87",
        "![Feedback thermodynamics](FIN_Programs_81_90_Asymptotic_Operational_Completion_Figures/program88_feedback_thermodynamics.png)": "FINFIG88",
    }
    for markdown, token in replacements.items():
        text = text.replace(markdown, f"\n\n{token}\n\n")
    text = re.sub(r"```[a-zA-Z]*\n(.*?)```", r"\1", text, flags=re.DOTALL)
    text = re.sub(r"`([^`]+)`", r"\1", text)
    text = re.sub(
        r"^(#{1,4})\s+\d+(?:\.\d+)*\.\s+",
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
        ("[Proven]", r"\statusProven{}"),
        ("[Strong evidence]", r"\statusStrong{}"),
        ("[Moderate evidence]", r"\statusModerate{}"),
        ("[Conditional]", r"\statusConditional{}"),
        ("[Refuted]", r"\statusRefuted{}"),
        ("[Open]", r"\statusOpen{}"),
    ]:
        body = body.replace(plain, macro)

    for filename in [
        "fin_programs_81_90_asymptotic_operational_completion.py",
        "test_fin_programs_81_90_asymptotic_operational_completion.py",
        "FIN_Programs_81_90_Asymptotic_Operational_Completion_Results.json",
        "FIN_Programs_81_90_External_Data_Intake_Template.json",
        "FIN_Programs_81_90_Asymptotic_Operational_Completion_Figures/",
    ]:
        body = body.replace(
            converter.escape_normal(filename), f"\\path{{{filename}}}"
        )

    for digest in [
        "7c533e7433fc2cfb34b3ce6f25068e9cd2ceaa3b73333d97ca67a31a6ddd7b98"
    ]:
        digest_tex = (
            "\\texttt{"
            + "}\\allowbreak\\texttt{".join(
                digest[i : i + 16] for i in range(0, len(digest), 16)
            )
            + "}"
        )
        body = body.replace(digest, digest_tex)

    TARGET.write_text(PREAMBLE + "\n" + body + "\n" + POSTAMBLE, encoding="utf-8")
    print(f"Wrote {TARGET} ({TARGET.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
