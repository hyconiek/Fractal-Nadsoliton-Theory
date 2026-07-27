#!/usr/bin/env python3
"""Build archival LaTeX for FIN Programs 101--112."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Programs_101_112_Fractional_Source_Completion_Monograph.md")
TARGET = Path("FIN_Programs_101_112_Fractional_Source_Completion_Monograph.tex")

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
  pdftitle={FIN Programs 101--112: Fractional Limits, Source Tests, and Operational Completion},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Schur nonclosure, graphon and fractional limits, adaptive dynamics, operational models, exponent source},
  pdfkeywords={FIN, fractional operator, Schur complement, graphon, stable process, information geometry, process tensor}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs 101--112}
\fancyhead[R]{\small Release 10.11}
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
{\Large\bfseries FIN Research Monograph --- Release 10.11\par}
\vspace{11mm}
{\Huge\bfseries FIN Programs 101--112\par}
\vspace{6mm}
{\Large Fractional Limits, Source Tests,\\
and Operational Completion\par}
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
\textbf{Scope.} Twelve executed studies of independent strict-kernel
certification, graphon/local/fractional continuum limits, long-range
propagation, adaptive information geometry, inverse variational sources,
selector boundaries, preregistration, apparatus memory, and the conditional
source skeleton for the exponent \(9/5\).

\medskip
\textbf{Central result.} The infinite strict lattice is a symmetric
fractional generator of order \(4/5\), while a local Laplacian requires an
additional singular localization and clock scaling.  The identity
\(4\ln2\mapsto2^{-4/5}\mapsto d^{-9/5}\) is exact but remains conditional on
an unsourced fivefold quotient and damping coupling.

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
This release reports Programs 101--112 and recommends Programs 113--124.
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
\path{FIN_Programs_101_112_Fractional_Source_Completion_Results.json}.
The executable, twenty-two tests, nine figures, and frozen process-tensor
preregistration reproduce the release.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Programs 101--112: Fractional Limits,
Source Tests, and Operational Completion}
(FIN Research Monograph, Release 10.11; Version 1.0.0) [Preprint].
\end{document}
"""

FIGURES = {
    "FINFIG101": (
        "FIN_Programs_101_112_Fractional_Source_Completion_Figures/program101_independent_interval_certificate.png",
        "Independent directed-interval separation of the native and Schur scale-invariant spectral ratios.",
    ),
    "FINFIG102": (
        "FIN_Programs_101_112_Fractional_Source_Completion_Figures/program102_nonclosure_region.png",
        "Certified ratio separation throughout the declared continuous mass interval.",
    ),
    "FINFIG103": (
        "FIN_Programs_101_112_Fractional_Source_Completion_Figures/program103_graphon_error_theorem.png",
        "Analytic operator-norm bound and observed low-mode convergence to the bounded graphon operator.",
    ),
    "FINFIG104": (
        "FIN_Programs_101_112_Fractional_Source_Completion_Figures/program104_singular_localizing_limit.png",
        "Modewise convergence of the singularly localized profile to the local Laplacian.",
    ),
    "FINFIG105": (
        "FIN_Programs_101_112_Fractional_Source_Completion_Figures/program105_fractional_tail_universality.png",
        r"Small-wave-number convergence to the fractional law \(C|q|^{4/5}\).",
    ),
    "FINFIG106": (
        "FIN_Programs_101_112_Fractional_Source_Completion_Figures/program106_long_range_semigroup.png",
        "One-big-jump far-tail ratio for the strict heat semigroup.",
    ),
    "FINFIG107": (
        "FIN_Programs_101_112_Fractional_Source_Completion_Figures/program107_adaptive_information_manifold.png",
        "Simplex-preserving Fisher gradient flow and monotone information functional.",
    ),
    "FINFIG111": (
        "FIN_Programs_101_112_Fractional_Source_Completion_Figures/program111_correlated_apparatus_ledger.png",
        "Apparatus-inclusive information-work excess with correlated detector memory.",
    ),
    "FINFIG112": (
        "FIN_Programs_101_112_Fractional_Source_Completion_Figures/program112_eta_source_skeleton.png",
        r"Discharged identity and open source obligations in the conditional \(\eta=9/5\) skeleton.",
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
    replacements = {
        "![Independent interval certificate](FIN_Programs_101_112_Fractional_Source_Completion_Figures/program101_independent_interval_certificate.png)": "FINFIG101",
        "![Nonclosure region](FIN_Programs_101_112_Fractional_Source_Completion_Figures/program102_nonclosure_region.png)": "FINFIG102",
        "![Graphon theorem](FIN_Programs_101_112_Fractional_Source_Completion_Figures/program103_graphon_error_theorem.png)": "FINFIG103",
        "![Localizing limit](FIN_Programs_101_112_Fractional_Source_Completion_Figures/program104_singular_localizing_limit.png)": "FINFIG104",
        "![Fractional universality](FIN_Programs_101_112_Fractional_Source_Completion_Figures/program105_fractional_tail_universality.png)": "FINFIG105",
        "![Long-range semigroup](FIN_Programs_101_112_Fractional_Source_Completion_Figures/program106_long_range_semigroup.png)": "FINFIG106",
        "![Adaptive manifold](FIN_Programs_101_112_Fractional_Source_Completion_Figures/program107_adaptive_information_manifold.png)": "FINFIG107",
        "![Correlated apparatus ledger](FIN_Programs_101_112_Fractional_Source_Completion_Figures/program111_correlated_apparatus_ledger.png)": "FINFIG111",
        "![Eta source skeleton](FIN_Programs_101_112_Fractional_Source_Completion_Figures/program112_eta_source_skeleton.png)": "FINFIG112",
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
            f"\\includegraphics[width=0.94\\textwidth]{{{path}}}\n"
            f"\\caption{{{caption}}}\n"
            "\\end{figure}"
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
    TARGET.write_text(PREAMBLE + "\n" + body + "\n" + POSTAMBLE, encoding="utf-8")
    print(TARGET)


if __name__ == "__main__":
    main()
