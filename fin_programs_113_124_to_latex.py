#!/usr/bin/env python3
"""Build archival LaTeX for FIN Programs 113--124."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Programs_113_124_Constructive_Completion_Monograph.md")
TARGET = Path("FIN_Programs_113_124_Constructive_Completion_Monograph.tex")

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
  pdftitle={FIN Programs 113--124: Constructive Completion Objects and Fractional Operational Physics},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Global Schur theorem, stable process, operational completion, Z12 fibre source, damping cocycle},
  pdfkeywords={FIN, Schur complement, stable process, fractional operator, operational physics, homology, characters, damping cocycle}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs 113--124}
\fancyhead[R]{\small Release 10.12}
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
{\Large\bfseries FIN Research Monograph --- Release 10.12\par}
\vspace{11mm}
{\Huge\bfseries FIN Programs 113--124\par}
\vspace{6mm}
{\Large Constructive Completion Objects\\
and Fractional Operational Physics\par}
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
\textbf{Scope.} Twelve constructive studies of global strict-kernel
nonclosure, stable fractional limits, operational process completion,
adaptive and variational source obstructions, correlated apparatus memory,
a new \(Z_{12}\) homological--character fibre object, and a conditional
damping cocycle.

\medskip
\textbf{Central construction.}
\(\mathcal F_p=\widetilde H_0(\ker m_p)\oplus\mathcal X_p^-\) has dimensions
\((1,2,2,2,2)\). Its uniform trace is \(9/5\), but the strict theory does not
yet force that trace. Conditional on it and a nonzero multiplicative tail,
\((\eta,\beta)=(9/5,1)\) follows.

\medskip
\textbf{Guardrail.} No strict selector, physical-unit source, completed
legacy--strict bridge, role transfer, role-bearing \(L_{\mathrm{total}}\), or
external physical validation is claimed.

\medskip
\textbf{License.} Creative Commons Attribution 4.0 International (CC BY 4.0).
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}

\pagenumbering{roman}
\chapter*{Publication statement}
\addcontentsline{toc}{chapter}{Publication statement}
This release reports Programs 113--124 and recommends Programs 125--137.
Proofs, computer-assisted certificates, conditional constructions, finite
exhaustions, synthetic operational tests, and open source obligations are
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
\path{FIN_Programs_113_124_Constructive_Completion_Results.json}.
The executable, twenty-eight tests, eleven figures, and frozen external-data
intake reproduce the release.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Programs 113--124: Constructive Completion
Objects and Fractional Operational Physics}
(FIN Research Monograph, Release 10.12; Version 1.0.0) [Preprint].
\end{document}
"""

FIGURES = {
    "FINFIG113": (
        "FIN_Programs_113_124_Constructive_Completion_Figures/program113_global_mass_theorem.png",
        "The exact positive native--Schur ratio gap for every positive mass.",
    ),
    "FINFIG114": (
        "FIN_Programs_113_124_Constructive_Completion_Figures/program114_continuous_parameter_box.png",
        "Positive analytic coefficient and normalization bounds on the continuous strict-parameter box.",
    ),
    "FINFIG115": (
        "FIN_Programs_113_124_Constructive_Completion_Figures/program115_effective_abelian_certificate.png",
        "Certified finite-wave-number remainder interval around the fractional Abelian law.",
    ),
    "FINFIG116": (
        "FIN_Programs_113_124_Constructive_Completion_Figures/program116_stable_invariance_principle.png",
        r"Finite characteristic-function check of convergence toward the symmetric \(4/5\)-stable law.",
    ),
    "FINFIG117": (
        "FIN_Programs_113_124_Constructive_Completion_Figures/program117_fractional_operational_process.png",
        "Dimensionless wave/diffusion record separation for the constructed fractional operational process.",
    ),
    "FINFIG118": (
        "FIN_Programs_113_124_Constructive_Completion_Figures/program118_local_fractional_crossover.png",
        "Scale-dependent local and fractional contributions to the crossover symbol.",
    ),
    "FINFIG119": (
        "FIN_Programs_113_124_Constructive_Completion_Figures/program119_inverse_fisher_potential.png",
        "The unique inverse strict Fisher potential and the target-independent envelope potential.",
    ),
    "FINFIG120": (
        "FIN_Programs_113_124_Constructive_Completion_Figures/program120_variational_grammar.png",
        "Exhaustive ranking of 440 target-free variational grammar candidates.",
    ),
    "FINFIG121": (
        "FIN_Programs_113_124_Constructive_Completion_Figures/program121_hidden_markov_apparatus.png",
        "Synthetic held-out inference with calibrated hidden-Markov apparatus memory.",
    ),
    "FINFIG122": (
        "FIN_Programs_113_124_Constructive_Completion_Figures/program122_homological_character_functor.png",
        r"Homology and character contributions to the exact fibre dimension vector \((1,2,2,2,2)\).",
    ),
    "FINFIG123": (
        "FIN_Programs_113_124_Constructive_Completion_Figures/program123_conditional_damping_cocycle.png",
        "Exponent and dyadic retention as functions of the remaining trace weight.",
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
        "![Global mass theorem](FIN_Programs_113_124_Constructive_Completion_Figures/program113_global_mass_theorem.png)": "FINFIG113",
        "![Continuous box](FIN_Programs_113_124_Constructive_Completion_Figures/program114_continuous_parameter_box.png)": "FINFIG114",
        "![Abelian certificate](FIN_Programs_113_124_Constructive_Completion_Figures/program115_effective_abelian_certificate.png)": "FINFIG115",
        "![Stable limit](FIN_Programs_113_124_Constructive_Completion_Figures/program116_stable_invariance_principle.png)": "FINFIG116",
        "![Fractional process](FIN_Programs_113_124_Constructive_Completion_Figures/program117_fractional_operational_process.png)": "FINFIG117",
        "![Crossover operator](FIN_Programs_113_124_Constructive_Completion_Figures/program118_local_fractional_crossover.png)": "FINFIG118",
        "![Inverse potential](FIN_Programs_113_124_Constructive_Completion_Figures/program119_inverse_fisher_potential.png)": "FINFIG119",
        "![Variational grammar](FIN_Programs_113_124_Constructive_Completion_Figures/program120_variational_grammar.png)": "FINFIG120",
        "![HMM apparatus](FIN_Programs_113_124_Constructive_Completion_Figures/program121_hidden_markov_apparatus.png)": "FINFIG121",
        "![Fibre functor](FIN_Programs_113_124_Constructive_Completion_Figures/program122_homological_character_functor.png)": "FINFIG122",
        "![Damping cocycle](FIN_Programs_113_124_Constructive_Completion_Figures/program123_conditional_damping_cocycle.png)": "FINFIG123",
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
