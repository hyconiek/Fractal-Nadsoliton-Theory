#!/usr/bin/env python3
"""Build archival LaTeX for FIN Programs 165--177."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Programs_165_177_Axiom_Falsification_Measurement_Monograph.md")
TARGET = Path("FIN_Programs_165_177_Axiom_Falsification_Measurement_Monograph.tex")

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
  pdftitle={FIN Programs 165--177: Axiom Falsification, Operational Identifiability, and Measurement Instruments},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Composition audit, AFIS independence, detector inference, phase completion, and finite measurement instruments},
  pdfkeywords={FIN, spectral theorem, axiom falsification, tensor composition, instruments, calibration, double slit, Dirichlet forms}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs 165--177}
\fancyhead[R]{\small Release 10.16}
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
\vspace*{8mm}
{\Large\bfseries FIN Research Monograph --- Release 10.16\par}
\vspace{10mm}
{\Huge\bfseries FIN Programs 165--177\par}
\vspace{6mm}
{\Large Axiom Falsification, Operational Identifiability,\\
and Measurement Instruments\par}
\vspace{16mm}
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
\textbf{Scope.} Thirteen executed studies of cancellation-aware kernel
analysis, tensor-composition falsification, formal axiom independence,
algorithmic nonresonance, categorical state laws, resource completeness,
finite detector inference, calibration controls, phase completion,
legacy--strict provenance, and operational double-slit measurement.

\medskip
\textbf{Central result.} The one-copy relation
\(\beta_F=\alpha_{\rm geo}/4=\log2\) fails tensor intensivity and is not a
derived thermal law. The spectral generator remains the minimal common
dynamical object; state, preparation, instrument, calibration and sourced
completion remain independent layers.

\medskip
\textbf{Guardrail.} No QW-2191 discharge, internal physical units, completed
legacy--strict bridge, role transfer, \(L_{\rm total}\), Theory-of-Everything
closure, or external physical validation is claimed.

\medskip
\textbf{License.} Creative Commons Attribution 4.0 International (CC BY 4.0).
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}

\pagenumbering{roman}
\chapter*{Publication statement}
\addcontentsline{toc}{chapter}{Publication statement}
This release reports Programs 165--177 and recommends Programs 178--190.
Analytic theorems, exact finite audits, synthetic evidence, conditional
constructions, refuted candidates and open physical obligations are kept
separate.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""

POSTAMBLE = r"""
\clearpage
\chapter*{Archival note}
\addcontentsline{toc}{chapter}{Archival note}
Machine-readable results are stored in
\path{FIN_Programs_165_177_Axiom_Falsification_Measurement_Results.json}.
The frozen composite protocol is stored in
\path{FIN_Programs_165_177_Composite_Preregistration.json}.
The executable, sixty-four tests, Lean finite-model source, and thirteen
figures reproduce the declared release scope.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Programs 165--177: Axiom Falsification,
Operational Identifiability, and Measurement Instruments}
(FIN Research Monograph, Release 10.16; Version 1.0.0) [Preprint]. Zenodo.
\end{document}
"""

FIGURES = {
    "program165_cancellation_tail.png": "Cancellation-aware Fourier--Dirichlet bounds reduce the oscillatory correction tail by more than three orders of magnitude.",
    "program166_A_ME_composition.png": r"The cardinality rule for \(\beta_F\) fails tensor intensivity beyond one copy.",
    "program167_AFIS_formal_matrix.png": "Capability matrix for all sixty-four finite AFIS axiom subsets.",
    "program168_algorithmic_modulus.png": "Finite continued-fraction moduli remain positive but do not imply an all-scale polynomial rate.",
    "program169_monoidal_valuation.png": "Tensor multiplicativity leaves a family of valuations; direct-sum additivity selects Hilbert dimension.",
    "program170_A_ME_stability.png": r"The conditional value \(\eta=9/5\) is smooth but not exact under generic perturbation.",
    "program171_full_resource_counterexample.png": "Equal transverse reflection monotone does not imply equal full-state asymmetry.",
    "program172_DKW_detector.png": "Distribution-free iid coverage survives bounded pixels but fails under undeclared block memory.",
    "program173_calibration_control.png": "A shared-gain reference control removes exponent bias in the declared nuisance model.",
    "program174_composite_confusion.png": "Held-out confusion matrix for the frozen six-class synthetic protocol.",
    "program175_phase_cocycle.png": "The exact phase correction imports the strict-minus-legacy frequency and phase differences.",
    "program176_completion_stop_rule.png": "The typed completion audit stops at the first unsourced damping edge.",
    "program177_double_slit_instrument.png": "A finite POVM, dephasing channel, detector map and record separate wave dynamics from measurement.",
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
        .replace("∎", r"\(\square\)")
    )
    for filename in FIGURES:
        pattern = re.compile(
            r"!\[[^\]]*\]\(FIN_Programs_165_177_Axiom_Falsification_Measurement_Figures/"
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
        path = "FIN_Programs_165_177_Axiom_Falsification_Measurement_Figures/" + filename
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
