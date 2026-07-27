#!/usr/bin/env python3
"""Build archival LaTeX for FIN Programs 151--164."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Programs_151_164_Axiomatic_Operational_Foundations_Monograph.md")
TARGET = Path("FIN_Programs_151_164_Axiomatic_Operational_Foundations_Monograph.tex")

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
  pdftitle={FIN Programs 151--164: Axiomatic Operational Foundations, State Selection, and Falsifiable Measurement},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Axiomatic FIN foundations, modular equilibrium, operational identifiability, detector statistics, phase obstruction, and role-transfer guardrails},
  pdfkeywords={FIN, axioms, Dirichlet forms, spectral theorem, modular equilibrium, resource theory, instruments, calibration, fractional dynamics}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs 151--164}
\fancyhead[R]{\small Release 10.15}
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
{\Large\bfseries FIN Research Monograph --- Release 10.15\par}
\vspace{10mm}
{\Huge\bfseries FIN Programs 151--164\par}
\vspace{6mm}
{\Large Axiomatic Operational Foundations,\\
State Selection, and Falsifiable Measurement\par}
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
\textbf{Scope.} Fourteen executed studies of validated fractional dynamics,
Diophantine control, functorial states, modular equilibrium, signed
preparation resources, detector statistics, semiparametric identifiability,
phase completion, legacy-role obstruction, and minimal axioms.

\medskip
\textbf{Central result.} A conservative spectral Dirichlet generator is the
minimal common mathematical object. State selection, signed preparation,
instrument, calibration, and legacy--strict completion remain independent
obligations. One explicit axiom \(A_{\rm ME}\) selects \(\eta=9/5\), but is
not strict-derived.

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
This release reports Programs 151--164 and recommends Programs 165--177.
Analytic theorems, computer-assisted certificates, conditional axioms,
synthetic evidence, refuted candidates, and open physical obligations are
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
\path{FIN_Programs_151_164_Axiomatic_Operational_Results.json}.
The extracted six-axiom system is stored in
\path{FIN_Programs_151_164_Minimal_Axiomatic_System.json}.
The executable, sixty-two tests, and fourteen figures reproduce the release.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Programs 151--164: Axiomatic Operational
Foundations, State Selection, and Falsifiable Measurement}
(FIN Research Monograph, Release 10.15; Version 1.0.0) [Preprint]. Zenodo.
\end{document}
"""

FIGURES = {
    "program151_tighter_interval_fft.png": "The doubled interval FFT improves the formal enclosure without yet reaching the three-percent target.",
    "program152_diophantine_axiom.png": r"Finite candidates for \(q^\nu\lVert q\theta\rVert\) do not certify an all-scale Diophantine axiom.",
    "program153_groupoid_measure_simplex.png": "Groupoid invariance leaves a two-dimensional simplex of orbit masses.",
    "program154_axiomatic_modular_equilibrium.png": r"The axiom \(\beta_F=\log2\) selects the uniform central state and \(\eta=9/5\).",
    "program155_reflection_resource_complete.png": r"Every reflection-covariant orbit-line channel contracts the complete monotone \(\lvert r\rvert\).",
    "program156_detector_bias.png": "Absolute Gaussian resolution biases the finite-time IQR spreading slope downward.",
    "program157_semiparametric_rank.png": "A linear log-time apparatus response exactly confounds the fractional exponent.",
    "program158_iqr_finite_sample.png": "Empirical IQR-slope uncertainty agrees with the asymptotic quantile theorem.",
    "program159_adversarial_protocol.png": "The frozen binary rule detects scaling but is not specific to FIN.",
    "program160_phase_obstruction.png": "The period-eight legacy character and infinite-order strict phase are representation-theoretically distinct.",
    "program161_energy_grammar.png": r"The energy \(d-1\) is the unique shortest target-realizing formula in the declared grammar.",
    "program162_role_obstruction_matrix.png": "Necessary completion edges for nine legacy roles; none currently licenses transfer.",
    "program163_external_readiness.png": "The intake schema is executable, while no external dataset is admitted.",
    "program164_axiom_capability_matrix.png": "Dependency of each mathematical or physical capability on the six AFIS axioms.",
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
            r"!\[[^\]]*\]\(FIN_Programs_151_164_Axiomatic_Operational_Figures/"
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
        path = "FIN_Programs_151_164_Axiomatic_Operational_Figures/" + filename
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
