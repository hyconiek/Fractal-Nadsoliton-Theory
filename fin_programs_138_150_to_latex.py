#!/usr/bin/env python3
"""Build archival LaTeX for FIN Programs 138--150."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Programs_138_150_State_Detector_Bridge_Monograph.md")
TARGET = Path("FIN_Programs_138_150_State_Detector_Bridge_Monograph.tex")

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
  pdftitle={FIN Programs 138--150: State Selection, Validated Fractional Dynamics, Detector Physics, and Bridge Architecture},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={KMS and entropy state-selection no-go theorems, validated fractional FFT, detector resolution, preparation resource theory, bridge diagram, and pre-data protocol},
  pdfkeywords={FIN, KMS states, maximum entropy, Morita equivalence, interval FFT, fractional dynamics, detector resolution, resource theory, operational physics}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs 138--150}
\fancyhead[R]{\small Release 10.14}
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
{\Large\bfseries FIN Research Monograph --- Release 10.14\par}
\vspace{11mm}
{\Huge\bfseries FIN Programs 138--150\par}
\vspace{6mm}
{\Large State Selection, Validated Fractional Dynamics,\\
Detector Physics, and Bridge Architecture\par}
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
\textbf{Scope.} Thirteen executed studies of state selection, KMS and
maximum-entropy obstructions, Morita stability, validated fractional
dynamics, detector resolution, joint system--instrument identification,
calibration-invariant observables, signed-preparation resources, completion
maps, and a frozen pre-data protocol.

\medskip
\textbf{Central result.} KMS theory, maximum entropy, and Morita equivalence
do not select the sector state producing \(9/5\). A formal interval FFT,
a reflection-asymmetry resource theorem, and a calibration-invariant
operational test are constructed.

\medskip
\textbf{Guardrail.} No canonical sector state, strict selector, internal
physical units, completed legacy--strict bridge, role transfer,
\(L_{\mathrm{total}}\), Theory-of-Everything closure, or external validation
is claimed.

\medskip
\textbf{License.} Creative Commons Attribution 4.0 International (CC BY 4.0).
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}

\pagenumbering{roman}
\chapter*{Publication statement}
\addcontentsline{toc}{chapter}{Publication statement}
This release reports Programs 138--150 and recommends Programs 151--163.
Analytic proofs, formal interval certificates, conditional constructions,
synthetic evidence, falsified candidates, and open physical obligations are
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
\path{FIN_Programs_138_150_State_Detector_Bridge_Results.json}.
The executable, forty-two tests, thirteen figures, and immutable pre-data
protocol reproduce the release.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Programs 138--150: State Selection,
Validated Fractional Dynamics, Detector Physics, and Bridge Architecture}
(FIN Research Monograph, Release 10.14; Version 1.0.0) [Preprint]. Zenodo.
\end{document}
"""

FIGURES = {
    "program138_modular_kms_family.png": r"The block-Gibbs modular family realizes \(9/5\) only after a gap \(\theta=\log 2\) is supplied.",
    "program139_entropy_reference_classification.png": "Maximum relative entropy returns the chosen reference measure.",
    "program140_morita_amplification.png": "Morita-equivalent block amplifications change normalized Hilbert central weights.",
    "program141_validated_interval_fft.png": "Formal interval-FFT remainder bounds on every frequency cell in the declared window.",
    "program142_diophantine_discrepancy.png": "The finite Denjoy--Koksma/Ostrowski discrepancy modulus and observed discrepancy.",
    "program143_weighted_wave_estimates.png": r"The Sobolev point-evaluation constant diverges as \(s\) approaches \(1/2\).",
    "program144_detector_resolution.png": "Resolution removal succeeds for a smooth preparation and fails for an ideal delta.",
    "program145_joint_identifiability.png": "Joint inversion of fractional exponent, detector error, and apparatus memory.",
    "program146_calibration_invariant_ratio.png": "The IQR log-slope cancels multiplicative calibration and estimates the stable scaling exponent.",
    "program147_preparation_resource.png": "Reflection-covariant operations monotonically destroy rather than create signed preparation.",
    "program148_coupled_variational_family.png": "The coupled state--damping action exposes one unsourced reference exponent.",
    "program149_completion_map_diagram.png": "The typed completion diagram retains absent phase/frequency and role-transfer arrows.",
    "program150_predata_protocol_power.png": "Synthetic power audit of the immutable pre-data operational protocol.",
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
            r"!\[[^\]]*\]\(FIN_Programs_138_150_State_Detector_Bridge_Figures/"
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
        path = "FIN_Programs_138_150_State_Detector_Bridge_Figures/" + filename
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
