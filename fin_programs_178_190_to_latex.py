#!/usr/bin/env python3
"""Build archival LaTeX for FIN Programs 178--190."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Programs_178_190_Composition_Process_Scale_Monograph.md")
TARGET = Path("FIN_Programs_178_190_Composition_Process_Scale_Monograph.tex")

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
  pdftitle={FIN Programs 178--190: Composition Laws, Process Memory, and Scale Obstructions},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Tensor state laws, compression torsors, covariant resources, process memory, environment nonuniqueness, phase and scale obstructions},
  pdfkeywords={FIN, spectral theorem, tensor composition, Dirichlet forms, resource theory, process tensor, calibration, scale torsor, phase cohomology}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs 178--190}
\fancyhead[R]{\small Release 10.17}
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
{\Large\bfseries FIN Research Monograph --- Release 10.17\par}
\vspace{10mm}
{\Huge\bfseries FIN Programs 178--190\par}
\vspace{6mm}
{\Large Composition Laws, Process Memory,\\
and Scale Obstructions\par}
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
\textbf{Scope.} Thirteen executed studies of tensor-intensive state laws,
strict compression scale, formal fractional certification, finite Dirichlet
dynamics, reflection-resource conversion, dependence-robust inference,
multi-control calibration, open-set falsification, process memory,
environment nonuniqueness, phase provenance, intrinsic scale and external
data intake.

\medskip
\textbf{Central result.} State, compression, environment and clock data live
on nontrivial moduli or torsors.  The spectral generator fixes none of their
representatives.  Complete finite resource and process objects can be
constructed conditionally without promoting them to strict FIN physics.

\medskip
\textbf{Guardrail.} No QW-2191 discharge, strict unit or beta source,
completed legacy--strict bridge, role transfer, \(L_{\rm total}\),
Theory-of-Everything closure, or external physical validation is claimed.

\medskip
\textbf{License.} Creative Commons Attribution 4.0 International (CC BY 4.0).
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}

\pagenumbering{roman}
\chapter*{Publication statement}
\addcontentsline{toc}{chapter}{Publication statement}
This release reports Programs 178--190 and recommends Programs 191--203.
Theorems, exact finite audits, synthetic evidence, conditional constructions,
refuted candidates, unavailable certifications, and open physical obligations
are kept separate.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""

POSTAMBLE = r"""
\clearpage
\chapter*{Archival note}
\addcontentsline{toc}{chapter}{Archival note}
Machine-readable results are stored in
\path{FIN_Programs_178_190_Composition_Process_Scale_Results.json}.
The frozen open-set challenge is stored in
\path{FIN_Programs_178_190_Open_Set_Preregistration.json}.
The executable, seventy-two tests, finite Lean source, and thirteen figures
reproduce the declared release scope.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Programs 178--190: Composition Laws,
Process Memory, and Scale Obstructions}
(FIN Research Monograph, Release 10.17; Version 1.0.0) [Preprint]. Zenodo.
\end{document}
"""

FIGURES = {
    "program178_tensor_intensive_classification.png": "Replication intensivity leaves multiple state laws; label splitting removes nonconstant ratios without selecting a constant.",
    "program179_compression_scale_torsor.png": r"All positive \(\beta\) values collapse to one universal profile after coordinate rescaling.",
    "program180_ball_certificate_ledger.png": "Only part of the desired fractional certificate currently resides in one directed-rounding chain.",
    "program181_dirichlet_core.png": "Unitarity, heat row sums and heat positivity for the exact four-cycle Dirichlet generator.",
    "program182_reflection_convertibility.png": "Complete Choi-feasibility region from the source Bloch state with transverse coordinate 0.6.",
    "program183_block_DKW.png": "Correct dependence coverage requires the number of independent blocks and produces much wider intervals.",
    "program184_multi_control.png": "Known controls remove shared nonlinear detector gain from the target exponent.",
    "program185_open_set.png": "Closed-set distance thresholds reject several unseen laws but cannot reject a feature-invariant ordering challenge.",
    "program186_process_tensor.png": "A three-component intervention vector separates four declared visibility-loss mechanisms.",
    "program187_environment_nonuniqueness.png": "One fixed system generator admits a continuum of environment couplings and dephasing visibilities.",
    "program188_phase_cohomology.png": "Legacy cyclotomic, strict infinite-order and correction phase values on the unit circle.",
    "program189_scale_orbit.png": "Raw spectral scale changes while normalized spectral predictions remain invariant.",
    "program190_external_intake.png": "None of three local external-lineage candidates passes the eleven-field operational intake gate.",
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
            r"!\[[^\]]*\]\(FIN_Programs_178_190_Composition_Process_Scale_Figures/"
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
        path = "FIN_Programs_178_190_Composition_Process_Scale_Figures/" + filename
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
