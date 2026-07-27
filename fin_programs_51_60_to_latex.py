#!/usr/bin/env python3
"""Build archival LaTeX for FIN Programs 51--60."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Programs_51_60_Green_Fractal_Physics_Monograph.md")
TARGET = Path("FIN_Programs_51_60_Green_Fractal_Physics_Monograph.tex")

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
  pdftitle={FIN Programs 51--60: Green Response and Fractal Information Compression},
  pdfauthor={Krzysztof \.Zuchowski},
  pdfsubject={Green response, fractal information compression, and the minimal bridge to physics},
  pdfkeywords={FIN, Green function, fractal compression, information geometry, operational physics}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs 51--60}
\fancyhead[R]{\small Release 10.6}
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
{\Large\bfseries FIN Research Monograph --- Release 10.6\par}
\vspace{11mm}
{\Huge\bfseries FIN Programs 51--60\par}
\vspace{6mm}
{\Large Green Response, Fractal Information Compression,\\and the Minimal Bridge to Physics\par}
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
\textbf{Scope.} Ten executed programs testing Green response, corrected
information geometry, support transfer, exact fractal Schur compression,
spectral dimension, chiral response, operational records, and the dimensional
bridge to physics.

\medskip
\textbf{Ontology and guardrail.} The nadsoliton is treated as primordial
information in a solitonic state. No lower informational layer, silent
legacy--strict identity, role transfer, strict selector, physical unit, or
Theory-of-Everything closure is claimed.

\medskip
\textbf{Reproducibility.}
\path{fin_programs_51_60_green_fractal_physics.py},
\path{test_fin_programs_51_60_green_fractal_physics.py}, JSON, and figures.

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
\path{FIN_Programs_51_60_Green_Fractal_Physics_Results.json}. The executable
and ten regression tests reproduce the finite results and all four figures.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Programs 51--60: Green Response,
Fractal Information Compression, and the Minimal Bridge to Physics}
(FIN Research Monograph, Release 10.6; Version 1.0.0) [Preprint].

\end{document}
"""

FIGURES = {
    "FINFIGGREEN": (
        "FIN_Programs_51_60_Green_Fractal_Physics_Figures/program52_green_geometry.png",
        r"Squared Green embedding \(R\) and the guaranteed Hilbert metric \(\sqrt R\) for the strict \(C_{12}\) realization.",
    ),
    "FINFIGCOMP": (
        "FIN_Programs_51_60_Green_Fractal_Physics_Figures/program55_fractal_compression.png",
        r"Exact binary Green-Schur compression preserves retained response but does not generate the strict profile.",
    ),
    "FINFIGDIM": (
        "FIN_Programs_51_60_Green_Fractal_Physics_Figures/program56_spectral_dimension.png",
        r"Finite local spectral dimension depends on support, diffusion window, and positivity repair.",
    ),
    "FINFIGRECORD": (
        "FIN_Programs_51_60_Green_Fractal_Physics_Figures/program59_operational_records.png",
        r"Unitary coherent, diffusive, and incoherent-control records from one declared strict generator.",
    ),
}


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    replacements = {
        "![Green geometry](FIN_Programs_51_60_Green_Fractal_Physics_Figures/program52_green_geometry.png)": "FINFIGGREEN",
        "![Fractal compression](FIN_Programs_51_60_Green_Fractal_Physics_Figures/program55_fractal_compression.png)": "FINFIGCOMP",
        "![Spectral dimension](FIN_Programs_51_60_Green_Fractal_Physics_Figures/program56_spectral_dimension.png)": "FINFIGDIM",
        "![Operational records](FIN_Programs_51_60_Green_Fractal_Physics_Figures/program59_operational_records.png)": "FINFIGRECORD",
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
        "fin_programs_51_60_green_fractal_physics.py",
        "test_fin_programs_51_60_green_fractal_physics.py",
        "FIN_Programs_51_60_Green_Fractal_Physics_Results.json",
        "FIN_Programs_51_60_Green_Fractal_Physics_Figures/",
    ]:
        body = body.replace(converter.escape_normal(filename), f"\\path{{{filename}}}")

    TARGET.write_text(PREAMBLE + "\n" + body + "\n" + POSTAMBLE, encoding="utf-8")
    print(f"Wrote {TARGET} ({TARGET.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
