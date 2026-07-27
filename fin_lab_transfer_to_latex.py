#!/usr/bin/env python3
"""Build archival LaTeX for the FIN P240--P242 laboratory transfer package."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Laboratory_Transfer_Package_P240_242.md")
TARGET = Path("FIN_Laboratory_Transfer_Package_P240_242.tex")

PREAMBLE = r"""\documentclass[11pt,a4paper,openany]{report}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
\usepackage{microtype}
\usepackage{amsmath,amssymb,amsthm,mathtools,bm}
\usepackage{booktabs,array,longtable,tabularx}
\usepackage{xcolor}
\usepackage[margin=23mm,headheight=23pt]{geometry}
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
\definecolor{fingray}{HTML}{4D4D4D}

\hypersetup{
  unicode=true,
  colorlinks=true,
  linkcolor=finblue,
  citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Laboratory Transfer Package: P240 Optimal Tomography, P241 Validator, P242 One-Shot Pipeline},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Scientific transfer package and executable specification for a blind twelve-state FIN laboratory experiment},
  pdfkeywords={FIN, heat-kernel tomography, continuous-time quantum walk, Markov semigroup, photonics, double slit, blind custody, preregistration, spectral fingerprint}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Laboratory Transfer Package}
\fancyhead[R]{\small Release 10.21}
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

\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{7mm}
{\Large\bfseries FIN Laboratory Transfer Package --- Release 10.21\par}
\vspace{10mm}
{\Huge\bfseries P240 Optimal Tomography\par}
\vspace{3mm}
{\Huge\bfseries + P241 Blind-Custody Validator\par}
\vspace{3mm}
{\Huge\bfseries + P242 One-Shot Pipeline\par}
\vspace{8mm}
{\Large A scientific transfer package and executable specification\\
for theoretical and experimental physicists\par}
\vspace{15mm}
{\Large Krzysztof \.Zuchowski\par}
\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
\vspace{2mm}
{\normalsize ORCID: \href{https://orcid.org/0009-0002-0909-3613}{0009-0002-0909-3613}\par}
\vspace{8mm}
{\large 27 July 2026\par}
{\normalsize Publication --- Preprint and executable specification; version 1.0.0\par}
\vfill
\begin{minipage}{0.91\textwidth}
\small
\textbf{Transfer target.} One finite, calibrated, twelve-state laboratory
test of the strict FIN generator, with an independently signed event bundle,
a sealed \(2\tau\) holdout and exactly one registered analysis.

\medskip
\textbf{Central mathematical object.}
\[
W=K_s,\qquad s=1.660307278766099,\qquad
A=sI-W\geq0,\qquad U_t=e^{-itA},\quad P_t=e^{-tA}.
\]

\textbf{Critical boundary.} This package supplies a design, validator and
analysis pipeline. It supplies no external events, physical clock, apparatus,
strict selector, legacy--strict role transfer, \(L_{\rm total}\), Standard
Model/gravity derivation or Theory-of-Everything closure.

\medskip
\textbf{Repository.}
\url{https://github.com/hyconiek/Fractal-Nadsoliton-Theory}

\textbf{Related project DOI supplied by the author.}
\url{https://doi.org/10.5281/zenodo.21435332}

\textbf{License.} Creative Commons Attribution 4.0 International (CC BY 4.0).
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}

\pagenumbering{roman}
\chapter*{Transfer statement}
\addcontentsline{toc}{chapter}{Transfer statement}
This document is intended to be handed, together with the hashed executable
files, to a laboratory team.  It summarizes the finite spectral theory,
the wave--diffusion duality, the missing operational object, the three
candidate apparatus classes, the registered acquisition design, the
independent custody contract and the one-shot analysis rule.

Programs 240--242 in this document are laboratory programme numbers from the
Release-10.20 roadmap.  They are not the older internal FAR packets that
happen to contain the same integers.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""

POSTAMBLE = r"""
\clearpage
\chapter*{Release reproducibility}
\addcontentsline{toc}{chapter}{Release reproducibility}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/matplotlib-finlab \
python3 fin_lab_p240_optimal_tomography.py --replicates 300

python3 fin_lab_p242_pipeline.py --write-lock-only

MPLCONFIGDIR=/tmp/matplotlib-finlab \
python3 -m unittest -v test_fin_lab_p240_242.py
\end{Verbatim}

Expected test result: 18 tests, all passing.  No P242 production execution is
part of this release because no signed external P241 bundle is supplied.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Laboratory Transfer Package: P240 Optimal
Tomography, P241 Blind-Custody Validator, and P242 One-Shot Pipeline}
(Release 10.21; Version 1.0.0) [Preprint and executable specification].
Independent Researcher --- Fractal Information Theory Project.
\end{document}
"""

FIGURES = {
    "duality_operational_map.png": (
        "The same positive generator supports coherent and diffusive "
        "functional calculi, while empirical probabilities require an "
        "operational completion."
    ),
    "platform_obligation_matrix.png": (
        "Fit of the three feasibility hypotheses to the five measurement "
        "duties. Independent custody remains a human and institutional gate."
    ),
    "platform_A_apparatus.png": (
        "Recommended direct P-A architecture. Coherent propagation and heat "
        "flow require separately verified physical realizations of the common "
        "target generator."
    ),
    "platform_C_double_slit.png": (
        "P-C event-level double-slit apparatus. Raw detections, four shutter "
        "configurations and calibrations are retained; a rendered image is "
        "not an admissible substitute."
    ),
    "custody_pipeline.png": (
        "Provider, registrar and analyst remain distinct. The sealed holdout "
        "is released only after hashes and analysis code are frozen."
    ),
    "p240_optimal_tomography.png": (
        "Program 240: the exact fastest-mode time optimum and deterministic "
        "finite-shot planning distributions."
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
        .replace("∎", r"\(\square\)")
    )
    for filename in FIGURES:
        pattern = re.compile(
            r"!\[[^\]]*\]\(FIN_Lab_P240_242_Transfer_Figures/"
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
    start = text.index("# Abstract")
    conversion_source = (
        "## Confidence convention\n\n"
        "Analytic proof, machine checking, synthetic planning, feasibility "
        "hypotheses, open obligations and scoped refutations are kept "
        "separate throughout this transfer package.\n\n"
        + text[start:]
    )
    body = converter.body_from_markdown(conversion_source)
    body = body.replace("Ż", r"\.Z").replace("ż", r"\.z")
    for filename, caption in FIGURES.items():
        token = f"FINFIG\\_{filename.replace('_', r'\_')}"
        path = f"FIN_Lab_P240_242_Transfer_Figures/{filename}"
        body = body.replace(
            token,
            "\\begin{figure}[htbp]\n"
            "\\centering\n"
            f"\\includegraphics[width=0.96\\textwidth]{{{path}}}\n"
            f"\\caption{{{caption}}}\n"
            "\\end{figure}",
        )
    TARGET.write_text(PREAMBLE + "\n" + body + "\n" + POSTAMBLE, encoding="utf-8")
    print(TARGET)


if __name__ == "__main__":
    main()
