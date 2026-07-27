#!/usr/bin/env python3
"""Build archival LaTeX for FIN Programs 204--216."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Programs_204_216_Categorical_Catalytic_External_Monograph.md")
TARGET = Path("FIN_Programs_204_216_Categorical_Catalytic_External_Monograph.tex")

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
  pdftitle={FIN Programs 204--216: Categorical Measures, Perfect References, and External Falsification Gates},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Categorical state selection, exponent cocycles, symmetry-reference catalysis, finite-shot inference and external spectral falsification},
  pdfkeywords={FIN, spectral theory, Morita equivalence, operator algebra, catalytic asymmetry, Dirichlet form, conformal e-process, moment problem, operational physics}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs 204--216}
\fancyhead[R]{\small Release 10.19}
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
\vspace*{8mm}
{\Large\bfseries FIN Research Monograph --- Release 10.19\par}
\vspace{10mm}
{\Huge\bfseries FIN Programs 204--216\par}
\vspace{6mm}
{\Large Categorical Measures, Perfect References,\\
and External Falsification Gates\par}
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
\textbf{Scope.} Thirteen executed studies of categorical central measures,
exponent cocycles, formal certification, catalytic symmetry references,
hidden mixing, sequential open-set detection, finite-shot process
tomography, environmental moment bounds, phase naturality, scale-free
spectral falsification, external intake, and prediction gating.

\medskip
\textbf{Central result.} Morita-permutation invariance conditionally selects
uniform central weights, while naturality under all unital homomorphisms is
impossible. A perfect asymmetric \(\mathbb Z_2\) reference enables exact
catalytic conversion but supplies, rather than derives, the missing
orientation bit.

\medskip
\textbf{Guardrail.} No QW-2191 discharge, strict selector, canonical physical
unit, target-independent strict exponent source, completed legacy--strict
bridge, role transfer, \(L_{\rm total}\), external physical validation, or
Theory-of-Everything closure is claimed.

\medskip
\textbf{License.} Creative Commons Attribution 4.0 International (CC BY 4.0).
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}

\pagenumbering{roman}
\chapter*{Publication statement}
\addcontentsline{toc}{chapter}{Publication statement}
This release reports Programs 204--216 and recommends Programs 217--229.
Analytic proofs, exact finite computation, machine-checked proof fragments,
synthetic evidence, conditional constructions, unavailable infrastructure,
and failed external-data gates are kept separate.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""

POSTAMBLE = r"""
\clearpage
\chapter*{Release reproducibility}
\addcontentsline{toc}{chapter}{Release reproducibility}
\begin{Verbatim}[fontsize=\small]
python3 fin_programs_204_216_categorical_catalytic_external.py
python3 -m unittest -v \
  test_fin_programs_204_216_categorical_catalytic_external.py
\end{Verbatim}

Machine-readable results are stored in
\path{FIN_Programs_204_216_Categorical_Catalytic_External_Results.json}.
The release contains six frozen auxiliary contracts or certificates, two Lean
sources, 77 tests, and thirteen figures. The dependency-free Lean core probe
compiles locally; the Mathlib-dependent general graph library does not.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Programs 204--216: Categorical Measures,
Perfect References, and External Falsification Gates}
(FIN Research Monograph, Release 10.19; Version 1.0.0) [Preprint]. Zenodo.
\end{document}
"""

FIGURES = {
    "program204_morita_central_measure.png": "Morita-permutation invariance conditionally selects uniform central weights; the represented-space trace remains a different choice.",
    "program205_eta_cocycle.png": r"The additive exponent cocycle reaches \(9/5\) along infinitely many rate--time pairs and therefore does not select the endpoint.",
    "program206_arb_environment.png": "The reproducibility contract is complete, but the local Arb engine and Docker server are unavailable.",
    "program207_lean_build.png": "The dependency-free Lean core compiles and the exact witness pack passes; the general Mathlib library remains uncompiled.",
    "program208_catalytic_conversion.png": r"Imperfect reference states fail the sampled necessary conversion test, whereas the perfect \(\mathbb Z_2\) reference reaches the exact boundary.",
    "program209_hidden_mixing.png": "Calibrated thinning restores valid coverage under hidden geometric mixing at a substantial information cost.",
    "program210_conformal_eprocess.png": "A fresh-block conformal e-process maintains time-uniform null control and rapidly rejects the declared temporal alternatives.",
    "program211_finite_shot_tomography.png": "Finite-shot simultaneous regions cover the Gaussian process parameters and detect the declared memory term.",
    "program212_environment_moments.png": r"Toeplitz positivity gives the sharp environmental moment interval \(7/25\le c_2\le1\).",
    "program213_phase_no_go.png": "Natural endomorphisms keep the legacy oscillatory datum in its finite torsion orbit and miss the strict non-torsion target.",
    "program214_scale_free_protocol.png": "The frozen projective fingerprint accepts exact rescalings and small perturbation, while rejecting large and structurally incompatible alternatives.",
    "program215_external_bundle.png": "No local candidate satisfies the independent eleven-field external intake gate.",
    "program216_external_prediction_lock.png": "The external prediction remains locked because the prerequisite independent bundle is absent.",
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
            r"!\[[^\]]*\]\(FIN_Programs_204_216_Categorical_Catalytic_External_Figures/"
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
    abstract_at = text.index("# Abstract")
    text = (
        "## Confidence convention\n\n"
        "The report distinguishes analytic proof, exact finite computation, "
        "machine-checked fragments, synthetic evidence, conditional "
        "construction, and open obligations. No unavailable certificate or "
        "external result is promoted by prose.\n\n"
        + text[abstract_at:]
    )
    body = converter.body_from_markdown(text)
    body = body.replace("Ż", r"\.Z").replace("ż", r"\.z")
    for filename, caption in FIGURES.items():
        token = f"FINFIG\\_{filename.replace('_', r'\_')}"
        path = "FIN_Programs_204_216_Categorical_Catalytic_External_Figures/" + filename
        body = body.replace(
            token,
            "\\begin{figure}[htbp]\n"
            "\\centering\n"
            f"\\includegraphics[width=0.94\\textwidth]{{{path}}}\n"
            f"\\caption{{{caption}}}\n"
            "\\end{figure}",
        )
    TARGET.write_text(PREAMBLE + "\n" + body + "\n" + POSTAMBLE, encoding="utf-8")
    print(TARGET)


if __name__ == "__main__":
    main()
