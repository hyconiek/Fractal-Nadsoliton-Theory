#!/usr/bin/env python3
"""Build archival LaTeX for FIN Programs 217--229."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Programs_217_229_Operational_Section_Spectral_Tomography_Monograph.md")
TARGET = Path("FIN_Programs_217_229_Operational_Section_Spectral_Tomography_Monograph.tex")

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
  pdftitle={FIN Programs 217--229: Operational Sections, Reference Costs, and Heat-Kernel Spectral Tomography},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Operational central states, exponent no-go, symmetry-reference costs, moment bounds, formal cores, and heat-kernel spectral tomography},
  pdfkeywords={FIN, spectral theory, operational tomography, heat kernel, bistochastic state, asymmetry resource, fidelity, e-process, trigonometric moment problem}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs 217--229}
\fancyhead[R]{\small Release 10.20}
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
{\Large\bfseries FIN Research Monograph --- Release 10.20\par}
\vspace{10mm}
{\Huge\bfseries FIN Programs 217--229\par}
\vspace{6mm}
{\Large Operational Sections, Reference Costs,\\
and Heat-Kernel Spectral Tomography\par}
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
\textbf{Scope.} Thirteen executed studies of operational central states,
natural exponent flows, formal certification, finite symmetry-reference
costs, dependence-aware inference, process tomography, environmental moments,
phase formalization, heat-kernel spectral reconstruction, independent
custody, and external prediction gating.

\medskip
\textbf{Central result.} Exact heat transition data identify the finite
generator by \(A=-\tau^{-1}\log P_\tau\), while projective eigenvalue ratios
remain identifiable without an absolute clock unit. Imperfect symmetry
references cannot be returned in uncorrelated product form after broadcasting
nonzero asymmetry.

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
This release reports Programs 217--229 and recommends Programs 230--242.
Analytic proofs, exact arithmetic, machine-checked finite fragments,
simulation evidence, conditional instruments, unavailable infrastructure,
and failed external-data gates remain separate.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""

POSTAMBLE = r"""
\clearpage
\chapter*{Release reproducibility}
\addcontentsline{toc}{chapter}{Release reproducibility}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/matplotlib-fin217 \
python3 fin_programs_217_229_operational_section_spectral_tomography.py

MPLCONFIGDIR=/tmp/matplotlib-fin217 \
python3 -m unittest -v \
  test_fin_programs_217_229_operational_section_spectral_tomography.py
\end{Verbatim}

Machine-readable results are stored in
\path{FIN_Programs_217_229_Operational_Section_Spectral_Tomography_Results.json}.
The release contains five auxiliary contracts or formal records, two
machine-compiled Lean core sources, 82 tests, and thirteen figures.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Programs 217--229: Operational Sections,
Reference Costs, and Heat-Kernel Spectral Tomography}
(FIN Research Monograph, Release 10.20; Version 1.0.0) [Preprint]. Zenodo.
\end{document}
"""

FIGURES = {
    "program217_operational_central_state.png": "Bistochastic naturality selects the uniform central state, while a preparation record can select and estimate a different state.",
    "program218_affine_eta_no_go.png": r"Every allowed affine-natural flow has \(F(4/5)=4\kappa/5\); only the zero flow fixes the target and then fixes every point.",
    "program219_arb_execution.png": "The immutable formal execution gate remains closed because no locally callable Arb/FLINT engine is available.",
    "program220_lean_build.png": "The dependency-free four-cycle certificate compiles; the general Mathlib graph library remains dependency-blocked.",
    "program221_reference_cost.png": "Imperfect product return violates fidelity monotonicity, whereas marginal return stores a measurable target--reference correlation.",
    "program222_robust_mixing_design.png": "A calibrated lower refresh bound and exact budget optimization determine the robust thinning design.",
    "program223_reusable_eprocess.png": "Sequential insertion ranks support an anytime-valid rank-filtration e-process while fixed calibration reuse lacks the same conditional law.",
    "program224_likelihood_region.png": "The likelihood-quadratic region is substantially narrower than the exact conservative region, with its asymptotic boundary kept explicit.",
    "program225_higher_moment.png": r"Toeplitz positivity gives the sharp interval \(-11/180\le c_3\le11/20\).",
    "program226_phase_formalization.png": "The finite torsion orbit is machine checked; automatic continuity and the transcendence theorem remain analytic.",
    "program227_spectral_tomography.png": "Finite heat-kernel transition counts reconstruct the projective strict spectral fingerprint with shot-dependent power.",
    "program228_external_custody.png": "No local process or event-level double-slit candidate passes the independent eleven-field custody gate.",
    "program229_prediction_lock.png": r"The held-out semigroup prediction \(P_{2\tau}=P_\tau^2\) remains locked pending independent data.",
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
            r"!\[[^\]]*\]\(FIN_Programs_217_229_Operational_Section_Spectral_Tomography_Figures/"
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
        "The report distinguishes analytic proof, exact arithmetic, "
        "machine-checked finite fragments, simulation evidence, conditional "
        "operational construction, and open infrastructure or data gates. "
        "No finite example or unavailable external test is promoted.\n\n"
        + text[abstract_at:]
    )
    body = converter.body_from_markdown(text)
    body = body.replace("Ż", r"\.Z").replace("ż", r"\.z")
    for filename, caption in FIGURES.items():
        token = f"FINFIG\\_{filename.replace('_', r'\_')}"
        path = "FIN_Programs_217_229_Operational_Section_Spectral_Tomography_Figures/" + filename
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
