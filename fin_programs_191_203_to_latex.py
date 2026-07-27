#!/usr/bin/env python3
"""Build archival LaTeX for FIN Programs 191--203."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Programs_191_203_Reference_Process_Prediction_Monograph.md")
TARGET = Path("FIN_Programs_191_203_Reference_Process_Prediction_Monograph.tex")

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
  pdftitle={FIN Programs 191--203: Reference States, Operational Quotients, and Conditional Prediction},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Reference-state naturality, scale quotients, reflection conversion, process tomography, environment identifiability and conditional prediction},
  pdfkeywords={FIN, operator algebra, spectral theory, Dirichlet form, Alberti-Uhlmann, process tensor, mixing, operational physics, scale quotient}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs 191--203}
\fancyhead[R]{\small Release 10.18}
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
{\Large\bfseries FIN Research Monograph --- Release 10.18\par}
\vspace{10mm}
{\Huge\bfseries FIN Programs 191--203\par}
\vspace{6mm}
{\Large Reference States, Operational Quotients,\\
and Conditional Prediction\par}
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
\textbf{Scope.} Thirteen executed studies of reference-state naturality,
compression gauge, rigorous certification, finite Dirichlet dynamics,
reflection-covariant conversion, mixing-aware inference, temporal open-set
detection, process tomography, environment equivalence, analytic phase
provenance, scale-free observables, experimental intake and one conditional
prediction.

\medskip
\textbf{Central result.} The qubit reflection-conversion boundary is reduced
to a finite closed formula, and single/plus/echo process tomography is proven
minimal for an unknown-blur Gaussian two-step model. States on non-simple
algebras, positive scale and microscopic environments remain nontrivial
moduli or quotient data.

\medskip
\textbf{Guardrail.} No QW-2191 discharge, strict selector, physical unit,
strict beta source, completed legacy--strict bridge, role transfer,
\(L_{\rm total}\), Theory-of-Everything closure or external physical
validation is claimed.

\medskip
\textbf{License.} Creative Commons Attribution 4.0 International (CC BY 4.0).
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}

\pagenumbering{roman}
\chapter*{Publication statement}
\addcontentsline{toc}{chapter}{Publication statement}
This release reports Programs 191--203 and recommends Programs 204--216.
Analytic proofs, exact finite audits, synthetic evidence, conditional
constructions, unavailable proof toolchains and external-data failures are
kept separate.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""

POSTAMBLE = r"""
\clearpage
\chapter*{Release reproducibility}
\addcontentsline{toc}{chapter}{Release reproducibility}
\begin{Verbatim}[fontsize=\small]
python3 fin_programs_191_203_reference_process_prediction.py
python3 -m unittest -v \
  test_fin_programs_191_203_reference_process_prediction.py
\end{Verbatim}

Machine-readable results are stored in
\path{FIN_Programs_191_203_Reference_Process_Prediction_Results.json}.
The release also contains the frozen order-sensitive protocol, the synthetic
operational dry-run bundle, the finite Lean source, sixty-five tests and
thirteen figures.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Programs 191--203: Reference States,
Operational Quotients, and Conditional Prediction}
(FIN Research Monograph, Release 10.18; Version 1.0.0) [Preprint]. Zenodo.
\end{document}
"""

FIGURES = {
    "program191_reference_state_naturality.png": "Inner-unitary naturality fixes normalized trace within blocks but leaves a central probability measure.",
    "program192_beta_gauge_quotient.png": r"Quotienting the positive \(\beta\) torsor removes coordinate scale while preserving the strict--legacy shape mismatch.",
    "program193_common_engine_certificate.png": "The common directed-rounding certificate remains blocked by missing engine closure for four components.",
    "program194_dirichlet_library.png": "Exact finite Dirichlet mathematics and numerical semigroup residuals, with proof-assistant compilation kept separate.",
    "program195_reflection_geometry.png": "Closed analytic reachability boundary for reflection-covariant conversion from the source state with transverse coordinate 0.6.",
    "program196_regenerative_mixing.png": "Observed regeneration supplies a valid iid subsample at a measurable information cost.",
    "program197_order_sensitive_open_set.png": "The frozen temporal feature protocol rejects four declared order-dependent alternatives while retaining iid specificity.",
    "program198_minimal_process_tomography.png": "The single, plus and echo log-visibility design has rank three and identifies blur, variance and covariance.",
    "program199_environment_equivalence.png": "Two distinct phase environments generate the same one-time dephasing channel but different two-time characteristics.",
    "program200_analytic_phase_sources.png": "Natural circle-group endomorphism images of the legacy phase miss the strict transcendental phase.",
    "program201_scale_free_observables.png": "The raw gap changes over twelve decades while projective heat and resolvent observables remain invariant.",
    "program202_operational_bundle.png": "The synthetic method-validation bundle passes ten of eleven intake fields and fails independent provenance.",
    "program203_conditional_prediction.png": "A held-out visibility is predicted in an explicit W0 plus CA plus OP synthetic dry run.",
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
            r"!\[[^\]]*\]\(FIN_Programs_191_203_Reference_Process_Prediction_Figures/"
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
        "synthetic evidence, conditional construction and open obligations. "
        "No unavailable certificate is promoted by prose.\n\n"
        + text[abstract_at:]
    )
    body = converter.body_from_markdown(text)
    body = body.replace("Ż", r"\.Z").replace("ż", r"\.z")
    for filename, caption in FIGURES.items():
        token = f"FINFIG\\_{filename.replace('_', r'\_')}"
        path = "FIN_Programs_191_203_Reference_Process_Prediction_Figures/" + filename
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
