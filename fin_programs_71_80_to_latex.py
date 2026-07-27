#!/usr/bin/env python3
"""Build archival LaTeX for FIN Programs 71--80."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Programs_71_80_Symbol_Operational_Bridge_Monograph.md")
TARGET = Path("FIN_Programs_71_80_Symbol_Operational_Bridge_Monograph.tex")

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
  pdftitle={FIN Programs 71--80: Symbol and Operational Bridge},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Fourier-symbol renormalization, quotient information geometry, and operational interfaces},
  pdfkeywords={FIN, Schur complement, Fourier symbol, information geometry, Landauer, process tomography}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs 71--80}
\fancyhead[R]{\small Release 10.8}
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
\newcommand{\statusStrong}{\textcolor{finblue}{[Strong evidence]}}
\newcommand{\statusModerate}{\textcolor{finviolet}{[Moderate evidence]}}
\newcommand{\statusConditional}{\textcolor{finorange}{[Conditional]}}
\newcommand{\statusRefuted}{\textcolor{finred}{[Refuted]}}
\newcommand{\statusOpen}{\textcolor{finviolet}{[Open]}}

\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{8mm}
{\Large\bfseries FIN Research Monograph --- Release 10.8\par}
\vspace{10mm}
{\Huge\bfseries FIN Programs 71--80\par}
\vspace{6mm}
{\Large Fourier-Symbol Renormalization, Quotient Information Geometry,\\
and the Operational Bridge to Physics\par}
\vspace{16mm}
{\Large Krzysztof \.Zuchowski\par}
\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
\vspace{2mm}
{\normalsize ORCID: \href{https://orcid.org/0009-0002-0909-3613}{0009-0002-0909-3613}\par}
\vspace{9mm}
{\large 27 July 2026\par}
{\normalsize Publication --- Preprint; record version 1.0.0\par}
\vfill
\begin{minipage}{0.88\textwidth}
\small
\textbf{Scope.} Ten executed programs on Fourier-symbol Schur elimination,
large-\(N\) naturality, quotient information geometry, locality, dimensional
calibration, state-dependent chirality, finite-time information
thermodynamics, noisy process tomography, a fixed legacy--strict bridge
falsification, and external-data preregistration.

\medskip
\textbf{Guardrail.} No strict selector, physical-unit source, causal/Lorentz
closure, legacy physical-role transfer, full legacy--strict completion, or
Theory-of-Everything closure is claimed.

\medskip
\textbf{Ontology.} The nadsoliton is retained as primordial fractal
information in a solitonic state; no lower information layer is introduced.

\medskip
\textbf{License.} Creative Commons Attribution 4.0 International (CC BY 4.0).
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}

\pagenumbering{roman}
\chapter*{Publication statement}
\addcontentsline{toc}{chapter}{Publication statement}
This self-contained release reports Programs 71--80 and recommends Programs
81--92. Analytic theorems, finite numerical evidence, conditional operational
models, failed bridges, and open obligations are labeled separately.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""

POSTAMBLE = r"""
\clearpage
\chapter*{Archival note}
\addcontentsline{toc}{chapter}{Archival note}
Machine-readable results are stored in
\path{FIN_Programs_71_80_Symbol_Operational_Bridge_Results.json}. The external
protocol is frozen in
\path{FIN_Programs_71_80_External_Data_Preregistration.json}. The executable,
twelve tests, and six figures reproduce the finite results.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Programs 71--80: Fourier-Symbol
Renormalization, Quotient Information Geometry, and the Operational Bridge to
Physics} (FIN Research Monograph, Release 10.8; Version 1.0.0) [Preprint].
\end{document}
"""

FIGURES = {
    "FINFIG72": (
        "FIN_Programs_71_80_Symbol_Operational_Bridge_Figures/program72_large_n_projective_scaling.png",
        "Large-N native-versus-Schur naturality defects under two explicitly distinct distance semantics.",
    ),
    "FINFIG73": (
        "FIN_Programs_71_80_Symbol_Operational_Bridge_Figures/program73_zero_mode_quotient_kl.png",
        r"Convergence of the shifted Gaussian relative entropy to the regulator-free quotient on \(\mathbf 1^\perp\).",
    ),
    "FINFIG74": (
        "FIN_Programs_71_80_Symbol_Operational_Bridge_Figures/program74_locality_truncation.png",
        "Operator fidelity and leading propagation order under strict-kernel range truncation.",
    ),
    "FINFIG77": (
        "FIN_Programs_71_80_Symbol_Operational_Bridge_Figures/program77_landauer_jarzynski.png",
        "Finite-time erasure work and numerical verification of the Jarzynski identity.",
    ),
    "FINFIG78": (
        "FIN_Programs_71_80_Symbol_Operational_Bridge_Figures/program78_process_tomography.png",
        "Shot-noise-limited discrimination of wave and diffusion records in a short-time regime.",
    ),
    "FINFIG79": (
        "FIN_Programs_71_80_Symbol_Operational_Bridge_Figures/program79_legacy_strict_schur_flow.png",
        "Fixed no-fit repaired-legacy Schur flow versus the native repaired-strict family.",
    ),
}


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    # The local LuaTeX installation falls back to T1 Latin Modern.  Normalize
    # typographic punctuation before conversion so every glyph is embedded.
    text = (
        text.replace("–", "--")
        .replace("—", "---")
        .replace("’", "'")
        .replace("“", '"')
        .replace("”", '"')
    )
    replacements = {
        "![Large-N projective scaling](FIN_Programs_71_80_Symbol_Operational_Bridge_Figures/program72_large_n_projective_scaling.png)": "FINFIG72",
        "![Zero-mode quotient KL](FIN_Programs_71_80_Symbol_Operational_Bridge_Figures/program73_zero_mode_quotient_kl.png)": "FINFIG73",
        "![Locality truncation](FIN_Programs_71_80_Symbol_Operational_Bridge_Figures/program74_locality_truncation.png)": "FINFIG74",
        "![Landauer and Jarzynski](FIN_Programs_71_80_Symbol_Operational_Bridge_Figures/program77_landauer_jarzynski.png)": "FINFIG77",
        "![Noisy process tomography](FIN_Programs_71_80_Symbol_Operational_Bridge_Figures/program78_process_tomography.png)": "FINFIG78",
        "![Legacy-to-strict Schur flow](FIN_Programs_71_80_Symbol_Operational_Bridge_Figures/program79_legacy_strict_schur_flow.png)": "FINFIG79",
    }
    for markdown, token in replacements.items():
        text = text.replace(markdown, f"\n\n{token}\n\n")
    text = re.sub(r"```[a-zA-Z]*\n(.*?)```", r"\1", text, flags=re.DOTALL)
    text = re.sub(r"`([^`]+)`", r"\1", text)
    text = re.sub(
        r"^(#{1,4})\s+\d+(?:\.\d+)*\.\s+",
        r"\1 ",
        text,
        flags=re.MULTILINE,
    )

    body = converter.body_from_markdown(text)
    body = body.replace("Ż", r"\.Z").replace("ż", r"\.z")
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
        ("[Open]", r"\statusOpen{}"),
    ]:
        body = body.replace(plain, macro)

    for filename in [
        "fin_programs_71_80_symbol_operational_bridge.py",
        "test_fin_programs_71_80_symbol_operational_bridge.py",
        "FIN_Programs_71_80_Symbol_Operational_Bridge_Results.json",
        "FIN_Programs_71_80_External_Data_Preregistration.json",
        "FIN_Programs_71_80_Symbol_Operational_Bridge_Figures/",
    ]:
        body = body.replace(
            converter.escape_normal(filename), f"\\path{{{filename}}}"
        )

    digest = "1a10b804a3e8265969fa66626f5152acc00830440841f2ddee6a8a59c018db4e"
    digest_tex = (
        "\\texttt{"
        + "}\\allowbreak\\texttt{".join(
            digest[i : i + 16] for i in range(0, len(digest), 16)
        )
        + "}"
    )
    body = body.replace(digest, digest_tex)

    TARGET.write_text(PREAMBLE + "\n" + body + "\n" + POSTAMBLE, encoding="utf-8")
    print(f"Wrote {TARGET} ({TARGET.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
