#!/usr/bin/env python3
"""Build the archival PDF source for the post-Programs-41--50 correction."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Post_41_50_Methodology_Correction_and_Mirror_Coupling.md")
TARGET = Path("FIN_Post_41_50_Methodology_Correction_and_Mirror_Coupling.tex")

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
\definecolor{fingray}{HTML}{4D4D4D}

\hypersetup{
  unicode=true,
  colorlinks=true,
  linkcolor=finblue,
  citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Post-41--50 Scientific Correction},
  pdfauthor={Krzysztof \.Zuchowski},
  pdfsubject={Methodology audit, mirror-coupling theorem, and canonical kernel policy},
  pdfkeywords={FIN, strict kernel, legacy kernel, mirror coupling, symmetry breaking, selector}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Post-41--50 Scientific Correction}
\fancyhead[R]{\small Release 10.5}
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
\vspace*{10mm}
{\Large\bfseries FIN Research Supplement --- Release 10.5\par}
\vspace{11mm}
{\Huge\bfseries FIN Post-41--50 Scientific Correction\par}
\vspace{6mm}
{\Large Methodology Audit, Mirror-Coupling Theorem, and Canonical Kernel Policy\par}
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
\textbf{Scope.} Scientific correction of the post-Programs-31--40 research,
including Programs 41--50 and P3172, followed by an executable mirror-coupling
audit and a canonical strict/legacy kernel policy.

\medskip
\textbf{Guardrail.} No silent legacy--strict identification, no role transfer,
no strict selector closure, no physical unit export, and no
Theory-of-Everything claim.

\medskip
\textbf{Reproducibility.}
\path{fin_post_41_50_correction_and_mirror.py},
\path{FIN_Post_41_50_Corrected_Results.json}, and the three generated figures.

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
The numerical record \path{FIN_Post_41_50_Corrected_Results.json} is generated
by \path{fin_post_41_50_correction_and_mirror.py}. Corrected Programs 41--50 are
recorded by \path{fin_programs_41_50_legacy_star.py}. The exact P3172 threshold
and full metric-axiom audit are independently covered by its updated tests.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Post-41--50 Scientific Correction:
Methodology Audit, Mirror-Coupling Theorem, and Canonical Kernel Policy}
(FIN Research Supplement, Release 10.5; Version 1.0.0) [Preprint].

\end{document}
"""

FIGURES = {
    "FINFIGKERNELPOLICY": (
        "FIN_Post_41_50_Correction_Figures/kernel_policy_profiles.png",
        "Canonical strict and legacy profiles together with the two noncanonical legacy-family benchmark freezes.",
    ),
    "FINFIGPSD": (
        "FIN_Post_41_50_Correction_Figures/psd_threshold_support_scan.png",
        "The Gram-positive threshold depends on the finite cycle support; it is not a universal constant.",
    ),
    "FINFIGMIRROR": (
        "FIN_Post_41_50_Correction_Figures/mirror_bias_current_entropy.png",
        "Directional current is odd in the mirror-bias parameter, whereas entropy production is even.",
    ),
}


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    replacements = {
        "![Kernel policy profiles](FIN_Post_41_50_Correction_Figures/kernel_policy_profiles.png)": "FINFIGKERNELPOLICY",
        "![PSD support scan](FIN_Post_41_50_Correction_Figures/psd_threshold_support_scan.png)": "FINFIGPSD",
        "![Mirror current](FIN_Post_41_50_Correction_Figures/mirror_bias_current_entropy.png)": "FINFIGMIRROR",
    }
    for markdown, token in replacements.items():
        text = text.replace(markdown, f"\n\n{token}\n\n")
    text = re.sub(r"```[a-zA-Z]*\n(.*?)```", r"\1", text, flags=re.DOTALL)
    text = re.sub(r"`([^`]+)`", r"\1", text)
    # LaTeX supplies section numbering; remove the manual Markdown prefixes.
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
    body = body.replace(
        "\\[\n\\boxed{\n"
        "\\text{it constructs the correct inversion-odd carrier, but not the law selecting its sign.}\n"
        "}\n\\]",
        "\\begin{center}\n"
        "\\fbox{\\parbox{0.86\\textwidth}{\\centering "
        "It constructs the correct inversion-odd carrier, but not the law selecting its sign.}}\n"
        "\\end{center}",
    )
    body = body.replace(
        "\\[\n\\boxed{\n"
        "K_{\\rm legacy\\_ont}\n"
        "\\quad\\text{(canonical intermediate bridge object)}\n"
        "\\qquad\\text{and}\\qquad\n"
        "K_{\\rm strict\\_gate}\n"
        "\\quad\\text{(primary strict working object)}.\n"
        "}\n\\]",
        "\\[\n\\boxed{\\begin{gathered}\n"
        "K_{\\rm legacy\\_ont}\\quad\\text{(canonical intermediate bridge object)},\\\\\n"
        "K_{\\rm strict\\_gate}\\quad\\text{(primary strict working object)}.\n"
        "\\end{gathered}}\n\\]",
    )
    # Let url.sty insert line breaks in archival filenames.
    for filename in [
        "fin_programs_41_50_legacy_star.py",
        "FIN_Programs_41_50_Legacy_Star_Results.json",
        "fundamental_action_reconstruction/p3172_s2122_legacy_star_operator_model_generator_potential_audit.py",
        "fin_post_41_50_correction_and_mirror.py",
        "FIN_Post_41_50_Corrected_Results.json",
        "FIN_Post_41_50_Correction_Figures/",
    ]:
        escaped = converter.escape_normal(filename)
        body = body.replace(escaped, f"\\path{{{filename}}}")

    TARGET.write_text(PREAMBLE + "\n" + body + "\n" + POSTAMBLE, encoding="utf-8")
    print(f"Wrote {TARGET} ({TARGET.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
