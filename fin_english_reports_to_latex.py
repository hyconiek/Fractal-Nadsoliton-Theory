#!/usr/bin/env python3
"""Build archival LaTeX sources for the two English FIN reports."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


COMMON_PACKAGES = r"""
\documentclass[11pt,a4paper,openany]{report}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
\usepackage{microtype}
\usepackage{amsmath,amssymb,amsthm,mathtools,bm}
\usepackage{booktabs,array,longtable,tabularx}
\usepackage{xcolor}
\usepackage[margin=22mm,headheight=23pt]{geometry}
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
"""

COMMON_SETUP = r"""
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
\newcommand{\statusWeak}{\textcolor{fingray}{[Weak evidence]}}
\newcommand{\statusSpeculative}{\textcolor{finorange}{[Speculative]}}
\newcommand{\statusRefuted}{\textcolor{finred}{[Refuted]}}
"""


ATLAS_SOURCE = Path("FIN_Second_Generation_Atlas_Derived_Nadsoliton_Objects_EN.md")
ATLAS_TARGET = Path("FIN_Second_Generation_Atlas_Derived_Nadsoliton_Objects_EN.tex")
ATLAS_FIGURES = {
    "combination_search_coverage.png": (
        "Complete coverage of singles, pairs, and triples, together with the "
        "output of the explicit construction grammar."
    ),
    "derived_object_score_matrix.png": (
        "Audit scores for O01--O15. Overinterpretation risk has the opposite "
        "semantics from the other scores."
    ),
    "generation2_generation3_map.png": (
        "Input dependencies between the O and G generations. Arrows do not "
        "assert physical equivalence."
    ),
}

PROGRAMS_SOURCE = Path("FIN_Programs_255_266_Research_Report_EN.md")
PROGRAMS_TARGET = Path("FIN_Programs_255_266_Research_Report_EN.tex")
PROGRAMS_FIGURES = {
    "p255_p258_stieltjes_chiral.png": (
        "Stieltjes-function limits and the norm of chiral memory "
        "susceptibility with its proved upper bound."
    ),
    "p259_p260_rg_information_ledger.png": (
        "Strictly decreasing context flow and the multitime "
        "distinguishability balance."
    ),
    "p261_p262_bridge_fingerprint.png": (
        "Amplitude-only completion obstruction and the synthetic "
        "fingerprint-first/calibration-second audit."
    ),
    "p264_false_positive_atlas.png": (
        "Pass rates of five frozen tests in the audited null ensembles."
    ),
    "p266_reservoir_benchmark.png": (
        "The FIN reservoir compared with 120 dimension- or spectrum-matched "
        "controls."
    ),
}

PUZZLE_SOURCE = Path("FIN_Nadsoliton_Z12_Puzzle_Atlas_Report_EN.md")
PUZZLE_TARGET = Path("FIN_Nadsoliton_Z12_Puzzle_Atlas_Report_EN.tex")
PUZZLE_FIGURES = {
    "single_puzzle_domain_atlas.png": (
        "Heuristic matching matrix for 12 puzzles and 15 fields. Scores are "
        "used only to search for analogies."
    ),
    "dynamic_schur_memory.png": (
        "Composition defect of reduced propagators and the error incurred by "
        "replacing dynamic reduction with one static Schur generator."
    ),
}


def common_body(text: str, figures: dict[str, str], figure_dir: str) -> str:
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
        .replace("→", "->")
        .replace("↔", "<->")
    )
    for filename in figures:
        pattern = re.compile(
            r"!\[[^\]]*\]\(" + re.escape(figure_dir + "/" + filename) + r"\)"
        )
        text = pattern.sub(f"\n\nFINFIG_{filename}\n\n", text)
    text = re.sub(r"`([^`]+)`", r"\1", text)
    body = converter.body_from_markdown(text)
    body = body.replace(
        r"\begingroup\footnotesize",
        r"\begingroup\fontsize{7}{8.4}\selectfont",
    )
    body = body.replace(
        r"\setlength{\tabcolsep}{3.5pt}",
        r"\setlength{\tabcolsep}{2.5pt}",
    )
    body = body.replace("Ż", r"\.Z").replace("ż", r"\.z")
    for filename, caption in figures.items():
        token = f"FINFIG\\_{filename.replace('_', r'\_')}"
        path = f"{figure_dir}/{filename}"
        body = body.replace(
            token,
            "\\begin{figure}[htbp]\n"
            "\\centering\n"
            f"\\includegraphics[width=0.97\\textwidth]{{{path}}}\n"
            f"\\caption{{{caption}}}\n"
            "\\end{figure}",
        )
    return body


def atlas_preamble() -> str:
    return COMMON_PACKAGES + r"""
\hypersetup{
  unicode=true,colorlinks=true,linkcolor=finblue,citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={Second-Generation Atlas: Derived Nadsoliton Objects},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Systematic search over FIN combinations C01--C12, objects O01--O15, and synthesis G01--G10},
  pdfkeywords={FIN, nadsoliton, Stieltjes functions, Schur complement, memory kernel, operator theory, information geometry}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small Second-Generation Atlas}
\fancyhead[R]{\small Release 10.23 --- English edition}
\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{12mm}
{\Large\bfseries FIN --- Release 10.23\par}
\vspace{12mm}
{\Huge\bfseries Second-Generation Atlas\par}
\vspace{5mm}
{\LARGE\bfseries Derived Nadsoliton Objects\par}
\vspace{8mm}
{\Large Systematic search over C01...C12,\\
construction of O01...O15, and synthesis of G01...G10\par}
\vspace{18mm}
{\Large Krzysztof \.Zuchowski\par}
\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
\vspace{2mm}
{\normalsize ORCID: \href{https://orcid.org/0009-0002-0909-3613}{0009-0002-0909-3613}\par}
\vspace{10mm}
{\large 28 July 2026\par}
\vfill
\begin{minipage}{0.9\textwidth}
\small
\textbf{Central object.}
\[
E\longmapsto \operatorname{Schur}_{V\setminus E}(zI+A),\qquad
\Sigma_E(z)=A_{EH}(zI+A_{HH})^{-1}A_{HE}.
\]
Contextual reductions compose, while self-energies are positive
operator-valued Stieltjes functions generating exact memory.

\medskip
\textbf{Boundary.} This report does not export a selector, physical scale,
legacy--strict bridge, role transfer, or laboratory validation.

\medskip
\textbf{Repository.}
\url{https://github.com/hyconiek/Fractal-Nadsoliton-Theory}

\textbf{License.} CC BY 4.0.
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}
\pagenumbering{roman}
\tableofcontents
\clearpage
\pagenumbering{arabic}
"""


def programs_preamble() -> str:
    return COMMON_PACKAGES + r"""
\hypersetup{
  unicode=true,colorlinks=true,linkcolor=finblue,citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Programs 255--266: Operator Memory Measures, Identifiability, Current Tomography and Adaptive-Analogy Falsification},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Executed FIN research Programs 255--266 and recommended Programs 267--280},
  pdfkeywords={FIN, operator-valued Stieltjes functions, Schur complement, memory kernels, process information, current tomography, identifiability, reservoir computing}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs 255--266}
\fancyhead[R]{\small Release 10.24 --- English edition}
\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{11mm}
{\Large\bfseries FIN --- Release 10.24\par}
\vspace{10mm}
{\Huge\bfseries Research Programs 255--266\par}
\vspace{5mm}
{\LARGE\bfseries Operator Memory Measures, Identifiability,\\
Current Tomography, and Falsification of Adaptive Analogies\par}
\vspace{16mm}
{\Large Krzysztof \.Zuchowski\par}
\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
\vspace{2mm}
{\normalsize ORCID: \href{https://orcid.org/0009-0002-0909-3613}{0009-0002-0909-3613}\par}
\vspace{9mm}
{\large 28 July 2026\par}
\vfill
\begin{minipage}{0.9\textwidth}
\small
\textbf{Central result.}
For every nontrivial strict \(Z_{12}\) context,
\[
\Sigma_E(z)=\int_{(0,\infty)}\frac{d\mathsf M_E(\mu)}{z+\mu},
\qquad d\mathsf M_E(\mu)\succeq0.
\]
Memory determines a minimal input--output realization class, not a unique
microscopic sector.

\medskip
\textbf{Boundary.}
This report does not export a selector, physical scale, complete
legacy--strict bridge, adaptive law, or laboratory validation.

\medskip
\textbf{Repository.}
\url{https://github.com/hyconiek/Fractal-Nadsoliton-Theory}

\textbf{License.} CC BY 4.0.
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}
\pagenumbering{roman}
\tableofcontents
\clearpage
\pagenumbering{arabic}
"""


def puzzle_preamble() -> str:
    return COMMON_PACKAGES + r"""
\hypersetup{
  unicode=true,colorlinks=true,linkcolor=finblue,citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={Nadsoliton Puzzle Atlas: Z12 Audit and Dynamic Reduction},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Controlled cross-disciplinary FIN analysis, Z12 simulator audit, and a theorem on memory generated by dynamic Schur reduction},
  pdfkeywords={FIN, Z12, spectral graph, Schur complement, memory kernel, Mori-Zwanzig, chirality, operator theory}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small Nadsoliton Puzzle Atlas}
\fancyhead[R]{\small Release 10.22 --- English edition}
\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{12mm}
{\Large\bfseries FIN --- Release 10.22\par}
\vspace{12mm}
{\Huge\bfseries Nadsoliton Puzzle Atlas\par}
\vspace{6mm}
{\Large Z12 Simulator Audit, Controlled Cross-Disciplinary Intuition,\\
and a New Memory Object after Reduction\par}
\vspace{18mm}
{\Large Krzysztof \.Zuchowski\par}
\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
\vspace{2mm}
{\normalsize ORCID: \href{https://orcid.org/0009-0002-0909-3613}{0009-0002-0909-3613}\par}
\vspace{10mm}
{\large 28 July 2026\par}
\vfill
\begin{minipage}{0.89\textwidth}
\small
\textbf{Central result.} Exact elimination of hidden strict \(Z_{12}\) nodes
generates a frequency-dependent self-energy and a process with memory.
The static Schur complement is only a resolvent limit, not the exact generator
of the complete reduced dynamics.

\medskip
\textbf{Boundary.} The result does not generate a selector, physical units,
a legacy--strict bridge, or experimental validation.

\medskip
\textbf{Repository.}
\url{https://github.com/hyconiek/Fractal-Nadsoliton-Theory}

\textbf{License.} CC BY 4.0.
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}
\pagenumbering{roman}
\tableofcontents
\clearpage
\pagenumbering{arabic}
"""


def atlas_postamble() -> str:
    return r"""
\clearpage
\chapter*{Reproduction certificate}
\addcontentsline{toc}{chapter}{Reproduction certificate}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/matplotlib-second \
python3 fin_nadsoliton_second_generation_atlas.py
\end{Verbatim}
Expected coverage: 305 combinations, 15 O-generation objects, and 10
G-generation objects. The construction grammar is a candidate filter, not
proof of equivalence or physical realization.
\end{document}
"""


def programs_postamble() -> str:
    return r"""
\clearpage
\chapter*{Reproduction certificate}
\addcontentsline{toc}{chapter}{Reproduction certificate}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/matplotlib-p255-266 \
python3 fin_programs_255_266.py

MPLCONFIGDIR=/tmp/matplotlib-p255-266 \
python3 -m unittest -v test_fin_programs_255_266.py
\end{Verbatim}
Expected result: 15 passed, 0 failed. P262, P264, P265, and P266 use synthetic
data and are not laboratory records.
\end{document}
"""


def puzzle_postamble() -> str:
    return r"""
\clearpage
\chapter*{Reproduction certificate}
\addcontentsline{toc}{chapter}{Reproduction certificate}
\begin{Verbatim}[fontsize=\small]
python3 audit_z12_sim.py
MPLCONFIGDIR=/tmp/matplotlib-puzzle \
python3 fin_nadsoliton_puzzle_atlas.py
\end{Verbatim}
Expected Z12 audit result: 15/15 tests passed. Atlas scores are search
heuristics; operator theorems are verified with separate residuals.
\end{document}
"""


def build_atlas() -> None:
    text = ATLAS_SOURCE.read_text(encoding="utf-8")
    text = re.sub(
        r"^(#{1,4})\s+\d+(?:\.\d+)*\.\s+",
        lambda match: match.group(1) + " ",
        text,
        flags=re.MULTILINE,
    )
    text = re.sub(
        r"^(#{2,4})(\s+)",
        lambda match: "#" * (len(match.group(1)) - 1) + match.group(2),
        text,
        flags=re.MULTILINE,
    )
    text = text.replace("# Confidence convention", "## Confidence convention", 1)
    body = common_body(
        text, ATLAS_FIGURES, "FIN_Second_Generation_Atlas_Figures"
    )
    body = re.sub(
        r"(C\d{2})\+(?=C\d{2})",
        lambda match: match.group(1) + r"+\allowbreak{}",
        body,
    )
    body = body.replace(
        r"\section{Confidence convention}",
        r"\chapter*{Confidence convention}"
        "\n"
        r"\addcontentsline{toc}{chapter}{Confidence convention}",
    )
    ATLAS_TARGET.write_text(
        atlas_preamble() + "\n" + body + "\n" + atlas_postamble(),
        encoding="utf-8",
    )


def build_programs() -> None:
    text = PROGRAMS_SOURCE.read_text(encoding="utf-8")
    text = re.sub(
        r"^(#{1,4})\s+\d+(?:\.\d+)*\.\s+",
        lambda match: match.group(1) + " ",
        text,
        flags=re.MULTILINE,
    )
    body = common_body(text, PROGRAMS_FIGURES, "FIN_Programs_255_266_Figures")
    body = body.replace(
        r"\section{Confidence convention}",
        r"\chapter*{Confidence convention}"
        "\n"
        r"\addcontentsline{toc}{chapter}{Confidence convention}",
    )
    PROGRAMS_TARGET.write_text(
        programs_preamble() + "\n" + body + "\n" + programs_postamble(),
        encoding="utf-8",
    )


def build_puzzle() -> None:
    text = PUZZLE_SOURCE.read_text(encoding="utf-8")
    text = re.sub(
        r"^(#{1,4})\s+\d+(?:\.\d+)*\.\s+",
        lambda match: match.group(1) + " ",
        text,
        flags=re.MULTILINE,
    )
    text = re.sub(
        r"^(#{2,4})(\s+)",
        lambda match: "#" * (len(match.group(1)) - 1) + match.group(2),
        text,
        flags=re.MULTILINE,
    )
    text = text.replace("# Confidence convention", "## Confidence convention", 1)
    body = common_body(
        text, PUZZLE_FIGURES, "FIN_Nadsoliton_Puzzle_Atlas_Figures"
    )
    body = body.replace(
        r"\section{Confidence convention}",
        r"\chapter*{Confidence convention}"
        "\n"
        r"\addcontentsline{toc}{chapter}{Confidence convention}",
    )
    PUZZLE_TARGET.write_text(
        puzzle_preamble() + "\n" + body + "\n" + puzzle_postamble(),
        encoding="utf-8",
    )


def main() -> None:
    build_atlas()
    build_programs()
    build_puzzle()
    print(ATLAS_TARGET)
    print(PROGRAMS_TARGET)
    print(PUZZLE_TARGET)


if __name__ == "__main__":
    main()
