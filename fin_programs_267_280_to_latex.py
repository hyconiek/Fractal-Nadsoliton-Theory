#!/usr/bin/env python3
"""Build the archival English PDF source for FIN Programs P267--P280."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body


SOURCE = Path("FIN_Current_State_After_Programs_267_280_EN.md")
TARGET = Path("FIN_Current_State_After_Programs_267_280_EN.tex")

FIGURES = {
    "p267_p269_measure_povm.png": (
        "Visible Stieltjes atoms, the severe noise sensitivity of inverse "
        "recovery, and the six-outcome current POVM."
    ),
    "p270_p272_rg_completion.png": (
        "Rank-defect obstruction for linear size-restoring RG and the "
        "positive-shift legacy completion test."
    ),
    "p274_p275_tomography_false_positive.png": (
        "Finite-flux bias--noise optimization and deterministic robustness "
        "certificates for the frozen null fingerprints."
    ),
    "p277_p278_adaptation_reservoir.png": (
        "Intervention-based adaptive-law discrimination and the nonlinear "
        "matched-reservoir task benchmark."
    ),
    "p279_p280_conversion_torsors.png": (
        "Conditional conversion from relative information to free energy and "
        "the independent scale--orientation product torsor."
    ),
}


PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{
  unicode=true,
  colorlinks=true,
  linkcolor=finblue,
  citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Current State after Programs P267--P280},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Spectral identifiability, operational instruments, falsification boundaries, and the conditional passage from information to physics},
  pdfkeywords={FIN, spectral operator, Stieltjes measure, Schur complement, quantum instrument, current tomography, renormalization, Landauer principle, operational physics}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs P267--P280}
\fancyhead[R]{\small Release 10.25}
\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{8mm}
{\Large\bfseries FIN --- Release 10.25\par}
\vspace{10mm}
{\Huge\bfseries Current State of the Theory\par}
\vspace{4mm}
{\LARGE\bfseries after Research Programs P267--P280\par}
\vspace{7mm}
{\Large Spectral identifiability, operational instruments,\\
falsification boundaries, and the conditional passage\\
from information to physics\par}
\vspace{15mm}
{\Large Krzysztof \.Zuchowski\par}
\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
\vspace{2mm}
{\normalsize ORCID:
\href{https://orcid.org/0009-0002-0909-3613}{0009-0002-0909-3613}\par}
\vspace{8mm}
{\large 28 July 2026\par}
{\normalsize Publication --- Preprint; Version 1.0.0\par}
\vfill
\begin{minipage}{0.91\textwidth}
\small
\textbf{Central conclusion.}
FIN has a rigorous finite spectral-operator core, exact context memory, and
constructible operational interfaces. It does not yet internally generate a
physical unit, an orientation choice, a laboratory instrument, or a unique
physical interpretation.

\medskip
\textbf{Minimal surviving architecture.}
\[
(A,\mathcal C)
\quad\longrightarrow\quad
(A,\mathcal C,\mathfrak O)
\quad\longrightarrow\quad
(A,\mathcal C,\mathfrak O,\mathfrak C,\mathfrak S).
\]
The last arrow remains conditional on operational, conversion, and
symmetry-breaking resources.

\medskip
\textbf{Repository.}
\url{https://github.com/hyconiek/Fractal-Nadsoliton-Theory}

\textbf{License.} CC BY 4.0.
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}
\pagenumbering{roman}
\chapter*{Publication statement}
\addcontentsline{toc}{chapter}{Publication statement}
This monograph reports all executed Programs P267--P280, integrates their
results with the current FIN state map, states the mathematical and physical
boundaries explicitly, and recommends Programs P281--P294. Synthetic audits
are never described as laboratory evidence. The strict and legacy kernels
remain separated, QW-2191 remains open, and no physical-role or Theory of
Everything closure is claimed.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""


POSTAMBLE = r"""
\clearpage
\chapter*{Reproduction certificate}
\addcontentsline{toc}{chapter}{Reproduction certificate}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/matplotlib-p267-280 \
python3 fin_programs_267_280.py

MPLCONFIGDIR=/tmp/matplotlib-p267-280 \
python3 -m unittest -v test_fin_programs_267_280.py
\end{Verbatim}
Expected result: 18 tests passed, 0 failed. Programs P274, P275, P277, and P278
use synthetic data. P273 certifies that no independent production P241 event
bundle was available; P242 was not authorized as a physical confirmatory run.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Current State of the Theory after Research
Programs P267--P280: Spectral Identifiability, Operational Instruments,
Falsification Boundaries, and the Conditional Passage from Information to
Physics} (Release 10.25; Version 1.0.0) [Preprint]. Zenodo.

\end{document}
"""


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    text = re.sub(
        r"^(#{1,4})\s+\d+(?:\.\d+)*\.\s+",
        lambda match: match.group(1) + " ",
        text,
        flags=re.MULTILINE,
    )
    body = common_body(text, FIGURES, "FIN_Programs_267_280_Figures")
    body = body.replace(
        r"\section{Confidence convention}",
        r"\chapter*{Confidence convention}"
        "\n"
        r"\addcontentsline{toc}{chapter}{Confidence convention}",
    )
    TARGET.write_text(PREAMBLE + "\n" + body + "\n" + POSTAMBLE, encoding="utf-8")
    print(TARGET)


if __name__ == "__main__":
    main()
