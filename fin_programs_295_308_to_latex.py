#!/usr/bin/env python3
"""Build the archival English PDF source for FIN Programs P295--P308."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body


SOURCE = Path("FIN_Programs_295_308_Legacy_Strict_Physics_Report_EN.md")
TARGET = Path("FIN_Programs_295_308_Legacy_Strict_Physics_Report_EN.tex")

FIGURES = {
    "p295_p298_inverse_completion.png": (
        "Multi-probe inverse recovery and the scalar/Bernstein completion "
        "obstruction."
    ),
    "p297_p299_mesh_memory.png": (
        "Exact complex-unitary current-mesh compilation and a multitime "
        "memory witness."
    ),
    "p302_p303_detector_rare_event.png": (
        "Joint detector nuisance design and strict-fingerprint rare-event "
        "geometry."
    ),
    "p304_p305_clock_adaptation.png": (
        "Finite-time clock ambiguity and intervention-assisted synthetic "
        "completion-coordinate recovery."
    ),
    "p300_p307_kernel_role_invariance.png": (
        "The legacy/strict kernel split and the coordinate-gauge failure of "
        "historical physical-role formulas."
    ),
}


PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{
  unicode=true,
  colorlinks=true,
  linkcolor=finblue,
  citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Research Programs P295--P308},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Legacy physics, strict spectral completion, and physical-role transfer},
  pdfkeywords={FIN, legacy kernel, strict kernel, positive Laplace mixture, Bernstein function, Naimark dilation, scale gauge, role transfer}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs P295--P308}
\fancyhead[R]{\small Release 10.27}
\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{7mm}
{\Large\bfseries FIN --- Release 10.27\par}
\vspace{9mm}
{\Huge\bfseries Research Programs P295--P308\par}
\vspace{6mm}
{\Large Legacy physics revisited:\\
positive path mixtures, strict spectral completion,\\
and physical-role transfer\par}
\vspace{13mm}
{\Large Krzysztof \.Zuchowski\par}
\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
\vspace{2mm}
{\normalsize ORCID:
\href{https://orcid.org/0009-0002-0909-3613}{0009-0002-0909-3613}\par}
\vspace{7mm}
{\large 30 July 2026\par}
{\normalsize Publication --- Preprint; Version 1.0.0\par}
\vfill
\begin{minipage}{0.92\textwidth}
\small
\textbf{Central result.}
The legacy hyperbolic envelope is an exact positive mixture of exponential
distance attenuations.  The strict envelope, phase signs, and spectral order
lie outside positive fixed-phase and monotone scalar-completion classes.

\medskip
\textbf{New boundary objects.}
The Typed Legacy--Strict Completion Span and the Role-Transfer Certificate
separate mathematical completion from physical interpretation.  No historical
electroweak, electromagnetic, or gravity role currently passes the
certificate.

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
This monograph reports all executed Programs P295--P308 and recommends
Programs P309--P322.  It corrects historical legacy-kernel calculations,
preserves the legacy/strict split, and keeps mathematical proof, synthetic
evidence, scoped refutation, and missing external evidence distinct.
QW-2191, full bridge completion, physical-role transfer, dimensional
emergence, laboratory validation, \(L_{\rm total}\), and Theory-of-Everything
closure remain open or blocked.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""


POSTAMBLE = r"""
\clearpage
\chapter*{Reproduction certificate}
\addcontentsline{toc}{chapter}{Reproduction certificate}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/mpl-fin \
python3 -m unittest -v test_fin_programs_295_308.py
\end{Verbatim}
Expected result: 17 tests passed, 0 failed.  The suite recompiles the
dependency-free Lean 4.28 certificate.  Detector, memory, adaptation, and
finite-count studies are synthetic.  No independent P241 laboratory bundle
or certified hardware reservoir record was admitted.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Research Programs P295--P308:
Legacy Physics Revisited---Positive Path Mixtures, Strict Spectral
Completion, and Physical-Role Transfer} (Release 10.27; Version 1.0.0)
[Preprint]. Zenodo.

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
    body = common_body(
        text,
        FIGURES,
        "FIN_Programs_295_308_Figures",
    )
    body = body.replace(
        r"\section{Confidence convention}",
        r"\chapter*{Confidence convention}"
        "\n"
        r"\addcontentsline{toc}{chapter}{Confidence convention}",
    )
    TARGET.write_text(
        PREAMBLE + "\n" + body + "\n" + POSTAMBLE,
        encoding="utf-8",
    )
    print(TARGET)


if __name__ == "__main__":
    main()
