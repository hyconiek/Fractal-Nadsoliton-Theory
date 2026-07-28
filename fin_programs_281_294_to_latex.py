#!/usr/bin/env python3
"""Build the archival English PDF source for FIN Programs P281--P294."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body


SOURCE = Path("FIN_Programs_281_294_Research_Report_EN.md")
TARGET = Path("FIN_Programs_281_294_Research_Report_EN.tex")

FIGURES = {
    "p281_p283_recovery_povm.png": (
        "Bounded Stieltjes recovery and the current POVM under a calibrated "
        "detector channel."
    ),
    "p284_p286_rg_bridge.png": (
        "Complement RG rank restoration and the failed one-atom fractal "
        "legacy completion."
    ),
    "p288_p289_detector_false_positive.png": (
        "Detector-aware flux design and the declared-null false-positive "
        "audit."
    ),
    "p290_p291_mechanism_intervention.png": (
        "Three-time mechanism discrimination and maximin intervention design."
    ),
    "p292_p294_replication_thermo_torsor.png": (
        "Reservoir replication, the thermodynamic resource ledger, and "
        "pointed-torsor counts."
    ),
}


PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{
  unicode=true,
  colorlinks=true,
  linkcolor=finblue,
  citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Research Programs P281--P294},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Regularized spectral inference, resource-explicit physics, and the pointed-torsor boundary},
  pdfkeywords={FIN, spectral operator, Stieltjes inverse, Schur complement, Naimark dilation, process ledger, mechanism identification, Landauer erasure, torsor}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs P281--P294}
\fancyhead[R]{\small Release 10.26}
\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{8mm}
{\Large\bfseries FIN --- Release 10.26\par}
\vspace{10mm}
{\Huge\bfseries Research Programs P281--P294\par}
\vspace{7mm}
{\Large Regularized spectral inference, resource-explicit physics,\\
and the pointed-torsor boundary\par}
\vspace{15mm}
{\Large Krzysztof \.Zuchowski\par}
\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
\vspace{2mm}
{\normalsize ORCID:
\href{https://orcid.org/0009-0002-0909-3613}{0009-0002-0909-3613}\par}
\vspace{8mm}
{\large 29 July 2026\par}
{\normalsize Publication --- Preprint; Version 1.0.0\par}
\vfill
\begin{minipage}{0.91\textwidth}
\small
\textbf{Central result.}
The strict finite spectral generator supports exact mathematical and
operational constructions, but physicalization still requires explicitly
supplied operational, dimensional, sector, and pointing resources.

\medskip
\textbf{New boundary object.}
\[
\mathfrak P=(T,t_0),
\]
the Pointed Physicalization Resource, removes an operational torsor ambiguity
after a point is supplied. Equivariance alone cannot derive that point.

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
This monograph reports all executed Programs P281--P294 and recommends
Programs P295--P308. Proof, simulation, conditional construction, refutation,
and missing empirical evidence are kept distinct. The strict and legacy
kernels remain separated, QW-2191 remains open, and no role transfer,
dimensional emergence, physical closure, or Theory of Everything claim is
made.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""


POSTAMBLE = r"""
\clearpage
\chapter*{Reproduction certificate}
\addcontentsline{toc}{chapter}{Reproduction certificate}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/matplotlib-p281-294 \
python3 -m unittest -v test_fin_programs_281_294.py
\end{Verbatim}
Expected result: 19 tests passed, 0 failed. The suite recompiles the
dependency-free Lean 4.28 rational Schur witness. Synthetic and conditional
studies are not laboratory evidence. P287 certifies that no independent
production P241 event bundle was available and therefore P242 was not
authorized.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Research Programs P281--P294:
Regularized Spectral Inference, Resource-Explicit Physics, and the
Pointed-Torsor Boundary} (Release 10.26; Version 1.0.0) [Preprint]. Zenodo.

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
    body = common_body(text, FIGURES, "FIN_Programs_281_294_Figures")
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
