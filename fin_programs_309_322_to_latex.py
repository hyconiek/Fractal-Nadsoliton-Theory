#!/usr/bin/env python3
"""Build the archival English PDF source for FIN Programs P309--P322."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body


SOURCE = Path("FIN_Programs_309_322_Bridge_Resource_Report_EN.md")
TARGET = Path("FIN_Programs_309_322_Bridge_Resource_Report_EN.tex")

FIGURES = {
    "p309_minimax_stieltjes.png": (
        "Finite-probe information improvement and the pole-collision minimax "
        "obstruction."
    ),
    "p311_lossy_mesh.png": (
        "Loss, phase drift, no-click completion, and detector response of the "
        "synthetic current readout."
    ),
    "p312_p313_parent_tail.png": (
        "A target-dependent joint-feature parent and the positive "
        "hyperbolic-mixture tail obstruction."
    ),
    "p315_p316_distance_negativity.png": (
        "Operational legacy/strict separation and signed path-resource "
        "bounds."
    ),
    "p318_p319_drift_chambers.png": (
        "Synthetic detector drift tracking and the spectral chamber complex."
    ),
    "p320_p322_clock_role.png": (
        "Restricted-clock design and carrier dependence of candidate physical "
        "role scalars."
    ),
}


PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{
  unicode=true,
  colorlinks=true,
  linkcolor=finblue,
  citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Research Programs P309--P322},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Bridge-resource lower bounds, operational separation, and sector-role obstruction},
  pdfkeywords={FIN, signed path resource, Hausdorff moment, nontorsion phase, common parent, operational distance, clock identifiability, sector pointing}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs P309--P322}
\fancyhead[R]{\small Release 10.28}
\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{7mm}
{\Large\bfseries FIN --- Release 10.28\par}
\vspace{9mm}
{\Huge\bfseries Research Programs P309--P322\par}
\vspace{6mm}
{\Large Bridge-resource lower bounds,\\
operational separation,\\
and the sector-role obstruction\par}
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
Strict attenuation lies outside the positive path-scale cone and requires
signed negative mass at least \(0.0715275581\) in the declared moment
representation.  Its frozen phase tuple additionally requires one
infinite-order generator not contained in the finite legacy phase group.

\medskip
\textbf{Bridge boundary.}
A cross-supported common parent exists after an explicit legacy positivity
repair, but it imports the strict operator.  The remaining bridge-resource
budget keeps signed path weight, nontorsion phase, cross law, pointing, scale,
and operational records separately typed.

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
This monograph reports all executed Programs P309--P322 and recommends
Programs P323--P336. Mathematical proof, target-dependent construction,
synthetic evidence, scoped refutation, and missing external evidence remain
separate. The legacy/strict split, QW-2191, dimensional-source, physical-role,
laboratory, \(L_{\rm total}\), and Theory-of-Everything boundaries are
preserved.

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
python3 -m unittest -v test_fin_programs_309_322.py
\end{Verbatim}
Expected result: 17 tests passed, 0 failed.  The suite recompiles the
dependency-free Lean 4.28 source.  Loss, detector drift, calibration, and
finite-count assumptions are synthetic.  P317 and P321 admit no external
laboratory or hardware evidence.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Research Programs P309--P322:
Bridge-Resource Lower Bounds, Operational Separation, and the Sector-Role
Obstruction} (Release 10.28; Version 1.0.0) [Preprint]. Zenodo.

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
        "FIN_Programs_309_322_Figures",
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
