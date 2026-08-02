#!/usr/bin/env python3
"""Build the archival TeX source for FIN checkpoint P454--P457."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body


SOURCE = Path("FIN_Local_Research_Checkpoint_P454_P457_EN.md")
TARGET = Path("FIN_Local_Research_Checkpoint_P454_P457_EN.tex")
FIGURES = {
    "p454_p457_dual_symmetry_and_refined_cover.png": (
        "The nested-dual central path and certified global bracket, the "
        "five-dimensional ordered-symmetry audit, and the refined globally "
        "licensed coarse-erasure cover."
    ),
}


PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{
  unicode=true,colorlinks=true,linkcolor=finblue,citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Local Research Checkpoint P454--P457},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Certified nested comb dual, ordered-symmetry obstruction, and six-decimal global cover},
  pdfkeywords={FIN, quantum comb, semidefinite duality, rational certificate, causal tester, symmetry obstruction, heralded erasure, interval cover}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Checkpoint P454--P457}
\fancyhead[R]{\small Release 10.41}
\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{8mm}
{\Large\bfseries FIN Local Research --- Release 10.41\par}
\vspace{10mm}
{\Huge\bfseries Checkpoint P454--P457\par}
\vspace{6mm}
{\Large A Certified Nested Comb Dual, an Ordered-Symmetry Obstruction,\
and a Six-Decimal Global Cover\par}
\vspace{14mm}
{\Large Krzysztof \.Zuchowski\par}
\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
\vspace{2mm}
{\normalsize ORCID: 0009-0002-0909-3613\par}
\vspace{8mm}
{\large 1 August 2026\par}
{\normalsize Publication --- Preprint; Version 1.0.0\par}
\vfill
\begin{minipage}{0.92\textwidth}
\small
\textbf{Central result.}
An explicit rational nested-dual certificate and the inherited coherent
primal witness confine the global three-slot causal half-distance to
\[
0.52332810026048937\le D_3^{\rm global}\le0.523334700252,
\]
a rigorous gap below $6.6\times10^{-6}$.

\medskip
\textbf{Additional results.}
Known ordered-comb symmetries leave a five-dimensional fixed space and cannot
force the three-dimensional coherent ansatz. A new globally licensed interval
cover confines the coarse erasure optimum to a gap below $10^{-6}$.

\medskip
\textbf{Boundary.}
No exact optimizer, optimizer uniqueness, selector, dimensional source,
laboratory record, complete legacy--strict bridge, role transfer, or physical
closure is claimed.

\medskip
\textbf{License.} CC BY 4.0.
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}
\pagenumbering{roman}
\chapter*{Publication statement}
\addcontentsline{toc}{chapter}{Publication statement}
This monograph reports local analytical and computer-assisted Programs P454,
P455, and P457. It uses no laboratory record, external audit, internet result,
remote computation, or physical validation.
\tableofcontents
\clearpage
\pagenumbering{arabic}
"""


POSTAMBLE = r"""
\clearpage
\chapter*{Reproduction certificate}
\addcontentsline{toc}{chapter}{Reproduction certificate}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/fin-mpl-1038 \
python3 fin_programs_454_455_457.py

MPLCONFIGDIR=/tmp/fin-mpl-1038 \
python3 -m unittest -v \
  test_fin_programs_454_455_457.py \
  test_fin_checkpoint_p454_p457.py

python3 fin_programs_454_455_457_to_latex.py
lualatex -interaction=nonstopmode -halt-on-error \
  FIN_Local_Research_Checkpoint_P454_P457_EN.tex

sha256sum -c FIN_PROGRAMS_454_455_457_RELEASE_MANIFEST.sha256
\end{Verbatim}

\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Local Research Checkpoint P454--P457:
A Certified Nested Comb Dual, an Ordered-Symmetry Obstruction, and a
Six-Decimal Global Cover} (Release 10.41; Version 1.0.0) [Preprint].
\end{document}
"""


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    text = re.sub(r"\\\((.*?)\\\)", r"$\1$", text)
    text = re.sub(
        r"^(#{1,4})\s+\d+(?:\.\d+)*(?:\.)?\s+",
        lambda match: match.group(1) + " ", text, flags=re.MULTILINE,
    )
    body = common_body(
        "## Confidence convention\n\n" + text,
        FIGURES,
        "FIN_Programs_454_455_457_Figures",
    )
    for plain, rendered in {
        "[Proven]": r"\textcolor{fingreen}{[Proven]}",
        "[Computer-assisted proof]": r"\textcolor{fingreen}{[Computer-assisted proof]}",
        "[Strong evidence]": r"\textcolor{finviolet}{[Strong evidence]}",
        "[Conditional]": r"\textcolor{finviolet}{[Conditional]}",
        "[Open]": r"\textcolor{finorange}{[Open]}",
        "[Refuted]": r"\textcolor{finred}{[Refuted]}",
        "[Blocked by external evidence]": r"\textcolor{finred}{[Blocked by external evidence]}",
    }.items():
        body = body.replace(plain, rendered)
    first_chapter = body.find(r"\chapter{Executive summary}")
    if first_chapter < 0:
        raise RuntimeError("could not locate executive summary")
    body = body[first_chapter:]
    TARGET.write_text(PREAMBLE + "\n" + body + "\n" + POSTAMBLE, encoding="utf-8")
    print(TARGET)


if __name__ == "__main__":
    main()
