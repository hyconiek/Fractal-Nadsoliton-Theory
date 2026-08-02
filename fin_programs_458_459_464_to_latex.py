#!/usr/bin/env python3
"""Build the archival TeX source for FIN checkpoint P458--P464."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body

SOURCE = Path("FIN_Local_Research_Checkpoint_P458_P464_EN.md")
TARGET = Path("FIN_Local_Research_Checkpoint_P458_P464_EN.tex")
FIGURES = {
    "p458_p459_p464_certificates.png": (
        "Certified palindromic curvature, fiberwise detector allocation, and "
        "the admitted nested-dual refinement."
    ),
}

PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{unicode=true,colorlinks=true,linkcolor=finblue,citecolor=fingreen,
urlcolor=finblue,pdftitle={FIN Local Research Checkpoint P458--P464},
pdfauthor={Krzysztof Żuchowski},
pdfsubject={Unique palindromic curvature, canonical allocation, and nested-dual refinement},
pdfkeywords={FIN, interval automatic differentiation, strict concavity, detector allocation, quantum comb, semidefinite duality, rational certificate}}
\pagestyle{fancy}\fancyhf{}
\fancyhead[L]{\small FIN Checkpoint P458--P464}
\fancyhead[R]{\small Release 10.42}\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}\hypersetup{pageanchor=false}
\begin{titlepage}\centering\vspace*{8mm}
{\Large\bfseries FIN Local Research --- Release 10.42\par}\vspace{10mm}
{\Huge\bfseries Checkpoint P458--P464\par}\vspace{6mm}
{\Large Unique Palindromic Curvature, Fiberwise Canonical Allocation,\
and a Thirty-Fold Dual-Gap Reduction\par}\vspace{14mm}
{\Large Krzysztof \.Zuchowski\par}\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
{\normalsize ORCID: 0009-0002-0909-3613\par}\vspace{8mm}
{\large 1 August 2026\par}{\normalsize Publication --- Preprint; Version 1.0.0\par}
\vfill
\begin{minipage}{0.92\textwidth}\small
\textbf{Central result.} A new exact-rational nested-dual certificate proves
\[0.52332810026048937\le D_3^{\rm global}\le0.52332832027103,\]
reducing the certified gap to $2.2001054063\times10^{-7}$.

\medskip\textbf{Additional theorems.} The declared palindromic coarse-erasure
maximizer is unique and isolated to width below $1.22\times10^{-13}$. Under
the explicit minimum-TV rule, every supplied detector-box fiber has one unique
strictly convex minimax allocation.

\medskip\textbf{Boundary.} No exact full-cone optimizer, full-simplex
uniqueness, selector, dimensional source, laboratory record, complete
legacy--strict bridge, role transfer, or physical closure is claimed.

\medskip\textbf{License.} CC BY 4.0.
\end{minipage}\vfill\end{titlepage}
\hypersetup{pageanchor=true}\pagenumbering{roman}
\chapter*{Publication statement}\addcontentsline{toc}{chapter}{Publication statement}
This monograph reports local analytical and computer-assisted Programs P458,
P459, and P464. It uses no laboratory record, external audit, internet result,
remote computation, or physical validation.
\tableofcontents\clearpage\pagenumbering{arabic}
"""

POSTAMBLE = r"""
\clearpage\chapter*{Reproduction certificate}
\addcontentsline{toc}{chapter}{Reproduction certificate}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/fin-mpl-1038 \
python3 fin_programs_458_459_464.py
python3 -m unittest -v \
  test_fin_programs_458_459_464.py \
  test_fin_checkpoint_p458_p464.py
python3 fin_programs_458_459_464_to_latex.py
lualatex -interaction=nonstopmode -halt-on-error \
  FIN_Local_Research_Checkpoint_P458_P464_EN.tex
sha256sum -c FIN_PROGRAMS_458_459_464_RELEASE_MANIFEST.sha256
\end{Verbatim}
\noindent\textbf{Suggested citation:}
\.Zuchowski, K. (2026). \emph{FIN Local Research Checkpoint P458--P464:
Unique Palindromic Curvature, Fiberwise Canonical Allocation, and a Thirty-Fold
Dual-Gap Reduction} (Release 10.42; Version 1.0.0) [Preprint].
\end{document}
"""


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    text = re.sub(r"\\\((.*?)\\\)", r"$\1$", text)
    text = re.sub(r"^(#{1,4})\s+\d+(?:\.\d+)*(?:\.)?\s+", lambda m: m.group(1)+" ", text, flags=re.MULTILINE)
    body = common_body("## Confidence convention\n\n"+text, FIGURES, "FIN_Programs_458_459_464_Figures")
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
    start = body.find(r"\chapter{Executive summary}")
    if start < 0:
        raise RuntimeError("could not locate executive summary")
    TARGET.write_text(PREAMBLE+"\n"+body[start:]+"\n"+POSTAMBLE, encoding="utf-8")
    print(TARGET)


if __name__ == "__main__":
    main()
