#!/usr/bin/env python3
"""Build the archival TeX source for FIN checkpoint P465--P469."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body


SOURCE = Path("FIN_Local_Research_Checkpoint_P465_P469_EN.md")
TARGET = Path("FIN_Local_Research_Checkpoint_P465_P469_EN.tex")
FIGURES = {
    "p465_p468_p469_certificates.png": (
        "Strict symmetrization equality, rational dual-gap collapse, and the "
        "floating full-ladder KKT residuals."
    ),
}

PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{unicode=true,colorlinks=true,linkcolor=finblue,citecolor=fingreen,
urlcolor=finblue,pdftitle={FIN Local Research Checkpoint P465--P469},
pdfauthor={Krzysztof Żuchowski},
pdfsubject={Strict full-simplex uniqueness and a support-ladder rational dual certificate},
pdfkeywords={FIN, strict symmetrization, equality cases, quantum comb, semidefinite duality, KKT, rational certificate}}
\pagestyle{fancy}\fancyhf{}
\fancyhead[L]{\small FIN Checkpoint P465--P469}
\fancyhead[R]{\small Release 10.43}\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}\hypersetup{pageanchor=false}
\begin{titlepage}\centering\vspace*{8mm}
{\Large\bfseries FIN Local Research --- Release 10.43\par}\vspace{10mm}
{\Huge\bfseries Checkpoint P465--P469\par}\vspace{6mm}
{\Large Strict Full-Simplex Uniqueness and a Support-Ladder\
Rational Dual Certificate\par}\vspace{14mm}
{\Large Krzysztof \.Zuchowski\par}\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
{\normalsize ORCID: 0009-0002-0909-3613\par}\vspace{8mm}
{\large 1 August 2026\par}{\normalsize Publication --- Preprint; Version 1.0.0\par}
\vfill
\begin{minipage}{0.92\textwidth}\small
\textbf{Central results.} Reversal symmetrization is proved strict outside
the palindromic fixed set, closing full-simplex uniqueness for the declared
coarse-erasure model. A new exact-rational nested dual proves
\[
0.52332810026048937\le D_3^{\rm global}\le0.52332810067104088,
\]
with exact gap $4.1055151\times10^{-10}$.

\medskip\textbf{Localized frontier.} The canonical O167 support ladder
satisfies the complete floating KKT residual to $1.33\times10^{-15}$, but an
exact interval root and exact attainment remain open.

\medskip\textbf{Boundary.} No selector, dimensional source, laboratory
record, complete legacy--strict bridge, role transfer, full physical closure,
or Theory of Everything is claimed.

\medskip\textbf{License.} CC BY 4.0.
\end{minipage}\vfill\end{titlepage}
\hypersetup{pageanchor=true}\pagenumbering{roman}
\chapter*{Publication statement}\addcontentsline{toc}{chapter}{Publication statement}
This monograph reports local analytical and computer-assisted Programs P465,
P468, and P469. It uses no laboratory record, external audit, internet result,
remote computation, or physical validation.
\tableofcontents\clearpage\pagenumbering{arabic}
"""

POSTAMBLE = r"""
\clearpage\chapter*{Reproduction certificate}
\addcontentsline{toc}{chapter}{Reproduction certificate}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/fin-mpl-1043 \
python3 fin_programs_465_468_469.py
python3 -m unittest -v \
  test_fin_programs_465_468_469.py \
  test_fin_checkpoint_p465_p469.py
python3 fin_programs_465_468_469_to_latex.py
lualatex -interaction=nonstopmode -halt-on-error \
  FIN_Local_Research_Checkpoint_P465_P469_EN.tex
sha256sum -c FIN_PROGRAMS_465_468_469_RELEASE_MANIFEST.sha256
\end{Verbatim}
\noindent\textbf{Suggested citation:}
\.Zuchowski, K. (2026). \emph{FIN Local Research Checkpoint P465--P469:
Strict Full-Simplex Uniqueness and a Support-Ladder Rational Dual Certificate}
(Release 10.43; Version 1.0.0) [Preprint].
\end{document}
"""


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    text = re.sub(r"\\\((.*?)\\\)", r"$\1$", text)
    text = re.sub(
        r"^(#{1,4})\s+\d+(?:\.\d+)*(?:\.)?\s+",
        lambda match: match.group(1) + " ",
        text,
        flags=re.MULTILINE,
    )
    body = common_body(
        "## Confidence convention\n\n" + text,
        FIGURES,
        "FIN_Programs_465_468_469_Figures",
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
    start = body.find(r"\chapter{Executive summary}")
    if start < 0:
        raise RuntimeError("could not locate executive summary")
    TARGET.write_text(PREAMBLE + "\n" + body[start:] + "\n" + POSTAMBLE, encoding="utf-8")
    print(TARGET)


if __name__ == "__main__":
    main()
