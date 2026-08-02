#!/usr/bin/env python3
"""Build the archival TeX source for FIN checkpoint P471--P473."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body


SOURCE = Path("FIN_Local_Research_Checkpoint_P471_P473_EN.md")
TARGET = Path("FIN_Local_Research_Checkpoint_P471_P473_EN.tex")
FIGURES = {
    "p471_p473_polynomial_krawczyk.png": (
        "Polynomial Jacobian spectrum, exact Krawczyk inclusion margins, and "
        "the certified global optimum interval."
    ),
}

PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{unicode=true,colorlinks=true,linkcolor=finblue,citecolor=fingreen,
urlcolor=finblue,pdftitle={FIN Local Research Checkpoint P471--P473},
pdfauthor={Krzysztof Żuchowski},
pdfsubject={Exact O167 primal-dual attainment by polynomial Krawczyk certification},
pdfkeywords={FIN, quantum comb, Riccati equation, Krawczyk operator, interval arithmetic, semidefinite duality, exact attainment}}
\pagestyle{fancy}\fancyhf{}
\fancyhead[L]{\small FIN Checkpoint P471--P473}
\fancyhead[R]{\small Release 10.44}\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}\hypersetup{pageanchor=false}
\begin{titlepage}\centering\vspace*{8mm}
{\Large\bfseries FIN Local Research --- Release 10.44\par}\vspace{10mm}
{\Huge\bfseries Checkpoint P471--P473\par}\vspace{6mm}
{\Large Exact O167 Primal--Dual Attainment by Polynomial\
Krawczyk Certification\par}\vspace{14mm}
{\Large Krzysztof \.Zuchowski\par}\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
{\normalsize ORCID: 0009-0002-0909-3613\par}\vspace{8mm}
{\large 1 August 2026\par}{\normalsize Publication --- Preprint; Version 1.0.0\par}
\vfill
\begin{minipage}{0.92\textwidth}\small
\textbf{Central theorem.} A strict rational Krawczyk inclusion proves one
exact positive O167 primal--dual contact for the declared reduced three-slot
ordered comb. Its global value satisfies
\[
0.5233281002710117\le D_3^{\rm global}\le0.5233281002710717,
\]
an interval of width $6\times10^{-14}$.

\medskip\textbf{Method.} The complete support-ladder condition is reduced
symbolically to 13 cubic polynomial residual orbits, enclosed with exact
rational arithmetic, and supplemented by exact box-positivity certificates.

\medskip\textbf{Boundary.} One exact global optimizer is proved; uniqueness
over the complete causal cone remains open. No selector, dimensional source,
laboratory record, complete legacy--strict bridge, role transfer, or physical
closure is claimed.

\medskip\textbf{License.} CC BY 4.0.
\end{minipage}\vfill\end{titlepage}
\hypersetup{pageanchor=true}\pagenumbering{roman}
\chapter*{Publication statement}\addcontentsline{toc}{chapter}{Publication statement}
This monograph reports local analytical and computer-assisted Programs P471,
P472, and P473. It uses no laboratory record, external audit, internet result,
remote computation, or physical validation.
\tableofcontents\clearpage\pagenumbering{arabic}
"""

POSTAMBLE = r"""
\clearpage\chapter*{Reproduction certificate}
\addcontentsline{toc}{chapter}{Reproduction certificate}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/fin-mpl-1044 \
python3 fin_programs_471_472_473.py
python3 -m unittest -v \
  test_fin_programs_471_472_473.py \
  test_fin_checkpoint_p471_p473.py
python3 fin_programs_471_472_473_to_latex.py
lualatex -interaction=nonstopmode -halt-on-error \
  FIN_Local_Research_Checkpoint_P471_P473_EN.tex
sha256sum -c FIN_PROGRAMS_471_472_473_RELEASE_MANIFEST.sha256
\end{Verbatim}
\noindent\textbf{Suggested citation:}
\.Zuchowski, K. (2026). \emph{FIN Local Research Checkpoint P471--P473:
Exact O167 Primal--Dual Attainment by Polynomial Krawczyk Certification}
(Release 10.44; Version 1.0.0) [Preprint].
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
        "FIN_Programs_471_472_473_Figures",
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
