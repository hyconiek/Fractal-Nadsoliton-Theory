#!/usr/bin/env python3
"""Build the archival TeX source for FIN checkpoint P451--P453."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body


SOURCE = Path("FIN_Local_Research_Checkpoint_P451_P453_EN.md")
TARGET = Path("FIN_Local_Research_Checkpoint_P451_P453_EN.tex")
FIGURES = {
    "p451_p453_coherence_symmetry_and_gauge_fixing.png": (
        "The certified causal-coherence advantage, the global coarse-erasure "
        "palindromic reduction, and the unique minimum-TV signed measure."
    ),
}


PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{
  unicode=true,colorlinks=true,linkcolor=finblue,citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Local Research Checkpoint P451--P453},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Causal coherence, global coarse-erasure symmetry, and canonical Jordan representation},
  pdfkeywords={FIN, quantum comb, causal tester, coherence, channel discrimination, heralded erasure, symmetrization, signed moments, strict complementarity, total variation}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Checkpoint P451--P453}
\fancyhead[R]{\small Release 10.40}
\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{8mm}
{\Large\bfseries FIN Local Research --- Release 10.40\par}
\vspace{10mm}
{\Huge\bfseries Checkpoint P451--P453\par}
\vspace{6mm}
{\Large Coherence Advantage, Global Coarse-Erasure Symmetry,\
and Canonical Jordan Representation\par}
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
An exact-rational coherent three-slot causal tester has certified half-distance
at least 0.52332810026048937, exceeding the proved global diagonal optimum by
at least 0.00056473095851463685.

\medskip
\textbf{Additional theorems.}
The actual coarse three-use erasure objective is globally reduced to the
palindromic line, and the P429 signed-moment optimizer is globally unique
under the minimum-negative-mass rule.

\medskip
\textbf{Boundary.}
No full 21-dimensional optimum, selector, dimensional source, laboratory
record, complete legacy--strict bridge, role transfer, or physical closure
is claimed.

\medskip
\textbf{License.} CC BY 4.0.
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}
\pagenumbering{roman}
\chapter*{Publication statement}
\addcontentsline{toc}{chapter}{Publication statement}
This monograph reports local analytical and computer-assisted Programs P451,
P452, and P453. It uses no laboratory record, external audit, internet result,
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
python3 fin_programs_451_453.py

python3 -m unittest -v \
  test_fin_programs_451_453.py \
  test_fin_checkpoint_p451_p453.py

lualatex -interaction=nonstopmode -halt-on-error \
  FIN_Local_Research_Checkpoint_P451_P453_EN.tex

sha256sum -c FIN_PROGRAMS_451_453_RELEASE_MANIFEST.sha256
\end{Verbatim}

\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Local Research Checkpoint P451--P453:
Coherence Advantage, Global Coarse-Erasure Symmetry, and Canonical Jordan
Representation} (Release 10.40; Version 1.0.0) [Preprint].
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
        "FIN_Programs_451_453_Figures",
    )
    for plain, rendered in {
        "[Proven]": r"\textcolor{fingreen}{[Proven]}",
        "[Computer-assisted proof]": r"\textcolor{fingreen}{[Computer-assisted proof]}",
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
