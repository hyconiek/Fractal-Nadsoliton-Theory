#!/usr/bin/env python3
"""Build the archival TeX source for FIN checkpoint P448--P450."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body


SOURCE = Path("FIN_Local_Research_Checkpoint_P448_P450_EN.md")
TARGET = Path("FIN_Local_Research_Checkpoint_P448_P450_EN.tex")
FIGURES = {
    "p448_p450_global_majorant_and_gauge_obstruction.png": (
        "The global fine-erasure majorant, the three-slot causal-cone audit, "
        "and the exact moment-null representation-gauge obstruction."
    ),
}


PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{
  unicode=true,colorlinks=true,linkcolor=finblue,citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Local Research Checkpoint P448--P450},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={A full-simplex erasure majorant, a three-slot causal echo counterexample, and a representation-gauge no-go},
  pdfkeywords={FIN, quantum comb, causal tester, channel discrimination, echo history, heralded erasure, concavity, interval arithmetic, signed moments, representation gauge}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Checkpoint P448--P450}
\fancyhead[R]{\small Release 10.39}
\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{8mm}
{\Large\bfseries FIN Local Research --- Release 10.39\par}
\vspace{10mm}
{\Huge\bfseries Checkpoint P448--P450\par}
\vspace{6mm}
{\Large A Full-Simplex Erasure Majorant,\
a Three-Slot Causal Echo Counterexample,\
and a Representation-Gauge No-Go\par}
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
An explicit rational three-slot causal echo-history tester has certified
half-distance at least 0.49063899018433244, exceeding the three-slot GHZ
value by at least 0.017612669538553616.

\medskip
\textbf{Additional theorems.}
A concave fine-erasure majorant gives the first rigorous full-simplex P446
upper bound.  An exact twelve-moment null cycle proves that the O163 sampler
is not representation independent.

\medskip
\textbf{Boundary.}
No global three-slot optimum, selector, dimensional source, laboratory
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
This monograph reports local analytical and computer-assisted Programs P448,
P449, and P450.  It uses no laboratory record, external audit, internet
result, remote computation, or physical validation.
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
python3 fin_programs_448_450.py

python3 -m unittest -v \
  test_fin_programs_448_450.py \
  test_fin_checkpoint_p448_p450.py

lualatex -interaction=nonstopmode -halt-on-error \
  FIN_Local_Research_Checkpoint_P448_P450_EN.tex

sha256sum -c FIN_PROGRAMS_448_450_RELEASE_MANIFEST.sha256
\end{Verbatim}

\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Local Research Checkpoint P448--P450:
A Full-Simplex Erasure Majorant, a Three-Slot Causal Echo Counterexample,
and a Representation-Gauge No-Go} (Release 10.39; Version 1.0.0) [Preprint].
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
        "FIN_Programs_448_450_Figures",
    )
    for plain, rendered in {
        "[Computer-assisted proof]": r"\textcolor{fingreen}{[Computer-assisted proof]}",
        "[Conditional]": r"\textcolor{finviolet}{[Conditional]}",
        "[Open]": r"\textcolor{finorange}{[Open]}",
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
