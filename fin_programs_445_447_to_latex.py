#!/usr/bin/env python3
"""Build the archival TeX source for FIN checkpoint P445--P447."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body


SOURCE = Path("FIN_Local_Research_Checkpoint_P445_P447_EN.md")
TARGET = Path("FIN_Local_Research_Checkpoint_P445_P447_EN.tex")
FIGURES = {
    "p445_p447_exact_reduction_and_intervals.png": (
        "The exact P445 two-slot envelope, the scoped P446 palindromic "
        "certificate, and the full P429-box detector-allocation tube."
    ),
}


PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{
  unicode=true,colorlinks=true,linkcolor=finblue,citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Local Research Checkpoint P445--P447},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Exact two-slot noisy-comb closure, scoped global heralded-code certification, and detector-allocation interval propagation},
  pdfkeywords={FIN, quantum comb, process tester, dephasing, GHZ, heralded erasure, representation theory, interval branch and bound, Jordan sampling}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Checkpoint P445--P447}
\fancyhead[R]{\small Release 10.38}
\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{8mm}
{\Large\bfseries FIN Local Research --- Release 10.38\par}
\vspace{10mm}
{\Huge\bfseries Checkpoint P445--P447\par}
\vspace{6mm}
{\Large Exact Closure of a Two-Slot Noisy Comb,\
a Scoped Global Code Certificate,\
and Full Atom-Box Propagation\par}
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
\textbf{Central theorem.}
For every causal two-slot tester of the declared reduced phase/dephasing
channel at q=4/5 and theta=pi/8, the exact optimal half-distance is
eight times the square root of two divided by twenty-five.
The GHZ history law is globally optimal, including arbitrary intermediate
memory.

\medskip
\textbf{Boundary.}
The three-use full code simplex is not closed: only a global \(10^{-3}\)
palindromic certificate and strong full-simplex evidence are obtained.  No
laboratory data, detector calibration, selector, unit source, complete
legacy--strict bridge, role transfer, or physical closure is claimed.

\medskip
\textbf{License.} CC BY 4.0.
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}
\pagenumbering{roman}
\chapter*{Publication statement}
\addcontentsline{toc}{chapter}{Publication statement}
This monograph reports local analytical and computer-assisted Programs P445,
P446, and P447.  It uses no laboratory record, external audit, internet
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
python3 fin_programs_445_447.py

python3 -m unittest -v \
  test_fin_programs_445_447.py \
  test_fin_checkpoint_p445_p447.py

lualatex -interaction=nonstopmode -halt-on-error \
  FIN_Local_Research_Checkpoint_P445_P447_EN.tex

sha256sum -c FIN_PROGRAMS_445_447_RELEASE_MANIFEST.sha256
\end{Verbatim}

\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Local Research Checkpoint P445--P447:
Exact Closure of a Two-Slot Noisy Comb, a Scoped Global Code Certificate,
and Full Atom-Box Propagation} (Release 10.38; Version 1.0.0) [Preprint].
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
        "FIN_Programs_445_447_Figures",
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
