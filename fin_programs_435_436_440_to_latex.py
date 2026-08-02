#!/usr/bin/env python3
"""Build the archival TeX source for FIN checkpoint P435--P440."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body


SOURCE = Path("FIN_Local_Research_Checkpoint_P435_P440_EN.md")
TARGET = Path("FIN_Local_Research_Checkpoint_P435_P440_EN.tex")
FIGURES = {
    "p435_p436_p440_certificates.png": (
        "Exact one-slot process-tester equality, the rigorously separated "
        "P436 code intervals, and the detector-envelope sampling allocation."
    ),
}


PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{
  unicode=true,
  colorlinks=true,
  linkcolor=finblue,
  citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Local Research Checkpoint P435--P440},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Exact process-tester admission, certified heralded-code gain, and detector-envelope minimax Jordan sampling},
  pdfkeywords={FIN, quantum comb, process tester, SDP, heralded erasure, trace norm, interval arithmetic, Jordan sampling, detector efficiency, dark counts}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Checkpoint P435--P440}
\fancyhead[R]{\small Release 10.37}
\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{8mm}
{\Large\bfseries FIN Local Research --- Release 10.37\par}
\vspace{10mm}
{\Huge\bfseries Checkpoint P435--P440\par}
\vspace{6mm}
{\Large Exact Process-Tester Admission, a Certified Heralded-Code Gain,\
and Detector-Envelope Minimax Jordan Sampling\par}
\vspace{15mm}
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
\textbf{Strongest theorem.}
A rational three-use symmetric code has certified heralded phase/dephasing
gain at least
\[
0.022572776021405654
\]
over both product and GHZ baselines in the declared nonideal cell.

\medskip
\textbf{Operational boundary.}
The one-slot process-tester SDP is exact.  The genuine two-slot comb is
exported but has no matching local dual certificate.  Detector response is a
supplied mathematical envelope, not apparatus calibration or laboratory data.

\medskip
No selector, dimensional source, complete legacy--strict bridge, role
transfer, laboratory validation, Standard Model, gravity, or
Theory-of-Everything closure is claimed.

\medskip
\textbf{License.} CC BY 4.0.
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}
\pagenumbering{roman}
\chapter*{Publication statement}
\addcontentsline{toc}{chapter}{Publication statement}
This monograph reports local analytical and computer-assisted Programs P435,
P436, and P440.  It uses no laboratory data, external audit, internet
resource, remote computation, or physical validation.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""


POSTAMBLE = r"""
\clearpage
\chapter*{Reproduction certificate}
\addcontentsline{toc}{chapter}{Reproduction certificate}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/fin-mpl-1037 \
python3 fin_programs_435_436_440.py

python3 -m unittest -v \
  test_fin_programs_435_436_440.py \
  test_fin_checkpoint_p435_p440.py

lualatex -interaction=nonstopmode -halt-on-error \
  FIN_Local_Research_Checkpoint_P435_P440_EN.tex

sha256sum -c \
  FIN_PROGRAMS_435_436_440_RELEASE_MANIFEST.sha256
\end{Verbatim}

\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Local Research Checkpoint P435--P440:
Exact Process-Tester Admission, a Certified Heralded-Code Gain, and
Detector-Envelope Minimax Jordan Sampling} (Release 10.37; Version 1.0.0)
[Preprint].

\end{document}
"""


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    text = re.sub(
        r"^(#{1,4})\s+\d+(?:\.\d+)*(?:\.)?\s+",
        lambda match: match.group(1) + " ",
        text,
        flags=re.MULTILINE,
    )
    # The shared archival converter normally starts at a historical
    # "Confidence convention" marker.  Inject a disposable marker so this
    # checkpoint can retain its own title and section vocabulary; the slice
    # below still removes all front matter before the executive summary.
    body = common_body(
        "## Confidence convention\n\n" + text,
        FIGURES,
        "FIN_Programs_435_436_440_Figures",
    )
    for plain, rendered in {
        "[Computer-assisted proof]": r"\textcolor{fingreen}{[Computer-assisted proof]}",
        "[Constructed]": r"\textcolor{finblue}{[Constructed]}",
        "[Conditional]": r"\textcolor{finviolet}{[Conditional]}",
        "[Blocked locally]": r"\textcolor{finorange}{[Blocked locally]}",
        "[Blocked by external evidence]": r"\textcolor{finred}{[Blocked by external evidence]}",
    }.items():
        body = body.replace(plain, rendered)
    # The Markdown repeats bibliographic metadata for machine-readable use;
    # suppress that front matter because the archival title page contains it.
    first_chapter = body.find(r"\chapter{Executive summary}")
    if first_chapter < 0:
        raise RuntimeError("could not locate the first report chapter")
    body = body[first_chapter:]
    TARGET.write_text(PREAMBLE + "\n" + body + "\n" + POSTAMBLE, encoding="utf-8")
    print(TARGET)


if __name__ == "__main__":
    main()
