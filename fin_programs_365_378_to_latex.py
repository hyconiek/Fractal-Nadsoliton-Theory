#!/usr/bin/env python3
"""Build the archival English PDF source for FIN Programs P365--P378."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body


SOURCE = Path(
    "FIN_Programs_365_378_Adaptive_Operational_Universality_Report_EN.md"
)
TARGET = Path(
    "FIN_Programs_365_378_Adaptive_Operational_Universality_Report_EN.tex"
)

FIGURES = {
    "p365_p366_certificate_tightening.png": (
        "Independent exact certificate checking and the tightened "
        "oscillatory signed-moment enclosure."
    ),
    "p368_p369_naturality_identifiability.png": (
        "The maximal radial-kernel naturality category and component-level "
        "photonic calibration identifiability."
    ),
    "p369_p371_loss_comb.png": (
        "Synthetic component-loss sensitivity and the exact ideal adaptive "
        "discrimination law through the first perfect time."
    ),
    "p373_p374_semantics_outer.png": (
        "Operational data-processing semantics and the noncanonical "
        "zero-free outer analytic extension domain."
    ),
}


PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{
  unicode=true,
  colorlinks=true,
  linkcolor=finblue,
  citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Research Programs P365--P378},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Exact certificate diversity, adaptive discrimination, operational resource semantics, and the universality boundary},
  pdfkeywords={FIN, exact certificate, Krawczyk method, Bernstein certificate, radial kernel naturality, adaptive quantum discrimination, data processing, outer function, operational resource theory}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs P365--P378}
\fancyhead[R]{\small Release 10.32}
\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{7mm}
{\Large\bfseries FIN --- Release 10.32\par}
\vspace{9mm}
{\Huge\bfseries Research Programs P365--P378\par}
\vspace{6mm}
{\Large Exact certificate diversity, adaptive discrimination,\\
operational resource semantics, and the universality boundary\par}
\vspace{13mm}
{\Large Krzysztof \.Zuchowski\par}
\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
\vspace{2mm}
{\normalsize ORCID:
\href{https://orcid.org/0009-0002-0909-3613}{0009-0002-0909-3613}\par}
\vspace{7mm}
{\large 31 July 2026\par}
{\normalsize Publication --- Preprint; Version 1.0.0\par}
\vfill
\begin{minipage}{0.92\textwidth}
\small
\textbf{Tightened oscillatory resource.}
\[
0.7072784585518420
\leq\mathcal N_{\rm osc}
\leq0.7073534683998260.
\]

\textbf{Ideal adaptive discrimination, first perfect branch.}
\[
\frac12\left\|
\mathcal U_{\rm strict}^{\otimes n}
-\mathcal U_{\rm legacy}^{\otimes n}
\right\|_{\diamond}
=
\sin\!\left(\frac{nt\Delta}{2}\right),
\qquad
0\leq t\leq\frac{\pi}{n\Delta}.
\]

The signed quantity is a classical representation cost. The channel theorem
is restricted to the declared commuting, ideal-unitary comparison and its
first perfect time. No selector, dimensional standard, laboratory
realization, physical role transfer, Standard Model, gravity, or
Theory-of-Everything closure is claimed.

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
This monograph reports every executed Program P365--P378 and recommends
Programs P379--P393. Exact proof, high-precision computation, synthetic
apparatus modelling, conditional operational semantics, and absent external
evidence are kept distinct. The legacy/strict split and all current
strict-core guardrails remain binding.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""


POSTAMBLE = r"""
\clearpage
\chapter*{Reproduction certificate}
\addcontentsline{toc}{chapter}{Reproduction certificate}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/mpl-fin-1032 \
python3 fin_programs_365_378.py

MPLCONFIGDIR=/tmp/mpl-fin-1032 \
python3 -m unittest -v test_fin_programs_365_378.py
\end{Verbatim}
Expected result: 14 tests passed, 0 failed. The dependency-free Lean 4.28
structural core must compile. P367, P370, P372, P375, P376, and P378 remain
external admission gates; no synthetic record is accepted as laboratory,
standards, hold-out, or custody evidence.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Research Programs P365--P378:
Exact Certificate Diversity, Adaptive Discrimination, Operational Resource
Semantics, and the Universality Boundary} (Release 10.32; Version 1.0.0)
[Preprint]. Zenodo.

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
    body = common_body(
        text,
        FIGURES,
        "FIN_Programs_365_378_Figures",
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
