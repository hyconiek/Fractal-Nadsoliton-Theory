#!/usr/bin/env python3
"""Build the archival English PDF source for FIN Programs P411--P427."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body


SOURCE = Path("FIN_Programs_411_427_Reference_Erasure_Report_EN.md")
TARGET = Path("FIN_Programs_411_427_Reference_Erasure_Report_EN.tex")

FIGURES = {
    "p411_p412_formal_contact.png": (
        "Exact reflected Taylor widths and the simultaneous seven-contact "
        "KKT candidate."
    ),
    "p414_p415_damping_photonic.png": (
        "The normalized nonlinear damping-source obstruction and sampled "
        "photonic Jacobian regularity."
    ),
    "p417_p418_noisy_erasure.png": (
        "The remaining noisy-comb certificate gap and gains from the new "
        "heralded-erasure-aware symmetric code."
    ),
    "p422_p425_estimator_scale.png": (
        "Variance-optimal Jordan sampling and inequivalent conditional "
        "outer-scale sections."
    ),
}


PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{
  unicode=true,
  colorlinks=true,
  linkcolor=finblue,
  citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Research Programs P411--P427},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Formal cosine boundary, simultaneous contact geometry, erasure-aware discrimination, and operational reference resources},
  pdfkeywords={FIN, Lean, Taylor theorem, KKT, signed moment problem, nonlinear recurrence, photonic identifiability, quantum discrimination, erasure, Jordan sampling, phase reference, orientation torsor, scale section}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs P411--P427}
\fancyhead[R]{\small Release 10.35}
\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{7mm}
{\Large\bfseries FIN --- Release 10.35\par}
\vspace{9mm}
{\Huge\bfseries Research Programs P411--P427\par}
\vspace{6mm}
{\Large Formal cosine boundary, simultaneous contact geometry,\\
erasure-aware discrimination, and operational reference resources\par}
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
\textbf{Constructive operational result.}
The heralded-erasure-aware symmetric code improves both product and GHZ
baselines in all 16 declared noise cells, with maximum gain approximately
\(0.02694\).

\medskip
\textbf{Exact estimator result.}
\[
q_i^\star\propto |w_i|
\sqrt{\sum_{k=0}^{11}x_i^{2k}},
\]
reducing the declared twelve-moment variance by approximately \(11.57\%\).

\medskip
No non-premise selector, dimensional source, laboratory record, complete
legacy-to-strict bridge, role transfer, Standard Model, gravity, or
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
This monograph reports every executed Program P411--P427 and recommends
Programs P428--P444. Mathematical proof, numerical evidence, conditional
operational structure, and absent external evidence are kept distinct.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""


POSTAMBLE = r"""
\clearpage
\chapter*{Reproduction certificate}
\addcontentsline{toc}{chapter}{Reproduction certificate}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/mpl-fin-1035 \
python3 fin_programs_411_427.py

python3 -m unittest -v test_fin_programs_411_427.py
\end{Verbatim}
Expected result: 16 tests passed, zero failed. P416, P420, P421, P426, and
P427 remain external gates. Synthetic or repository-generated records are
not accepted as physical evidence.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Research Programs P411--P427: Formal Cosine
Boundary, Simultaneous Contact Geometry, Erasure-Aware Discrimination, and
Operational Reference Resources} (Release 10.35; Version 1.0.0) [Preprint].
Zenodo.

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
    body = common_body(text, FIGURES, "FIN_Programs_411_427_Figures")
    body = body.replace(
        r"\section{Confidence convention}",
        r"\chapter*{Confidence convention}"
        "\n"
        r"\addcontentsline{toc}{chapter}{Confidence convention}",
    )
    TARGET.write_text(PREAMBLE + "\n" + body + "\n" + POSTAMBLE, encoding="utf-8")
    print(TARGET)


if __name__ == "__main__":
    main()
