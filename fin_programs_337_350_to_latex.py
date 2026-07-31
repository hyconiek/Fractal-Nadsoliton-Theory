#!/usr/bin/env python3
"""Build the archival English PDF source for FIN Programs P337--P350."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body


SOURCE = Path("FIN_Programs_337_350_Naturality_Resource_Report_EN.md")
TARGET = Path("FIN_Programs_337_350_Naturality_Resource_Report_EN.tex")

FIGURES = {
    "p337_p338_moment_certificates.png": (
        "The exact rational envelope lower certificate and the full "
        "phase-bearing strict moment sequence."
    ),
    "p339_p340_refreeze_naturality.png": (
        "Continuous QW-2038 sensitivity and failure of the C12-trained scalar "
        "map on other carriers."
    ),
    "p341_p343_photonic_comb.png": (
        "Synthetic short-time inverse calibration and ideal parallel-use "
        "channel discrimination."
    ),
    "p344_p345_clock_invariants.png": (
        "Nonparametric clock uncertainty and carrier-dependent operator "
        "invariants."
    ),
    "p346_p347_resources_couplings.png": (
        "Signed-measure monotonicity under Markov maps and the tested "
        "resource-coupling matrix."
    ),
    "p348_p349_conditional_metrology.png": (
        "The conditional electroweak vector and rank-three dimensional "
        "conversion metrology."
    ),
}


PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{
  unicode=true,
  colorlinks=true,
  linkcolor=finblue,
  citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Research Programs P337--P350},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Certified signed resources, carrier naturality, operational compilation, and conditional physics interfaces},
  pdfkeywords={FIN, signed moment, Bernstein certificate, oscillatory kernel, carrier naturality, photonic mesh, quantum channel, clock identifiability, resource theory, dimensional metrology}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs P337--P350}
\fancyhead[R]{\small Release 10.30}
\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{7mm}
{\Large\bfseries FIN --- Release 10.30\par}
\vspace{9mm}
{\Huge\bfseries Research Programs P337--P350\par}
\vspace{6mm}
{\Large Certified signed resources, carrier naturality,\\
operational compilation, and conditional physics interfaces\par}
\vspace{13mm}
{\Large Krzysztof \.Zuchowski\par}
\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
\vspace{2mm}
{\normalsize ORCID:
\href{https://orcid.org/0009-0002-0909-3613}{0009-0002-0909-3613}\par}
\vspace{7mm}
{\large 30 July 2026\par}
{\normalsize Publication --- Preprint; Version 1.0.0\par}
\vfill
\begin{minipage}{0.92\textwidth}
\small
\textbf{Certified envelope result.}
Exact rational Bernstein and algebraic-input interval arithmetic proves
\(\mathcal N_{\rm env}\geq0.4067063265581317\).

\medskip
\textbf{Full oscillatory result.}
For the twelve phase-bearing strict moments,
\[
0.7072784747\lesssim\mathcal N_{\rm osc}
\lesssim0.7073599749
\]
is the strong numerical enclosure.  It is a classical signed-path resource,
not a physical-negativity theorem.

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
This monograph reports all executed Programs P337--P350 and recommends
Programs P351--P364.  Exact proof, high-precision evaluation, synthetic
apparatus design, conditional sector physics, dimensional calibration, and
missing external evidence remain separate.  The legacy/strict split,
QW-2191, dimensional-source, role-transfer, laboratory, \(L_{\rm total}\),
Standard-Model/gravity, and Theory-of-Everything boundaries are preserved.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""


POSTAMBLE = r"""
\clearpage
\chapter*{Reproduction certificate}
\addcontentsline{toc}{chapter}{Reproduction certificate}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/mpl-fin \
python3 fin_programs_337_350.py

MPLCONFIGDIR=/tmp/mpl-fin \
python3 -m unittest -v test_fin_programs_337_350.py
\end{Verbatim}
Expected result: 18 tests passed, 0 failed.  The suite recompiles the
dependency-free Lean 4.28 structural core.  Refreeze search, photonic
calibration, detector models, clocks, and metrology designs retain the
boundaries stated in the text.  P342 and P350 admit no external evidence.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Research Programs P337--P350:
Certified Signed Resources, Carrier Naturality, Operational Compilation, and
Conditional Physics Interfaces} (Release 10.30; Version 1.0.0) [Preprint].
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
    body = common_body(
        text,
        FIGURES,
        "FIN_Programs_337_350_Figures",
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
