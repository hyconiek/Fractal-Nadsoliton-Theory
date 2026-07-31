#!/usr/bin/env python3
"""Build the archival English PDF source for FIN Programs P351--P364."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body


SOURCE = Path("FIN_Programs_351_364_Rigorous_Operational_Completion_Report_EN.md")
TARGET = Path("FIN_Programs_351_364_Rigorous_Operational_Completion_Report_EN.tex")

FIGURES = {
    "p351_p352_rigorous_enclosures.png": (
        "The theorem-grade primal-dual enclosures for the attenuation and "
        "full oscillatory signed-moment resources."
    ),
    "p354_p355_naturality_loss.png": (
        "The exact graph-morphism naturality boundary and the declared "
        "loss-aware photonic digital twin."
    ),
    "p357_p358_comb_clock.png": (
        "Convex-zero multi-use channel witnesses and the bounded-curvature "
        "clock-design envelope."
    ),
    "p360_p362_phase_metrology.png": (
        "Equal-modulus analytic inner-factor phases and conditional "
        "conversion-frame uncertainty."
    ),
}


PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{
  unicode=true,
  colorlinks=true,
  linkcolor=finblue,
  citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Research Programs P351--P364},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Rigorous moment certificates, graph naturality, operational resource completion, and external physics gates},
  pdfkeywords={FIN, signed moments, Krawczyk method, Bernstein certificate, graph naturality, quantum comb, photonic mesh, resource category, Blaschke factor, operational physics}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs P351--P364}
\fancyhead[R]{\small Release 10.31}
\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{7mm}
{\Large\bfseries FIN --- Release 10.31\par}
\vspace{9mm}
{\Huge\bfseries Research Programs P351--P364\par}
\vspace{6mm}
{\Large Rigorous moment certificates, graph naturality,\\
operational resource completion, and external physics gates\par}
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
\textbf{Certified attenuation resource.}
\[
0.4067063265581317
\leq\mathcal N_{\rm env}
\leq0.4067063346922584.
\]

\textbf{Certified oscillatory resource.}
\[
0.7072782253162320
\leq\mathcal N_{\rm osc}
\leq0.7073599748862595.
\]

Both are classical signed-representation costs. No physical-negativity,
information-loss, selector, unit, role-transfer, laboratory, total-action,
Standard-Model, gravity, or Theory-of-Everything closure is claimed.

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
This monograph reports every executed Program P351--P364 and recommends
Programs P365--P378. Exact proof, strong numerical evidence, synthetic
apparatus design, conditional physics, and missing external evidence are
kept distinct. The legacy/strict split and all current strict-core
guardrails remain binding.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""


POSTAMBLE = r"""
\clearpage
\chapter*{Reproduction certificate}
\addcontentsline{toc}{chapter}{Reproduction certificate}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/mpl-fin-1031 \
python3 fin_programs_351_364.py

MPLCONFIGDIR=/tmp/mpl-fin-1031 \
python3 -m unittest -v test_fin_programs_351_364.py
\end{Verbatim}
Expected result: 16 tests passed, 0 failed. The dependency-free Lean 4.28
structural core must compile. P353, P356, P361, and P364 remain external
admission gates; no synthetic record is accepted as laboratory evidence.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Research Programs P351--P364:
Rigorous Moment Certificates, Graph Naturality, Operational Resource
Completion, and External Physics Gates} (Release 10.31; Version 1.0.0)
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
        "FIN_Programs_351_364_Figures",
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
