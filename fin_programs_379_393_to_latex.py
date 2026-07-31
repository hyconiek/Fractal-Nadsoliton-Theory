#!/usr/bin/env python3
"""Build the archival English PDF source for FIN Programs P379--P393."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body


SOURCE = Path("FIN_Programs_379_393_Proof_Reflected_Operational_Report_EN.md")
TARGET = Path("FIN_Programs_379_393_Proof_Reflected_Operational_Report_EN.tex")

FIGURES = {
    "p379_p380_reflection_contact.png": (
        "Lean-reflected terminal predicates and the reduction of the exact "
        "oscillatory primal-dual gap."
    ),
    "p381_p383_injectivity_gauge.png": (
        "Global strict radial injectivity diagnostics and the photonic "
        "component-phase coordinate gauge."
    ),
    "p385_p386_loss_clock.png": (
        "The declared noisy-comb lower/upper sandwich and propagated "
        "clock-calibration tubes."
    ),
    "p387_p389_realization_scale.png": (
        "Positive Jordan-sampling preparation probabilities and the "
        "noncanonical outer-scale dilation orbits."
    ),
}


PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{
  unicode=true,
  colorlinks=true,
  linkcolor=finblue,
  citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={FIN Research Programs P379--P393},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Proof-reflected certificates, global radial injectivity, noisy dynamics, and a first explicit operational realization},
  pdfkeywords={FIN, Lean, exact certificate, Bernstein polynomial, Lindemann-Weierstrass, radial injectivity, graph isometry, quantum comb, importance sampling, operational realization, dilation semigroup}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small FIN Programs P379--P393}
\fancyhead[R]{\small Release 10.33}
\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{7mm}
{\Large\bfseries FIN --- Release 10.33\par}
\vspace{9mm}
{\Huge\bfseries Research Programs P379--P393\par}
\vspace{6mm}
{\Large Proof-reflected certificates, global radial injectivity,\\
noisy dynamics, and a first explicit operational realization\par}
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
\textbf{Near-contact oscillatory certificate.}
\[
0.7073534379053974
\leq\mathcal N_{\rm osc}
\leq0.7073534683998260,
\qquad
\operatorname{width}<3.05\times10^{-8}.
\]

\textbf{Global strict radial theorem.}
\[
d\neq e,\quad d,e\in\mathbb N_0
\quad\Longrightarrow\quad
K_{\rm strict}(d)\neq K_{\rm strict}(e).
\]

Thirty-five exact terminal predicates are checked by Lean 4.28. One
positive-probability mathematical realization is constructed for the signed
resource. No selector, dimensional standard, calibrated laboratory
realization, role transfer, Standard Model, gravity, or Theory-of-Everything
closure is claimed.

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
This monograph reports every executed Program P379--P393 and recommends
Programs P394--P410. Exact proof, proof-kernel reflection, bounded numerical
evidence, conditional operational construction, and absent external evidence
are kept distinct. The legacy/strict split and all current strict-core
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
MPLCONFIGDIR=/tmp/mpl-fin-1033 \
python3 fin_programs_379_393.py

MPLCONFIGDIR=/tmp/mpl-fin-1033 \
python3 -m unittest -v test_fin_programs_379_393.py
\end{Verbatim}
Expected result: 16 tests passed, 0 failed. Both Lean 4.28 sources must
compile. P384 and P390--P393 remain external admission gates; no synthetic
record is accepted as hardware, standards, hold-out, reservoir, electroweak,
or custody evidence.

\medskip
\noindent\textbf{Suggested citation:}

\.Zuchowski, K. (2026). \emph{FIN Research Programs P379--P393:
Proof-Reflected Certificates, Global Radial Injectivity, Noisy Dynamics, and
a First Explicit Operational Realization} (Release 10.33; Version 1.0.0)
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
        "FIN_Programs_379_393_Figures",
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
