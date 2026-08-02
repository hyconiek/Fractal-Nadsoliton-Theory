#!/usr/bin/env python3
"""Build the archival TeX source for FIN checkpoint P474/P479/P483."""

from pathlib import Path
import re

from fin_english_reports_to_latex import COMMON_PACKAGES, COMMON_SETUP, common_body


SOURCE = Path("FIN_Local_Research_Checkpoint_P474_P479_P483_EN.md")
TARGET = Path("FIN_Local_Research_Checkpoint_P474_P479_P483_EN.tex")
FIGURES = {
    "p474_p483_optimal_face_and_parameter_tube.png": (
        "Positivity along the candidate P474 complex optimal direction and the "
        "exactly admitted/rejected P483 uniform Krawczyk tubes."
    ),
}

PREAMBLE = COMMON_PACKAGES + r"""
\hypersetup{unicode=true,colorlinks=true,linkcolor=finblue,citecolor=fingreen,
urlcolor=finblue,pdftitle={FIN Local Research Checkpoint P474-P479-P483},
pdfauthor={Krzysztof Żuchowski},
pdfsubject={Complex optimal-face evidence and exact parametric O181 continuation},
pdfkeywords={FIN, quantum comb, optimal face, Riccati equation, Krawczyk operator, interval arithmetic, parametric continuation}}
\pagestyle{fancy}\fancyhf{}
\fancyhead[L]{\small FIN Checkpoint P474/P479/P483}
\fancyhead[R]{\small Release 10.45}\fancyfoot[C]{\thepage}
""" + COMMON_SETUP + r"""
\begin{document}\hypersetup{pageanchor=false}
\begin{titlepage}\centering\vspace*{8mm}
{\Large\bfseries FIN Local Research --- Release 10.45\par}\vspace{10mm}
{\Huge\bfseries Checkpoint P474--P479--P483\par}\vspace{6mm}
{\Large Complex Optimal-Face Evidence and an Exact Parametric O181 Tube\par}\vspace{14mm}
{\Large Krzysztof \.Zuchowski\par}\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
{\normalsize ORCID: 0009-0002-0909-3613\par}\vspace{8mm}
{\large 1 August 2026\par}{\normalsize Publication --- Preprint; Version 1.0.0\par}
\vfill
\begin{minipage}{0.92\textwidth}\small
\textbf{Exact result.} A uniform rational Krawczyk inclusion proves a unique
positive O181 polynomial root in one common box for every
\[
|q-4/5|\le3\times10^{-9},\qquad
|\theta-\pi/8|\le3\times10^{-9}.
\]
The corresponding O167 point attains the global value of the declared
reduced three-slot ordered-comb SDP throughout this rectangle.

\medskip\textbf{Optimal-face result.} Two independent finite linearizations
show one candidate complex flat direction through the real optimum. This is
strong numerical evidence, not yet an exact nonuniqueness theorem.

\medskip\textbf{Boundary.} No selector, dimensional source, laboratory
record, complete legacy--strict bridge, role transfer, or physical closure is
claimed. The exported Lean source was not machine checked because no local
toolchain is installed.

\medskip\textbf{License.} CC BY 4.0.
\end{minipage}\vfill\end{titlepage}
\hypersetup{pageanchor=true}\pagenumbering{roman}
\chapter*{Publication statement}\addcontentsline{toc}{chapter}{Publication statement}
This monograph reports local analytical and computer-assisted Programs P474,
P479, and P483. It uses no laboratory record, external audit, internet
result, remote computation, or physical validation.
\tableofcontents\clearpage\pagenumbering{arabic}
"""

POSTAMBLE = r"""
\clearpage\chapter*{Reproduction certificate}
\addcontentsline{toc}{chapter}{Reproduction certificate}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/fin-mpl-1045 \
python3 fin_programs_474_479_483.py
python3 -m unittest -v test_fin_programs_474_479_483.py
python3 fin_programs_474_479_483_to_latex.py
lualatex -interaction=nonstopmode -halt-on-error \
  FIN_Local_Research_Checkpoint_P474_P479_P483_EN.tex
sha256sum -c FIN_PROGRAMS_474_479_483_RELEASE_MANIFEST.sha256
\end{Verbatim}
\noindent\textbf{Suggested citation:}
\.Zuchowski, K. (2026). \emph{FIN Local Research Checkpoint
P474--P479--P483: Complex Optimal-Face Evidence and an Exact Parametric O181
Tube} (Release 10.45; Version 1.0.0) [Preprint].
\end{document}
"""


def repair_inline_math(source: str) -> str:
    """Recover inline TeX delimiters lost while the Markdown was assembled.

    Display blocks are left untouched.  Outside them, only parenthesized
    spans carrying an unmistakable TeX/math marker (or one declared scalar
    name) are promoted, so ordinary prose parentheses remain prose.
    """

    output: list[str] = []
    in_display = False
    scalar = re.compile(r"^[KNQqrt]$")
    candidate = re.compile(r"\(([^()\n]+)\)")

    def replace(match: re.Match[str]) -> str:
        content = match.group(1)
        math_marker = (
            "\\" in content or "_" in content or "^" in content
            or "=" in content or content.startswith((">", "<"))
            or scalar.fullmatch(content) is not None
        )
        return rf"\({content}\)" if math_marker else match.group(0)

    for line in source.splitlines(keepends=True):
        stripped = line.strip()
        if stripped == r"\[":
            in_display = True
            output.append(line)
        elif stripped == r"\]":
            in_display = False
            output.append(line)
        elif in_display or "![" in line:
            output.append(line)
        else:
            output.append(candidate.sub(replace, line))
    return "".join(output)


def main() -> None:
    text = repair_inline_math(SOURCE.read_text(encoding="utf-8"))
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
        "FIN_Programs_474_479_483_Figures",
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
