#!/usr/bin/env python3
"""Convert the FIN twenty-program Markdown monograph to archival LaTeX.

The converter is deliberately small and deterministic.  It uses markdown-it
only for structural parsing and preserves TeX mathematics already delimited by
``\\(...\\)`` and ``\\[...\\]``.  It supports the exact constructs used by the
source monograph: headings, paragraphs, emphasis, links, lists, block quotes,
horizontal rules, one code block, and pipe tables.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

from markdown_it import MarkdownIt


SOURCE = Path("FIN_Ten_Research_Programs_Monograph.md")
TARGET = Path("FIN_Release_10_1_Can_FIN_Become_Predictive_Physics.tex")


PREAMBLE = r"""\documentclass[11pt,a4paper,openany]{report}

\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{microtype}
\usepackage{amsmath,amssymb,amsthm,mathtools,bm}
\usepackage{booktabs,array,longtable,tabularx}
\usepackage{xcolor}
\usepackage[margin=24mm,headheight=23pt]{geometry}
\usepackage{enumitem}
\usepackage{fancyhdr}
\usepackage{fancyvrb}
\usepackage{url}
\usepackage{hyperref}

\definecolor{finblue}{HTML}{1F5A99}
\definecolor{fingreen}{HTML}{19733A}
\definecolor{finorange}{HTML}{D55E00}
\definecolor{finred}{HTML}{A61B1B}
\definecolor{finviolet}{HTML}{6A3D9A}
\definecolor{fingray}{HTML}{4D4D4D}

\hypersetup{
  unicode=true,
  colorlinks=true,
  linkcolor=finblue,
  citecolor=fingreen,
  urlcolor=finblue,
  pdftitle={Can FIN Become Predictive Physics?},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Comprehensive mathematical synthesis of FIN: twenty executed research programs, the wave-diffusion observer problem, and ten next recommendations — Release 10.1},
  pdfkeywords={FIN, spectral graph theory, operator theory, unitary dynamics, Markov semigroup, operational physics, process tensors, mathematical physics}
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small\nouppercase{\leftmark}}
\fancyhead[R]{\small FIN Release 10.1}
\fancyfoot[C]{\thepage}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.55em}
\setlist{nosep,leftmargin=2em}
\setcounter{tocdepth}{2}
\setcounter{secnumdepth}{3}
\emergencystretch=2em

\newcommand{\one}{\bm 1}
\newcommand{\C}{\mathbb C}
\newcommand{\R}{\mathbb R}
\newcommand{\Z}{\mathbb Z}
\newcommand{\statusProven}{\textcolor{fingreen}{[Proven]}}
\newcommand{\statusStrong}{\textcolor{finblue}{[Strong evidence]}}
\newcommand{\statusModerate}{\textcolor{finviolet}{[Moderate evidence]}}
\newcommand{\statusWeak}{\textcolor{fingray}{[Weak evidence]}}
\newcommand{\statusSpeculative}{\textcolor{finorange}{[Speculative]}}
\newcommand{\statusRefuted}{\textcolor{finred}{[Refuted]}}

\begin{document}
\hypersetup{pageanchor=false}
\begin{titlepage}
\centering
\vspace*{15mm}
{\Large\bfseries FIN Theory Synthesis — Release 10.1\par}
\vspace{12mm}
{\Huge\bfseries Can FIN Become Predictive Physics?\par}
\vspace{6mm}
{\Large Twenty Executed Research Programs, the Wave--Diffusion Observer Problem, and Ten Next Recommendations\par}
\vspace{18mm}
{\Large Krzysztof Żuchowski\par}
\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
\vspace{2mm}
{\normalsize ORCID: \href{https://orcid.org/0009-0002-0909-3613}{0009-0002-0909-3613}\par}
\vspace{10mm}
{\large 19 July 2026\par}
{\normalsize Publication --- Preprint; record version 1.0.0\par}
\vfill
\begin{minipage}{0.88\textwidth}
\small
\textbf{Scope.} Mathematical physics only. Cosmological, metaphysical and Theory-of-Everything interpretations are excluded. The purpose is falsification, reconstruction and identification of the smallest operational and mathematical structures required for predictive physics.

\medskip
\textbf{Repository.} \url{https://github.com/hyconiek/Fractal-Nadsoliton-Theory}

\textbf{Reproducibility.} The release includes deterministic core and observer-dynamics suites.

\textbf{License.} Creative Commons Attribution 4.0 International (CC BY 4.0).

\textbf{Zenodo DOI.} To be assigned to Release 10.1.
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}

\pagenumbering{roman}
\chapter*{Publication statement}
\addcontentsline{toc}{chapter}{Publication statement}
This release is a single self-contained synthesis of the present FIN mathematical theory. It reconstructs the minimal spectral core, the strict--legacy relation, variational and informational structures, adaptive dynamics, observer/apparatus requirements, unitary and diffusive temporal semantics, continuum obstructions and the bridge still required for predictive physics. It reports twenty executed mathematical research programs, a cross-field falsification search, ten recommended next programs, a dependency roadmap, and a final scientific judgment. Confidence labels distinguish proof, numerical evidence, conjecture and refutation throughout.

The exact source is \path{FIN_Ten_Research_Programs_Monograph.md}. Numerical claims are reproducible with \path{fin_ten_programs_experiments.py} and \path{fin_observer_dynamics_experiments.py}. The scripts use explicit finite kernel data and NumPy; they do not import generated FIN verdicts.

\tableofcontents
\clearpage
\pagenumbering{arabic}
"""


POSTAMBLE = r"""
\clearpage
\chapter*{Release reproducibility}
\addcontentsline{toc}{chapter}{Release reproducibility}
The numerical suites accompanying this monograph can be executed with:
\begin{Verbatim}[fontsize=\small]
python3 fin_ten_programs_experiments.py
python3 fin_observer_dynamics_experiments.py
\end{Verbatim}
All times in the finite observer calculations are dimensionless unless a calibration is explicitly supplied. Floating-point residuals must be interpreted with the tolerances declared by the scripts, not as physical effects.

\medskip
\noindent\textbf{Suggested citation before DOI assignment:}

Żuchowski, K. (2026). \emph{Can FIN Become Predictive Physics? Twenty Executed Research Programs, the Wave--Diffusion Observer Problem, and Ten Next Recommendations} (Release 10.1; Version 1.0.0) [Preprint]. Zenodo. DOI to be assigned.

\end{document}
"""


STATUS = {
    "[Proven]": r"\statusProven{}",
    "[Strong evidence]": r"\statusStrong{}",
    "[Moderate evidence]": r"\statusModerate{}",
    "[Weak evidence]": r"\statusWeak{}",
    "[Speculative]": r"\statusSpeculative{}",
    "[Refuted]": r"\statusRefuted{}",
}

MATH_REPLACEMENTS: dict[str, str] = {}
MATH_TOKEN_RE = re.compile(r"FINMATH\d{6}TOKEN")


def protect_math(source: str) -> str:
    """Replace TeX math spans before CommonMark consumes backslash escapes."""
    MATH_REPLACEMENTS.clear()
    pattern = re.compile(r"\\\[.*?\\\]|\\\(.*?\\\)", re.DOTALL)

    def replace(match: re.Match[str]) -> str:
        token = f"FINMATH{len(MATH_REPLACEMENTS):06d}TOKEN"
        MATH_REPLACEMENTS[token] = match.group(0)
        return token

    return pattern.sub(replace, source)


def escape_normal(text: str) -> str:
    """Escape prose while preserving confidence markers."""
    out: list[str] = []
    i = 0
    markers = sorted(STATUS, key=len, reverse=True)
    while i < len(text):
        math_match = MATH_TOKEN_RE.match(text, i)
        if math_match:
            token = math_match.group(0)
            out.append(MATH_REPLACEMENTS[token])
            i = math_match.end()
            continue
        marker = next((m for m in markers if text.startswith(m, i)), None)
        if marker:
            out.append(STATUS[marker])
            i += len(marker)
            continue
        ch = text[i]
        replacements = {
            "&": r"\&",
            "%": r"\%",
            "$": r"\$",
            "#": r"\#",
            "_": r"\_",
            "{": r"\{",
            "}": r"\}",
            "~": r"\textasciitilde{}",
            "^": r"\textasciicircum{}",
            "\\": r"\textbackslash{}",
            "/": r"/\allowbreak{}",
        }
        out.append(replacements.get(ch, ch))
        i += 1
    return "".join(out)


def escape_url(url: str) -> str:
    return url.replace("%", r"\%").replace("#", r"\#")


def render_inline(token) -> str:
    out: list[str] = []
    mode = "normal"
    link_stack: list[str] = []

    def text_piece(text: str) -> None:
        nonlocal mode
        i = 0
        buffer: list[str] = []

        def flush() -> None:
            if buffer:
                out.append(escape_normal("".join(buffer)))
                buffer.clear()

        while i < len(text):
            if mode == "normal" and text.startswith(r"\(", i):
                flush()
                out.append("$")
                mode = "inline_math"
                i += 2
            elif mode == "normal" and text.startswith(r"\[", i):
                flush()
                out.append(r"\[")
                mode = "display_math"
                i += 2
            elif mode == "inline_math" and text.startswith(r"\)", i):
                out.append("$")
                mode = "normal"
                i += 2
            elif mode == "display_math" and text.startswith(r"\]", i):
                out.append(r"\]")
                mode = "normal"
                i += 2
            elif mode == "normal":
                buffer.append(text[i])
                i += 1
            else:
                out.append(text[i])
                i += 1
        flush()

    for child in token.children or []:
        kind = child.type
        if kind == "text":
            text_piece(child.content)
        elif kind == "strong_open":
            out.append(r"\textbf{")
        elif kind == "strong_close":
            out.append("}")
        elif kind == "em_open":
            out.append(r"\emph{")
        elif kind == "em_close":
            out.append("}")
        elif kind == "link_open":
            href = child.attrGet("href") or ""
            link_stack.append(href)
            out.append(r"\href{" + escape_url(href) + "}{")
        elif kind == "link_close":
            if link_stack:
                link_stack.pop()
            out.append("}")
        elif kind == "softbreak":
            out.append("\n" if mode == "display_math" else " ")
        elif kind == "hardbreak":
            out.append(r"\\" + "\n")
        else:
            raise ValueError(f"Unsupported inline token: {kind}")
    if mode != "normal":
        raise ValueError(f"Unclosed math delimiter in inline token: {token.content[:80]!r}")
    return "".join(out)


def column_spec(rows: list[list[str]]) -> str:
    count = len(rows[0])

    def visible(cell: str) -> str:
        text = cell.replace(r"\allowbreak{}", "")
        text = re.sub(r"\\[A-Za-z]+", "", text)
        return text.replace("{", "").replace("}", "").replace("$", "")

    weights: list[float] = []
    for index in range(count):
        cells = [visible(row[index]) for row in rows]
        longest_word = max((len(word) for cell in cells for word in re.split(r"[\s/]+", cell) if word), default=1)
        longest_cell = max((len(cell) for cell in cells), default=1)
        weights.append(max(3.0, longest_word + 0.15 * min(longest_cell, 50)))

    target = 0.86
    minimum = 0.045 if count >= 7 else 0.055
    maximum = 0.42
    widths: list[float | None] = [None] * count
    active = set(range(count))
    remaining = target
    while active:
        total_weight = sum(weights[index] for index in active)
        changed = False
        for index in list(active):
            proposal = remaining * weights[index] / total_weight
            if proposal < minimum:
                widths[index] = minimum
                remaining -= minimum
                active.remove(index)
                changed = True
            elif proposal > maximum:
                widths[index] = maximum
                remaining -= maximum
                active.remove(index)
                changed = True
        if not changed:
            total_weight = sum(weights[index] for index in active)
            for index in active:
                widths[index] = remaining * weights[index] / total_weight
            break
    final_widths = [float(width) for width in widths]
    return "@{}" + "".join(f"p{{{w:.3f}\\textwidth}}" for w in final_widths) + "@{}"


def render_table(tokens, start: int) -> tuple[str, int]:
    rows: list[list[str]] = []
    header_rows = 0
    current: list[str] | None = None
    in_header = False
    i = start + 1
    while i < len(tokens):
        token = tokens[i]
        if token.type == "table_close":
            break
        if token.type == "thead_open":
            in_header = True
        elif token.type == "thead_close":
            in_header = False
        elif token.type == "tr_open":
            current = []
        elif token.type == "tr_close":
            if current is not None:
                rows.append(current)
                if in_header:
                    header_rows += 1
            current = None
        elif token.type == "inline" and current is not None:
            current.append(render_inline(token))
        i += 1
    if not rows:
        return "", i + 1
    count = max(len(row) for row in rows)
    for row in rows:
        row.extend([""] * (count - len(row)))
    table_size = r"\fontsize{6.5}{7.8}\selectfont" if count >= 7 else (r"\footnotesize" if count >= 5 else r"\small")
    lines = [
        r"\begin{center}",
        r"\begingroup" + table_size,
        r"\setlength{\tabcolsep}{3.5pt}",
        r"\renewcommand{\arraystretch}{1.16}",
        r"\begin{longtable}{" + column_spec(rows) + "}",
        r"\toprule",
    ]
    for row_index, row in enumerate(rows):
        rendered = [r"\textbf{" + cell + "}" if row_index < header_rows else cell for cell in row]
        lines.append(" & ".join(rendered) + r" \\")
        if row_index + 1 == header_rows:
            lines.extend([r"\midrule", r"\endhead"])
    lines.extend([
        r"\bottomrule",
        r"\end{longtable}",
        r"\endgroup",
        r"\end{center}",
        "",
    ])
    return "\n".join(lines), i + 1


def ascii_code(text: str) -> str:
    translation = str.maketrans(
        {
            "├": "+", "┤": "+", "└": "+", "┌": "+", "┐": "+", "┘": "+",
            "┬": "+", "┴": "+", "┼": "+", "─": "-", "│": "|",
            "▼": "v", "→": "->", "←": "<-", "—": "--", "–": "-",
        }
    )
    return text.translate(translation)


def body_from_markdown(source: str) -> str:
    # The archival title page replaces the short Markdown metadata header.
    marker = "## Confidence convention"
    if marker not in source:
        raise ValueError("Expected confidence-convention marker not found")
    source = marker + source.split(marker, 1)[1]
    source = protect_math(source)
    md = MarkdownIt("commonmark").enable("table")
    tokens = md.parse(source)
    lines: list[str] = []
    i = 0
    list_depth = 0
    while i < len(tokens):
        token = tokens[i]
        kind = token.type
        if kind == "table_open":
            rendered, i = render_table(tokens, i)
            lines.append(rendered)
            continue
        if kind == "heading_open":
            level = int(token.tag[1])
            inline = tokens[i + 1]
            title = render_inline(inline)
            if "FINMATH" in inline.content:
                plain = inline.content
                for key, value in MATH_REPLACEMENTS.items():
                    if key in plain:
                        plain = plain.replace(key, value[2:-2])
                plain = re.sub(r"\\[A-Za-z]+", "", plain).replace("{", "").replace("}", "")
                title = r"\texorpdfstring{" + title + "}{" + escape_normal(plain) + "}"
            command = {1: "chapter", 2: "section", 3: "subsection", 4: "subsubsection"}.get(level, "paragraph")
            lines.append(f"\\{command}{{{title}}}\n")
            i += 3
            continue
        if kind == "paragraph_open":
            inline = tokens[i + 1]
            lines.append(render_inline(inline) + "\n")
            i += 3
            continue
        if kind == "bullet_list_open":
            lines.append(r"\begin{itemize}")
            list_depth += 1
        elif kind == "bullet_list_close":
            lines.append(r"\end{itemize}" + "\n")
            list_depth -= 1
        elif kind == "ordered_list_open":
            lines.append(r"\begin{enumerate}")
            list_depth += 1
        elif kind == "ordered_list_close":
            lines.append(r"\end{enumerate}" + "\n")
            list_depth -= 1
        elif kind == "list_item_open":
            lines.append(r"\item ")
        elif kind == "list_item_close":
            pass
        elif kind == "blockquote_open":
            lines.append(r"\begin{quote}\itshape")
        elif kind == "blockquote_close":
            lines.append(r"\end{quote}" + "\n")
        elif kind == "hr":
            lines.append(r"\par\medskip\hrule\medskip" + "\n")
        elif kind == "code_block":
            lines.extend([
                r"\begin{Verbatim}[fontsize=\small]",
                ascii_code(token.content).rstrip(),
                r"\end{Verbatim}",
                "",
            ])
        elif kind in {
            "inline", "heading_close", "paragraph_close", "table_close",
            "thead_open", "thead_close", "tbody_open", "tbody_close",
            "tr_open", "tr_close", "th_open", "th_close", "td_open", "td_close",
        }:
            pass
        else:
            raise ValueError(f"Unsupported block token: {kind}")
        i += 1
    if list_depth != 0:
        raise ValueError(f"Unbalanced list nesting: {list_depth}")
    return "\n".join(lines)


def main() -> None:
    source = SOURCE.read_text(encoding="utf-8")
    body = body_from_markdown(source)
    TARGET.write_text(PREAMBLE + "\n" + body + "\n" + POSTAMBLE, encoding="utf-8")
    print(f"wrote {TARGET} ({TARGET.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
