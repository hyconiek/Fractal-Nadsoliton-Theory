#!/usr/bin/env python3
"""Convert the FIN second-generation puzzle atlas to archival LaTeX."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Atlas_Drugiej_Generacji_Pochodne_Obiekty_Nadsolitona_PL.md")
TARGET = Path("FIN_Atlas_Drugiej_Generacji_Pochodne_Obiekty_Nadsolitona_PL.tex")

PREAMBLE = r"""\documentclass[11pt,a4paper,openany]{report}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
\usepackage{microtype}
\usepackage{amsmath,amssymb,mathtools,bm}
\usepackage{booktabs,array,longtable,tabularx}
\usepackage{xcolor}
\usepackage[margin=22mm,headheight=23pt]{geometry}
\usepackage{enumitem}
\usepackage{fancyhdr}
\usepackage{fancyvrb}
\usepackage{graphicx}
\usepackage{url}
\usepackage{xurl}
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
  pdftitle={Atlas Drugiej Generacji: Pochodne Obiekty Nadsolitona},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Systematyczne przeszukanie kombinacji FIN C01-C12, obiekty O01-O15 i synteza G01-G10},
  pdfkeywords={FIN, nadsoliton, Stieltjes functions, Schur complement, memory kernel, operator theory, information geometry}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small Atlas Drugiej Generacji}
\fancyhead[R]{\small Release 10.23}
\fancyfoot[C]{\thepage}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.55em}
\setlist{nosep,leftmargin=2em}
\setcounter{tocdepth}{2}
\setcounter{secnumdepth}{3}
\emergencystretch=3em
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
\vspace*{12mm}
{\Large\bfseries FIN --- Release 10.23\par}
\vspace{12mm}
{\Huge\bfseries Atlas Drugiej Generacji\par}
\vspace{5mm}
{\LARGE\bfseries Pochodne Obiekty Nadsolitona\par}
\vspace{8mm}
{\Large Systematyczne przeszukanie C01...C12,\\
konstrukcja O01...O15 i synteza G01...G10\par}
\vspace{18mm}
{\Large Krzysztof \.Zuchowski\par}
\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
\vspace{2mm}
{\normalsize ORCID: \href{https://orcid.org/0009-0002-0909-3613}{0009-0002-0909-3613}\par}
\vspace{10mm}
{\large 28 lipca 2026\par}
\vfill
\begin{minipage}{0.9\textwidth}
\small
\textbf{Główny obiekt.}
\[
E\longmapsto \operatorname{Schur}_{V\setminus E}(zI+A),\qquad
\Sigma_E(z)=A_{EH}(zI+A_{HH})^{-1}A_{HE}.
\]
Kontekstowe redukcje składają się, a samoenergie są dodatnimi operatorowymi
funkcjami Stieltjesa generującymi dokładną pamięć.

\medskip
\textbf{Granica.} Raport nie eksportuje selektora, skali fizycznej,
mostu legacy--strict, role transfer ani walidacji laboratoryjnej.

\medskip
\textbf{Repozytorium.}
\url{https://github.com/hyconiek/Fractal-Nadsoliton-Theory}

\textbf{Licencja.} CC BY 4.0.
\end{minipage}
\vfill
\end{titlepage}
\hypersetup{pageanchor=true}
\pagenumbering{roman}
\tableofcontents
\clearpage
\pagenumbering{arabic}
"""

POSTAMBLE = r"""
\clearpage
\chapter*{Reprodukcja}
\addcontentsline{toc}{chapter}{Reprodukcja}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/matplotlib-second \
python3 fin_nadsoliton_second_generation_atlas.py
\end{Verbatim}
Oczekiwane pokrycie: 305 kombinacji, 15 obiektów generacji O i 10 obiektów
generacji G.  Gramatyka konstrukcji jest filtrem kandydatów, nie dowodem
równoważności ani fizyczności.
\end{document}
"""

FIGURES = {
    "combination_search_coverage.png": (
        "Pełne pokrycie singli, par i trójek oraz wynik jawnej gramatyki konstrukcji."
    ),
    "derived_object_score_matrix.png": (
        "Audyt punktowy O01--O15; ryzyko nadinterpretacji ma przeciwną semantykę do pozostałych ocen."
    ),
    "generation2_generation3_map.png": (
        "Zależności wejściowe między obiektami generacji O i G; strzałki nie oznaczają równoważności fizycznej."
    ),
}


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    text = re.sub(
        r"\$\$\s*(.*?)\s*\$\$",
        lambda match: "\\[\n" + match.group(1).strip() + "\n\\]",
        text,
        flags=re.DOTALL,
    )
    text = re.sub(
        r"(?<!\\)\$([^\n$]+?)(?<!\\)\$",
        lambda match: "\\(" + match.group(1) + "\\)",
        text,
    )
    text = (
        text.replace("–", "--")
        .replace("—", "---")
        .replace("’", "'")
        .replace("„", '"')
        .replace("“", '"')
        .replace("”", '"')
        .replace("→", "->")
        .replace("↔", "<->")
    )
    for filename in FIGURES:
        pattern = re.compile(
            r"!\[[^\]]*\]\(FIN_Second_Generation_Atlas_Figures/"
            + re.escape(filename)
            + r"\)"
        )
        text = pattern.sub(f"\n\nFINFIG_{filename}\n\n", text)
    text = re.sub(r"`([^`]+)`", r"\1", text)
    text = re.sub(
        r"^(#{2,4})(\s+)",
        lambda match: "#" * (len(match.group(1)) - 1) + match.group(2),
        text,
        flags=re.MULTILINE,
    )
    text = text.replace("# Confidence convention", "## Confidence convention", 1)
    body = converter.body_from_markdown(text)
    # The atlas contains several semantically dense 5--6 column audit tables.
    # Keep them inside the text block without changing their mathematical
    # content; the generic converter's footnote-sized default is too wide here.
    body = body.replace(
        r"\begingroup\footnotesize",
        r"\begingroup\fontsize{7}{8.4}\selectfont",
    )
    body = body.replace(
        r"\setlength{\tabcolsep}{3.5pt}",
        r"\setlength{\tabcolsep}{2.5pt}",
    )
    body = re.sub(
        r"(C\d{2})\+(?=C\d{2})",
        lambda match: match.group(1) + r"+\allowbreak{}",
        body,
    )
    body = body.replace(
        r"\section{Confidence convention}",
        r"\chapter*{Confidence convention}"
        "\n"
        r"\addcontentsline{toc}{chapter}{Confidence convention}",
    )
    body = body.replace("Ż", r"\.Z").replace("ż", r"\.z")
    for filename, caption in FIGURES.items():
        token = f"FINFIG\\_{filename.replace('_', r'\_')}"
        path = f"FIN_Second_Generation_Atlas_Figures/{filename}"
        body = body.replace(
            token,
            "\\begin{figure}[htbp]\n"
            "\\centering\n"
            f"\\includegraphics[width=0.97\\textwidth]{{{path}}}\n"
            f"\\caption{{{caption}}}\n"
            "\\end{figure}",
        )
    TARGET.write_text(PREAMBLE + "\n" + body + "\n" + POSTAMBLE, encoding="utf-8")
    print(TARGET)


if __name__ == "__main__":
    main()
