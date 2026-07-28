#!/usr/bin/env python3
"""Convert the Polish FIN Z12 puzzle-atlas report to archival LaTeX."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Atlas_Puzzli_Nadsolitona_Z12_Raport_PL.md")
TARGET = Path("FIN_Atlas_Puzzli_Nadsolitona_Z12_Raport_PL.tex")

PREAMBLE = r"""\documentclass[11pt,a4paper,openany]{report}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
\usepackage{microtype}
\usepackage{amsmath,amssymb,mathtools,bm}
\usepackage{booktabs,array,longtable,tabularx}
\usepackage{xcolor}
\usepackage[margin=23mm,headheight=23pt]{geometry}
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
  pdftitle={Atlas puzzli nadsolitona: audyt Z12 i dynamiczna redukcja},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Kontrolowana analiza międzydziedzinowa FIN, audyt symulatora Z12 i twierdzenie o pamięci generowanej przez dynamiczny Schur},
  pdfkeywords={FIN, Z12, spectral graph, Schur complement, memory kernel, Mori-Zwanzig, chirality, operator theory}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small Atlas puzzli nadsolitona}
\fancyhead[R]{\small Release 10.22}
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
{\Large\bfseries FIN --- Release 10.22\par}
\vspace{12mm}
{\Huge\bfseries Atlas puzzli nadsolitona\par}
\vspace{6mm}
{\Large Audyt symulatora Z12, kontrolowana intuicja międzydziedzinowa\\
i nowy obiekt pamięci po redukcji\par}
\vspace{18mm}
{\Large Krzysztof \.Zuchowski\par}
\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
\vspace{2mm}
{\normalsize ORCID: \href{https://orcid.org/0009-0002-0909-3613}{0009-0002-0909-3613}\par}
\vspace{10mm}
{\large 28 lipca 2026\par}
\vfill
\begin{minipage}{0.89\textwidth}
\small
\textbf{Główny wynik.} Dokładna eliminacja ukrytych węzłów strict \(Z_{12}\)
generuje zależną od częstotliwości samoenergię i proces z pamięcią.
Statyczne dopełnienie Schura jest tylko granicą resolwentową, a nie dokładnym
generatorem całej zredukowanej dynamiki.

\medskip
\textbf{Granica.} Wynik nie generuje selektora, jednostek fizycznych,
mostu legacy--strict ani walidacji eksperymentalnej.

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
python3 audit_z12_sim.py
MPLCONFIGDIR=/tmp/matplotlib-puzzle \
python3 fin_nadsoliton_puzzle_atlas.py
\end{Verbatim}
Oczekiwany wynik audytu Z12: 15/15 testów. Punktacja atlasu jest heurystyką
wyszukiwania; twierdzenia operatorowe są weryfikowane osobnymi residualami.
\end{document}
"""

FIGURES = {
    "single_puzzle_domain_atlas.png": (
        "Heurystyczna macierz dopasowania 12 puzzli do 15 dziedzin. "
        "Punktacja służy wyłącznie do wyszukiwania analogii."
    ),
    "dynamic_schur_memory.png": (
        "Defekt składania zredukowanych propagatorów i błąd zastąpienia "
        "dynamicznej redukcji jednym statycznym generatorem Schura."
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
    )
    for filename in FIGURES:
        pattern = re.compile(
            r"!\[[^\]]*\]\(FIN_Nadsoliton_Puzzle_Atlas_Figures/"
            + re.escape(filename)
            + r"\)"
        )
        text = pattern.sub(f"\n\nFINFIG_{filename}\n\n", text)
    text = re.sub(r"`([^`]+)`", r"\1", text)
    text = re.sub(
        r"^(#{2,4})\s+\d+(?:\.\d+)*(?:\.)?\s+",
        r"\1 ",
        text,
        flags=re.MULTILINE,
    )
    # The Markdown source uses H2 as its top report level because H1 is the
    # human-readable title. The archival title page replaces H1, so promote
    # the remaining hierarchy by one level for proper chapter numbering.
    text = re.sub(
        r"^(#{2,4})(\s+)",
        lambda match: "#" * (len(match.group(1)) - 1) + match.group(2),
        text,
        flags=re.MULTILINE,
    )
    text = text.replace("# Confidence convention", "## Confidence convention", 1)
    body = converter.body_from_markdown(text)
    body = body.replace(
        r"\section{Confidence convention}",
        r"\chapter*{Confidence convention}"
        "\n"
        r"\addcontentsline{toc}{chapter}{Confidence convention}",
    )
    body = body.replace("Ż", r"\.Z").replace("ż", r"\.z")
    for filename, caption in FIGURES.items():
        token = f"FINFIG\\_{filename.replace('_', r'\_')}"
        path = f"FIN_Nadsoliton_Puzzle_Atlas_Figures/{filename}"
        body = body.replace(
            token,
            "\\begin{figure}[htbp]\n"
            "\\centering\n"
            f"\\includegraphics[width=0.96\\textwidth]{{{path}}}\n"
            f"\\caption{{{caption}}}\n"
            "\\end{figure}",
        )
    TARGET.write_text(PREAMBLE + "\n" + body + "\n" + POSTAMBLE, encoding="utf-8")
    print(TARGET)


if __name__ == "__main__":
    main()
