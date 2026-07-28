#!/usr/bin/env python3
"""Convert the Programs 255--266 research report to archival LaTeX."""

from pathlib import Path
import re

import fin_monograph_to_latex as converter


SOURCE = Path("FIN_Programy_255_266_Raport_Badawczy_PL.md")
TARGET = Path("FIN_Programy_255_266_Raport_Badawczy_PL.tex")

PREAMBLE = r"""\documentclass[11pt,a4paper,openany]{report}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
\usepackage{microtype}
\usepackage{amsmath,amssymb,amsthm,mathtools,bm}
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
  pdftitle={FIN Programs 255--266: Operator Memory Measures, Identifiability, Current Tomography and Adaptive-Analogy Falsification},
  pdfauthor={Krzysztof Żuchowski},
  pdfsubject={Executed FIN research Programs 255--266 and recommended Programs 267--280},
  pdfkeywords={FIN, operator-valued Stieltjes functions, Schur complement, memory kernels, process information, current tomography, identifiability, reservoir computing}
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small Programy FIN 255--266}
\fancyhead[R]{\small Release 10.24}
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
\vspace*{11mm}
{\Large\bfseries FIN --- Release 10.24\par}
\vspace{10mm}
{\Huge\bfseries Programy badawcze 255--266\par}
\vspace{5mm}
{\LARGE\bfseries Miary operatorowe pamięci, identyfikowalność,\\
tomografia prądów i falsyfikacja analogii adaptacyjnych\par}
\vspace{16mm}
{\Large Krzysztof \.Zuchowski\par}
\vspace{2mm}
{\normalsize Independent Researcher --- Fractal Information Theory Project\par}
\vspace{2mm}
{\normalsize ORCID: \href{https://orcid.org/0009-0002-0909-3613}{0009-0002-0909-3613}\par}
\vspace{9mm}
{\large 28 lipca 2026\par}
\vfill
\begin{minipage}{0.9\textwidth}
\small
\textbf{Centralny wynik.}
Dla każdego nietrywialnego kontekstu strict \(Z_{12}\):
\[
\Sigma_E(z)=\int_{(0,\infty)}\frac{d\mathsf M_E(\mu)}{z+\mu},
\qquad d\mathsf M_E(\mu)\succeq0.
\]
Pamięć określa minimalną klasę realizacji wejście--wyjście, a nie unikalny
sektor mikroskopowy.

\medskip
\textbf{Granica.}
Raport nie eksportuje selektora, skali fizycznej, pełnego mostu
legacy--strict, prawa adaptacji ani walidacji laboratoryjnej.

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
\chapter*{Certyfikat reprodukcji}
\addcontentsline{toc}{chapter}{Certyfikat reprodukcji}
\begin{Verbatim}[fontsize=\small]
MPLCONFIGDIR=/tmp/matplotlib-p255-266 \
python3 fin_programs_255_266.py

MPLCONFIGDIR=/tmp/matplotlib-p255-266 \
python3 -m unittest -v test_fin_programs_255_266.py
\end{Verbatim}
Oczekiwany wynik testów: 15 passed, 0 failed.  Dane P262, P264, P265 i P266
są syntetyczne; nie są rekordem laboratoryjnym.
\end{document}
"""

FIGURES = {
    "p255_p258_stieltjes_chiral.png": (
        "Granice funkcji Stieltjesa oraz norma chiralnej podatności pamięci wraz z udowodnionym ograniczeniem."
    ),
    "p259_p260_rg_information_ledger.png": (
        "Ściśle malejący przepływ kontekstów i wieloczasowy bilans rozróżnialności."
    ),
    "p261_p262_bridge_fingerprint.png": (
        "Obstruction amplitude-only completion oraz syntetyczny audyt fingerprint-first/calibration-second."
    ),
    "p264_false_positive_atlas.png": (
        "Częstości przejścia pięciu testów w zamrożonych zespołach zerowych."
    ),
    "p266_reservoir_benchmark.png": (
        "Rezerwuar FIN na tle 120 kontroli dopasowanych wymiarem lub widmem."
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
    text = text.replace("## Konwencja pewności", "## Confidence convention", 1)
    for filename in FIGURES:
        pattern = re.compile(
            r"!\[[^\]]*\]\(FIN_Programs_255_266_Figures/"
            + re.escape(filename)
            + r"\)"
        )
        text = pattern.sub(f"\n\nFINFIG_{filename}\n\n", text)
    text = re.sub(r"`([^`]+)`", r"\1", text)
    body = converter.body_from_markdown(text)
    body = body.replace(
        r"\section{Confidence convention}",
        r"\chapter*{Konwencja pewności}"
        "\n"
        r"\addcontentsline{toc}{chapter}{Konwencja pewności}",
    )
    body = body.replace(
        r"\begingroup\footnotesize",
        r"\begingroup\fontsize{7}{8.4}\selectfont",
    )
    body = body.replace(
        r"\setlength{\tabcolsep}{3.5pt}",
        r"\setlength{\tabcolsep}{2.5pt}",
    )
    body = body.replace("Ż", r"\.Z").replace("ż", r"\.z")
    for filename, caption in FIGURES.items():
        token = f"FINFIG\\_{filename.replace('_', r'\_')}"
        path = f"FIN_Programs_255_266_Figures/{filename}"
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
