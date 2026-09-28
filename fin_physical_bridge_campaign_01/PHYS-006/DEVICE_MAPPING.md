# PHYS-006 — physical realization contract

## Wybrana platforma
Najmniej dodatkowych założeń wymaga programowalny elektroniczny/mixed-signal sampler kategoryczny: kontroler przechowuje 12-state labels, liczy siedem pól zbiorowych `M=sum xi_sigma`, odejmuje własny `xi_old`, tworzy 12 logitów
`ell_j=(g/N) xi_j·(M-xi_old)+vartheta 1[j=0]`, a skalibrowany fizyczny element losowy realizuje kategorię softmax. Zegar aktualizacji jest jawny.

## Minimalny prototyp
N=1,g=0: kontrola jednostajnego 12-state losowania i bias/readout.  
N=2: jedna etykieta może być kotwicą; druga ma realizować dokładnie `P(j|i)=softmax[(g/2)A7[i,j]]`. To wystarcza do pierwszego statycznego pair testu bez bariery/metastability.

## Literatura a projekcja
FB-MOSFET Potts p-bits (Advanced Materials 2026, PMCID PMC12994325) demonstrują fizyczną stochastyczność, one-hot multi-state sampling i Boltzmann-like systemy, ale publikacja nie jest demonstracją dokładnego q=12 FIN/A7 leave-one-out sampler. Autorzy pokazują głównie q=4 hardware i kontrolowane wielostanowe probabilistyczne jednostki.  
https://pmc.ncbi.nlm.nih.gov/articles/PMC12994325/

Coupled CMOS ring-oscillator Potts machine (arXiv:2504.11376) jest przede wszystkim architekturą optimization/graph-coloring. Samo znajdowanie niskiej energii nie certyfikuje Gibbs distribution ani wymaganych conditional rates.  
https://arxiv.org/abs/2504.11376

## Kalibracja konieczna przed jakimkolwiek EV
Trzeba zmierzyć `P(new label | old environment)` dla reprezentatywnych 12-logitowych wektorów, całkowity attempt rate, korelacje RNG, latency/jitter i rzeczywisty schedule. Target H pozostaje funkcjonałem samplera, dopóki osobno nie wykazano fizycznej relacji `beta_th H` do energii/temperatury urządzenia.
