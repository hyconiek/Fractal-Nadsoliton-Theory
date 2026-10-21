# Jeden bit krytyczny: budżet informacji to nie produkcja entropii

Data: 2026-10-21. Warunkowe badanie analityczne istniejącego rodzimego
modelu HB, bez nowego postulatu i bez niezależnej recenzji naukowej.

## Wynik

Zbadałem tylko bit 1{Z>0}, bez obserwowania pełnych etykiet. Udowodniłem,
że warunkowa niepewność następnego bitu przy znajomości pełnego bieżącego
stanu ma w granicy krytycznej postać

    c(t)=C_bit sqrt(t)+o(sqrt(t)), C_bit>0.

Stała jest jawną całką entropii progowanego Gaussa, pomnożoną przez
sqrt(2)rho(0). Jej dodatniość i skończoność mają opłacone oszacowania,
nie są wnioskiem z kilku punktów kwadratury. Bit nie jest traktowany jak
dwustanowy proces Markowa. To dolne ograniczenie rzeczywistej entropii
zapisu, nie jej rzekomo dokładny współczynnik.

## Co jest nowym wyprowadzeniem, a co ponownym użyciem

Przeniosłem tę wielkość z rzeczywistego skończonego procesu liczności
przez przyjętą silną granicę operatorów NT. Nieciągłość znaku, zera/tie,
moduł ciągłości entropii i zaokrąglanie dodatniego kroku czasu są opłacone.
Następnie wyprowadziłem krótkoczasową warstwę przy progu0, z kontrolą
dalekich ogonów dla nieograniczonego dryfu sześciennego.

Fano i argument przerzedzania siatki pochodzą z wcześniejszego ID; nie
przedstawiam ich jako nowych. Po zastosowaniu do tej nowej jednobitowej
wielkości źródło z ustalonym skończonym budżetem entropii H w czasie
krytycznym ma dodatnie ograniczenie błędu pełnego stacjonarnego sprzężenia.
Jawna funkcja tego koniecznego ograniczenia zachowuje się dla dużego H jak

    e_bit(H) ~ C_bit^2/[4H log H].

Nie jest to osiągalna optymalna funkcja rate-distortion ani projekt źródła.
Dla H=0 ten stacjonarny próg wynosi1/2. Twierdzenie wymaga całych historii
i zadanej mapy zegara; nie dotyczy dowolnego jednego skończonego okna.

## Rzeczywisty rodzic iid i informacja początkowa

Dla faktycznego, niepostselekcjonowanego późnego rodzica uzyskałem osobną
nierówność na skończone bloki. Płaci ona błąd TV przez liczebność alfabetu
zapisu i jawnie zawiera budżet I_N informacji początkowej źródła. Czasy
odczytu to tau_N+sqrt(N)jt; czynnika sqrt(N) nie pominięto.

Wynika stąd konkretny nowy test źródła: na N^2 przedziałach odczytu
(N^2+1 bitach) żaden deterministyczny model z najwyżej 12^N możliwymi
pełnymi stanami początkowymi
nie uzyska zanikającego średniego błędu bitowego. Dotyczy to nawet dowolnych
ustalonych odczytów zależnych od czasu; nie wymaga stacjonarności źródła.
To wykluczenie zamkniętej klasy o ograniczonej informacji początkowej,
nie dowód fundamentalnego indeterminizmu lub skończoności nadsolitonu.

Zerowa entropia KS nie wystarcza, aby wykluczyć odtworzenie długiego
skończonego zapisu. Dokładny cykl de Bruijna daje kontrolę tego błędu:
ma zerową stopę entropii, lecz początkowa losowa faza płaci skończone
niezależne bity. Jeśli losuje się całą historię, budżet obejmuje także
wybór tej historii, nie tylko fazę. Zachowuje to wcześniejszy wynik FH.

## Znaczenie i granice

Zanik strzałki czasu w sensie KL, opłacony w TR, nie usuwa losowości i
informacyjnego kosztu odczytów. To różne wielkości. Nowy wynik precyzuje
wymagania wobec proponowanego źródła ewolucji; nie dostarcza jeszcze tego
źródła, fizycznego zegara, aparatury ani teorii obserwatora.

Nie zakładam, że nadsoliton ma skończoną entropię KS. Nie ogłaszam
kurczącego się kroku t_N->0 dla dowolnej wspólnej granicy N ani praktycznego
N0. Udzielone granice i metryki są jawnie rozdzielone. F01, jednostki,
geometria, kalibracja/EV i QM/GR/ToE pozostają otwarte.

[Pełny dowód](PROOF_OR_OPEN_LEMMA.md),
[skończony zapis i budżet początkowy](FINITE_RECORD_AND_INITIAL_BUDGET.md),
[wyniki kontroli](RESULTS.json), [odtworzenie](REPLAY.md).
AGENTS, wszystkie wcześniejsze źródła i istniejące zmiany zachowano.
Przechodzi 36 ograniczonych kontroli, związano 17 wejść i zachowano 73
wpisy pięciu manifestów poprzedników. Nie są to niezależne dowody
ani numeryczne certyfikaty krytycznej domeny przy małym N.
