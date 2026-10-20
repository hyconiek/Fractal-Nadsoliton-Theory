# Wyprowadzono krytyczny proces — nie tylko jego formalne współczynniki

Data: 2026-10-20. Wynik warunkowy w upoważnionym modelu hipotetycznym.
Nie zmieniono podstawowego FIN ani jego zegara HB.

## Nowy wynik

Z rzeczywistych skoków SQRT otrzymano granicę całych trajektorii parzystości:

    dZ=12(chi lambda Z-Z^3/3)ds+sqrt(24)dB,
    Z0~rho_lambda.

Dowód obejmuje dziesięć szybkich kierunków, poprawkę leave-one-out,
odrzucane wyjścia porównania odbitego, kontrolę maksimum i martyngałów oraz
ich wspólną filtrację. Szum nie został dopisany. Czynnik 12 pochodzi
z zadanych nowych szybkości, nie z dopasowania czasu do HB.

Wynik dotyczy najpierw dokładnego lokalnego startu stacjonarnego. Następnie
przeniesiono go na rzeczywisty, niepostselekcjonowany późny start iid,
korzystając z wcześniej udowodnionego okna lokalnego. Osobno opłacono
skończone źródło entropijne z L_N rzędu N^3. To słaba zbieżność ścieżek,
nie zgodność wszystkich idealnych znaczników czasu ani posteriorów w TV.

## Ważna zmiana dalszego kierunku

Granica SQRT powstała z rodziny z dodatkowymi jednostajnymi współrzędnymi.
Nie należy zastępować nią oryginalnego źródła hipotezy na N oddziałujących
współrzędnych. Dla tego oryginalnego źródła wyprowadziłem już długie
przebywanie w sektorze bez takich dodatków. Nie udowodniłem jego lokalnego
mieszania ani tej samej granicy krytycznej.

Następny konkretny brak to kontrola jego rzeczywistego dodatniego skalowania
v_N, określonego jednoznacznym równaniem. [Zlecenie](NEXT_INTRINSIC_SCALING_ORDER.md)
nie wymyśla kolejnego prawa; bada to już upoważnione, bez zastępowania
globalnego problemu nowym problemem warunkowym na lokalnej dziedzinie.

Źródła, AGENTS i historyczne plany zachowano. Kontrole dokładne są małymi
testami rachunku, nie dowodem asymptotycznym ani niezależną recenzją.
Przechodzą 42 ograniczone kontrole wymierne; związano 20 plików wejściowych
i dwanaście manifestów poprzedników. W następnym kroku wystarczy badać dolne
ograniczenie aktywności na krawędziach: osobny lemat opłaca już ujemne widmo
w czasie dyskretnym bez zakładania stałej leniwości lub zbędnego górnego
ograniczenia.
Naturalne źródło hipotezy, A7 i inicjalizacji, fizyczny czas, jednostki,
geometria, walidacja, QM/GR i ToE nadal nie są wyprowadzone; F01 pozostaje
otwarte. Szczegóły: [dowody](PROOF_OR_OPEN_LEMMA.md), [wyniki](RESULTS.json).
