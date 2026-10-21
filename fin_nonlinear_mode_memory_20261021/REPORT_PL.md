# Wyprowadzony wpływ zwrotny pola: jawna pamięć nieliniowa

Data: 2026-10-21. Konstruktywny warunkowy wynik dla istniejącego modelu
polowego F_RP; nie nowy postulat, nie rodzime HB i nie niezależna recenzja.

## Wynik

Z pełnych sprzężonych równań wyprowadziłem równanie amplitudy modu
związanego, z pamięcią po przyczynowym wyeliminowaniu drugiego pola:

    z_tt+omega^2 z+(18/35)g4 z^3
       -72 z(t) integral_0^t k(t-s)z(s)^2 ds=R_epsilon(t).

Współczynnik72 jest skutkiem dwóch rzeczywistych sił sprzężenia6 i12,
a nie dobraną stałą ruchliwości. Zachowano podpisane g4,g6, początkowy
kink i brak niezależnej fali wejściowej. Składnik kwintyczny g6 nie został
ustawiony na zero: jest opłacony w reszcie.

Jednocześnie skonstruowałem skalarne równanie Volterry z R=0 i dowiodłem
błędu jego przewidywania O(epsilon^5(1+t)^4) na przedziale spełniającym
epsilon^2(1+T)^3<=kappa. To rzeczywista kontrolowana redukcja tego samego
rozwiązania, nie tylko formalna zamiana pełnego pola na oscillator.

## Co zawiera jądro pamięci

Jądro obliczono dla właściwego profilu sprzężenia i normalizacji widmowej:

- mod kształtu daje trwałą oscylację o częstości sqrt3;
- kontinuum daje ogon sin(2t+3pi/4)/(16sqrt(pi)t^(3/2));
- statyczna podatność kanału wynosi1/70.

Pełne jądro nie jest całkowalne bezwzględnie. Kontinuum jest całkowalne,
ale jego pierwszy moment bezwzględny jest rozbieżny. Nie wolno więc użyć
nieopłaconego argumentu o „krótkiej pamięci” i od razu dostać tarcia Markowa.
Statyczna podatność nie usuwa oscylacji ani pamięci początkowego włączenia.

Odpowiedź liniowego kanału jest pasywna, zgodnie z energią równania falowego.
Jej dokładne odtworzenie przez skończenie wiele liniowych stanów o stałych
współczynnikach jest niemożliwe ze względu na niezerowy algebraiczny ogon.
To ograniczenie dokładnej liniowej realizacji, nie zakaz nieliniowego
ukrytego mechanizmu ani kontrolowanej aproksymacji na skończonym oknie.

## Znaczenie fundamentalne

To przykład rzeczywistego nieliniowego sprzężenia zwrotnego wyprowadzonego
z już zadanej ewolucji, a nie dopisanego jako nowa mikroskopowa reguła.
Metoda projekcji/pamięci jest znana; nowy wkład to jawny profil FIN,
rozkład widmowy pamięci i opłacony błąd nielokalnego predyktora.

Nie powstał niezależny szum Gaussowski. Przy aktualnych zerowych danych
pola eliminowanego jego wkład jednorodny jest zerowy. Dodatnia gęstość
Gamma lub pasywność same nie wybierają rozkładu początkowego ani prawa
stochastycznego. Nie twierdzę też, że redukcja działa na skali epsilon^-2
lub dowodzi całkowitego zaniku modu.

Badany model nadal ma zadaną akcję legacy, współczynniki i geometrię1+1.
Jego aktywny sektor nie rozpoznaje niekomutującego bloku spektatorów.
Nie jest to źródło naturalnego A7/HB, czasu, jednostek, cząstek lub QM/GR.
F01 pozostaje otwarte. Kink K(x) i pamięć k(t) nie zastępują przestrzennego
K_legacy/K_strict. Nie dodano warstwy informacji pod nadsolitonem.

[Pełny wywód](PROOF_OR_OPEN_LEMMA.md), [redukcja](NONLINEAR_REDUCTION.md),
[jawne jądro](SPECTRAL_MEMORY_KERNEL.md), [kontrole](RESULTS.json),
[replay](REPLAY.md). Źródła, manifesty, plany i AGENTS zachowano.
Przechodzą 34 ograniczone kontrole dokładnego rachunku; związano 13 wejść
i zachowano 63 wpisy manifestów. Nie są to niezależne dowody ani
komputerowa weryfikacja całego twierdzenia o PDE i długoczasowym jądrze.
