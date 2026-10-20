# Dryf entropijny rzeczywistego procesu usuwania przejść

Data: 2026-10-20. Samodzielny dowód warunkowy, nie formalizacja ani
niezależna recenzja. Podstawowe FIN HB, hipoteza HYP i stary kontrmodel CPL
pozostają niezmienione.

## Nowy wynik

Wyprowadziłem granicę wolnej trajektorii **istniejącego wcześniej procesu
porównawczego L-dagger na rzeczywistych licznościach**. Nie naprawiałem
po teście odrzuconej konstrukcji CPL. Dla stacjonarnego startu pi_D,N,
N podzielnych przez 12, dokładnego g_c i ustalonego horyzontu:

    dZ=[-Z^3/3+partial_z log p_P(Z)]ds+sqrt(2)dB,
    P(s)=P0, p_A(z)=1/2+exp(-z^2)/8, p_B=1-p_A.

Faza jest funkcją dotychczasowych liczności, nie dodanym ukrytym bitem.
Początkowa mieszanina faz wynika z niepostselekcjonowanego pi_D,N.
Zbieżność dotyczy całej wolnej ścieżki i fazy w J1, nie ścieżek szybkich.

## Co było brakującym krokiem

Samo usunięcie rzadkich krawędzi fazowych nie jest małą zmianą krytycznego
dryfu. Najpierw wyprowadziłem jego dokładny podpisany profil graniczny.
Następnie opłaciłem silniejszy krok: formy Dirichleta, rezolwenty i semigrupy
w zakresie wolnym. Ograniczona energia usuwa zależność od szybkich zmiennych
WEWNĄTRZ faz, a ich warunkowa masa p_P(z) pozostawia poprawkę logarytmiczną.

Dowód zachowuje poprawkę na prawdziwą warstwę wyjścia i nie zakłada wspólnej
luki spektralnej. Ścisłość trajektorii wynika z odwracalnego rozkładu na
martyngał w przód i wstecz, a nie z nieudowodnionego małego błędu dryfu
wzdłuż ścieżek. [Pełny wywód](PROOF_OR_OPEN_LEMMA.md).

## Co można teraz rzeczywiście porównać

Pierwotny HB i L-dagger mają tę samą skończoną miarę Gibbsa, a w granicy
tę samą miękką wariancję 2 i globalnie skompresowany dryf -z^3/3. Mimo to
ich zapisy w dwóch chwilach różnią się. Wyprowadziłem ograniczony kontrast
oraz jawne prawdopodobieństwo względne pomiędzy dwoma PRAWAMI GRANICZNYMI.
Ich KL jest co najwyżej S/30. Nie przenoszę tej granicy na entropię
skończonych ścieżek skokowych. [Konsekwencje](LIMIT_RECORD_CONSEQUENCES.md),
[rozdzielenie modeli](MODEL_SCOPE_TABLE.md).

## Kontrola i zakres

Przechodzi **41 ograniczonych kontroli**: współczynniki rzeczywistego grafu
K6,6, znaki kalibracji dualnej, mnożnik izometryczny, podpis odwróconego
martyngału, profil fazowy i małe przykłady pamięci. To nie 41 dowodów
asymptotycznych ani enumeracja dziedziny MT. Związano 17 wejść oraz
zachowano 182 wpisy dziewięciu manifestów. Dane spektralne i wersja A7
zgadzają się z zapisanymi naukowymi wejściami poprzedników.
[Wyniki](RESULTS.json), [odtworzenie](REPLAY.md).

**L-dagger nie zachowuje punktowo rodzimego dryfu**, więc wynik NIE jest
kontrprzykładem do KR pod jego rzeczywistymi założeniami. Pełne sprzężenie
z zadanym dryfem, cap CN i bezwarunkowe rodzime IBA pozostają otwarte.
Nie ma tu twierdzenia o starcie iid, fizycznym czasie, aparaturze lub
pochodzeniu maski przejść. F01, jednostki, geometria, walidacja, QM/GR/ToE
pozostają otwarte. Oryginały i AGENTS zachowano. [Handoff](HANDOFF.md).
