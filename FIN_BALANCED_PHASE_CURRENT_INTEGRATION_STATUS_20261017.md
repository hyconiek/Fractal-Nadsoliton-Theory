# Właściwy pakiet PC: potwierdzenie integracji

Data kontroli: 2026-10-17. HEAD podczas kontroli: `dc51caba`.
Status: **EXISTING_PC_INTEGRATION_REVALIDATED; NO_NEW_THEOREM**.

Wskazanym przez użytkownika obiektem jest
[fin_balanced_phase_current_campaign_20261016](fin_balanced_phase_current_campaign_20261016/HANDOFF.md),
nie wcześniejszy `fin_balanced_phase_cut_research_20261016` ani kampania SC.
W aktualnym repozytorium właściwa kampania i jej
[osobna recenzja](fin_balanced_phase_current_review_20261017/INTAKE_20261017.md)
są już śledzone. Wprowadził je istniejący commit `48add8fc`; późniejsze
aktualizacje recenzji obejmuje `68e01759`. Obecna kontrola nie tworzy commitu
i nie przypisuje sobie autorstwa tych wcześniejszych integracji.

## Co potwierdzono ponownie

- Manifest źródła: **26/26**, pełne pokrycie z zadeklarowanymi wyłączeniami;
  wyłączony sprawdzacz i metadane mają osobno sprawdzone hashe.
- Manifest recenzji: **13/13**; kontrola JSON, składni, odnośników i zależności
  przechodzi w trybie tylko do odczytu.
- Weryfikator recenzji z replayem: **17/17 ograniczonych kontroli**.
  Oryginalny skrypt naukowy uruchomiony wyłącznie na kopii w `/tmp` odtwarza
  `verification.json` identycznie bajtowo. Są to kontrole spójności i małych
  przykładów, nie siedemnaście dowodów asymptotycznych.
- Z ośmiu pierwotnych odwołań wejściowych siedem zgadza się z lokalnymi
  plikami. Historyczny obraz AGENTS pozostaje niedostępny; nie utożsamiono
  go z aktualnym plikiem. ZIP użyty w replayu jest lokalnie utworzoną kopią
  kontrolną, nie odzyskanym oryginalnym uploadem.

Pole `entry_HEAD` w wyniku istniejącego weryfikatora zachowuje historyczne
`4aa415bc`; nie oznacza HEAD tej ponownej kontroli. Pełny bieżący wynik
zapisano w `/tmp/fin-pc-corrected-target-revalidation-20261017.json`.
Trwałe pierwotne dowody kontroli pozostają w
[verification.json](fin_balanced_phase_current_review_20261017/verification.json).

## Przyjęty wynik i granica

Przyjęte pozostają lokalne porównanie wag Gibbsa oraz nierówność Halla dla
zbiorów spoza warstwy wyjścia, które mają wybór dziewięciu bezpiecznych
następników o nakładaniu najwyżej cztery. Dokładny zapas wynosi
`729/2000 - 357/1000 = 3/400`, pomnożony przez `N W(T)`.

**Nie przyjęto pełnego sprzężenia sąsiednich włókien.** Nie wykazano tej
własności nakładania dla wszystkich gęstych zbiorów ani kontroli wszystkich
stanów warstwy wyjścia. PC-003 pozostaje NOT_ACTIVATED. Pełny generator,
sklejenie na całym D_N i kontrrodzina obalająca IBA nie wynikają z kampanii PC.
Brak też nowego fizycznego źródła A7, aktywności, zegara lub jednostek.

## Następny krok bez powtarzania wykonanej pracy

Zachować [istniejący plan po PC](FIN_PC_PHYSICAL_COMPATIBILITY_NEXT_CHAIN_20261017.md)
jako historyczny punkt wyjścia. Katalog
[późniejszej kampanii DH](fin_dense_orbit_fibre_transport_campaign_20261017/HANDOFF.md)
już istnieje: nie zlecać jej ponownie jako niepodjętej pracy. Jej handoff
zgłasza wynik dla szczególnej klasy dopełnień, ale nadal pozostawia pełną
nierówność Halla otwartą. Ta wzmianka identyfikuje aktualne materiały;
niniejsza rewalidacja PC nie jest odrębnym przyjęciem twierdzeń DH.

Dalszy łańcuch matematyczny: audyt rzeczywiście dostarczonego DH -> brakująca
nierówność dla pozostałych gęstych sum orbit na pełnym H -> sprzężenie
wierszy i kolumn -> oddzielne sklejenie z resztą D_N. Każdy krok wymaga
własnego dowodu albo jawnego wyniku OPEN; nie wystarczą lokalne stopnie
wierzchołków lub kolejne stałe dla rzadkich zbiorów. Most fizyczny pozostaje
osobnym zadaniem wymagającym rzeczywistego, niekołowego obiektu źródłowego.

Nie zmieniono AGENTS, źródeł PC, recenzji ani zastanych prac.
Nie dopisano drugi raz już przyjętych wyników.

## Odtworzenie kontroli i związanie wejść

Z katalogu głównego repozytorium:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_balanced_phase_current_review_20261017/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_balanced_phase_current_review_20261017/verify_intake.py --replay --output /tmp/fin-pc-corrected-target-revalidation-20261017.json
```

SHA-256 podczas kontroli:

```text
ffb62af2671f62aa582c688a90edd998e4a9621b29e35fd240e1863038c9e0d9  AGENTS.md
d941ae114972a276fd2cc06dd78526154848bc0c4f67b97d3df3f14cdd49dfc1  fin_balanced_phase_current_campaign_20261016/MANIFEST.sha256
5e6ea3ad58c939c728ccecbbd17065a63e8e93f7d0fc1a894105a5fcfaa7da30  fin_balanced_phase_current_review_20261017/MANIFEST.sha256
```
