# Scalenie metodologicznie poprawnych wyników R7O2

Data: 2026-09-20. Źródło: `FIN_R7O2_CONTINUATION_HANDOFF_20260920`.
Oryginalny pakiet pozostawiono niezmieniony. Przyjęcie dotyczy poniższych
zakresów i nowych certyfikatów zastępczych, nie wszystkich deklaracji z historii.

## Przyjęte wyniki

1. **Metoda zależnego od wspólnych parametrów momentu drugiego.** Zachowuje
   zależności między wagami w `(r,s,t,y)` i ich wspólną normalizacją. Dla
   dokładnej wymiernej bazy B rzędu 3 i stałego środka c dowodzi dodatniości
   `K=(67/250)BᵀB−E[(BᵀF−c)(BᵀF−c)ᵀ]`. Ponieważ kowariancja jest nie
   większa od tego momentu drugiego, dodatniość K wystarcza do ograniczenia
   drugiej wartości własnej M4 przez 67/250. Jest to metoda matematyczna
   dla zadanej klasy modeli, a nie źródło fizycznego gainu.
2. **Pełny replay 7 340 aktywnych bezpiecznych liści.** Zachowano zapisane
   wymierne bazy B, niezależnie sprawdzono ich rangę i zapisano nowe jawne
   wymierne środki c. Wszystkie nowe certyfikaty przeszły kontrolę. Nie
   wykorzystano zapisanych floatów d1/d2/d3 jako dowodów dodatniości.
3. **Dokładna geometria 3 460 przetworzonych rodziców.** Sprawdzono podziały
   wszystkich aktywnych drzew, a nie tylko sumę ich objętości. Naprawy
   rodziców 332, 338, 357 i 363 zastępują wcześniejsze drzewa i nie są
   doliczane ponownie.
4. **Przyjęte jest częściowe pokrycie, nie globalne Target P.** Całkowicie
   zamknięto 3 406 rodziców. W 54 częściowo zamkniętych rodzicach pozostaje
   171 nierozstrzygniętych terminali. Dalszych 1 972 rodziców nie przetwarzano.

## Dwa poprawne sposoby rozliczania pozostałej dziedziny

Wszystkie poniższe wielkości odnoszą się do objętości zwartej dziedziny
współrzędnych, nie do prawdopodobieństwa fizycznego lub pewności twierdzenia.

| Wielkość | Ułamek objętości zwartej dziedziny |
|---|---:|
| Całkowicie zamknięte rodzice R7O2 | 0.23991936357483346 |
| Całe 54 rodzice oczekujące na naprawę | 0.00010424154100882935 |
| Tylko ich nierozstrzygnięte terminale | 0.000006331724266361647 |
| Nieprzetworzone 1 972 rodziców | 0.1230813421248372 |
| Ostrożny residual liczony całymi rodzicami | 0.12318558366584603 |
| Residual liczony faktycznymi terminalami | około 0.12308767384910356 |

Po uwzględnieniu wcześniej przyjętego pokrycia R7N, dowód obejmuje około
**87,6912% zwartej dziedziny**, a około **12,3088%** pozostaje nierozstrzygnięte.
Nie ogłasza się globalnego `lambda2(M4)<=67/250`. Ostrzejszy próg sigma,
pełny model X7 i źródła fizyczne również nie zostały domknięte.

## Korekty metodologiczne

- W źródłowych liściach brakowało `center_c`. Zamiast udawać dokładne odtworzenie
  nieobecnego świadka, wydano jawne certyfikaty zastępcze na **tych samych
  komórkach**, ze stałą zapisaną bazą i nowym zapisanym wymiernym środkiem.
- Sprawdzono dokładny podział rodziców na dzieci. Równość objętości sama nie
  wyklucza nakładania i luk.
- Sprawdzono wszystkie 3 100 wpisów manifestu oraz hashe aktywnych plików.
  Wszystkie 5 432 wejściowe komórki są identyczne z zaakceptowanym R7N.
- Ścieżki `/mnt/data/...` przekierowano **w pamięci procesu** do dostępnego
  poprzednika. Archiwalnego kodu i danych nie edytowano.
- Do pełnego replay użyto przedziałów binary64 z jawnym zaokrąglaniem na
  zewnątrz przez `nextafter`. Końce udostępniane pozostałemu kodowi są dokładnymi
  liczbami wymiernymi, więc obliczenia środków i promieni Taylora nie tracą
  kierunku zaokrąglania. To nie są zwykłe obliczenia float z tolerancją.
  Operacje niefinitywne i dzielenie przez przedział zawierający zero są odrzucane.
- Pierwszych 20 liści sprawdzono również wcześniejszym wymiernym backendem
  z zaokrąglaniem na siatkę 10^-12. Nowy interfejs pozwala odtwarzać dowód
  z już zapisaną bazą i środkiem, bez ponownego generowania propozycji.
  Nie twierdzimy, że cały replay wykonano
  dwukrotnie tymi dwoma backendami.
- Nie przyjęto eksploracyjnego residualu około 2,26% z wcześniejszej gałęzi R7O.
  Nie rozpoczęto obliczeń nad nieprzetworzonymi rodzicami ani nowych napraw.

## Uzasadnienie przedziałowej metody

Tożsamość
`E[(Z-c)(Z-c)ᵀ]−Cov(Z)=(E[Z]−c)(E[Z]−c)ᵀ >= 0`
obowiązuje dla każdego stałego c. B oraz c pozostają stałe podczas różniczkowania.
Weryfikowany model ma siedem zagregowanych wag
`1, 2st³, A²t⁴, 2A²st, 2At²y, 2Ast²y/z, 2Ast²zy`,
gdzie `A=sqrt(r)` i `z=t^sqrt(3)`. Przedziałowe wartości oraz pochodne z
obejmują tę samą rzeczywistą funkcję. Kotwica 1 zapewnia dodatnią normalizację.

Drugorzędowe jety propagują pochodne iloczynu i odwrotności. Rozwinięcie Taylora
w środku zewnętrznego pudełka `(A,1-s,1-t,y)` z ograniczeniem hesjanu na całym
pudełku obejmuje każdy wpis momentu drugiego. Dodatniość minorów Sylvestera
lub dolnego ograniczenia Gershgorina daje dodatnią formę na trójwymiarowej
podprzestrzeni obrazu B. Z zasady min–max wynika żądana kontrola lambda2(M4).
Ani optymalność numerycznej bazy, ani jej idealna ortonormalność nie są potrzebne.

Backend jest sprawdzany względem dokładnej arytmetyki wymiernej oraz przez
negatywne testy pokrycia. Jest to dowód analityczny wsparty rachunkiem
przedziałowym, nie formalizacja w asystencie dowodowym.

## Artefakty i odtwarzanie

- [Rejestr geometrii i pochodzenia](registry.json).
- [Aktywne liście źródłowe](safe_leaves.json).
- [Nowe świadki i pełny replay](leaf_replay.json).
- [Weryfikacja końcowa](verification.json).
- [Checker](review.py), [backend przedziałowy](intervals_fast.py),
  [testy kontrolne](test_review.py).

```bash
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7o2_review/review.py registry
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7o2_review/review.py replay --workers 4
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7o2_review/finalize.py
```

Ostatni krok sprawdza zapisane świadectwa i integralność; nie zastępuje
przeliczenia wszystkich liści. Wszystkie wyniki trafiają tylko do tego katalogu.
Źródłowy ZIP poprzednika nie jest wymagany do rachunku, jeśli dostępne jest
zweryfikowane rozpakowane drzewo; nie odtworzono ani nie przywrócono usuniętych ZIP-ów.
