# A1 Minimal Action Ansatz

Status: `A1_SPEC_READY_SINGLE_NADSOLITON_GUIDED`
As of: `2026-03-06`

## Cel

Zapisac najmniejszy recenzencko-obronny ansatz dzialania, ktory:
- jest lokalny,
- ma maksymalnie drugi rzad pochodnych w warstwie bazowej,
- dopuszcza tlo supersolitonowe jako przyszly target etapu `A2`,
- zostawia miejsce na sektor gauge i grawitacyjny,
- nie udaje jeszcze pelnego SM+GR closure.

## Ontologiczna wskazowka programu

Ten tor konstrukcyjny jest prowadzony pod zalozeniem roboczym:
- warstwa pierwotna ma charakter informacyjny,
- fundamentalnym obiektem konstrukcyjnym jest jeden nadsoliton,
- pozostale sektory fizycznego opisu moga byc emergentne albo efektywne wzgledem tej jednej struktury.

To nie jest jeszcze theorem-level closure.
To jest wskazowka konstrukcyjna dla `A1..A10`, a nie recenzencki dowod, ze single-nadsoliton ontology jest juz domknieta.

## Minimalna zawartosc pol

1. `Psi^A`
   - jedyne pole traktowane na tym etapie jako ontologicznie fundamentalne,
   - wieloskladowy nosnik struktury nadsolitonowej.
2. `Phi`
   - efektywny parametr porzadku / vacuum-order parameter,
   - dopuszczony jako pole pomocnicze lub funkcjonal emergentny od `Psi`,
   - nie jest traktowany jako drugi wspolfundamentalny byt.
3. `A_mu^I`
   - efektywna rezerwa dla przyszlego sektora cechowania,
   - nie jest jeszcze traktowana jako rownorzedne pole bazowe,
   - grupa gauge nie jest jeszcze ustalona.
4. `g_{mu nu}`
   - efektywna warstwa geometryczno-grawitacyjna,
   - dopuszczona jako przyszly sektor emergentny lub sprzezony,
   - nie jest jeszcze rownorzednym fundamentem ontologicznym.

Sektor fermionowy nie jest uznany za domkniety w `A1`.
Zostaje odlozony do pozniejszego rozgalezienia `A5`.

## Ansatzu nie wolno mylic z finalnym lagranzianem

`A1` nie twierdzi, ze finalna teoria ma juz postac ponizej.
To jest najmniejsza klasa kandydatow, od ktorej sensownie zaczac dalsza rekonstrukcje przy zachowaniu ontologii jednego nadsolitonu.

## Bazowy ansatz

Najpierw zapisujemy strukture warstwowa:

\[
\mathcal{L}_{A1} =
\mathcal{L}^{\mathrm{base}}_{\Psi}
+ \mathcal{L}^{\mathrm{eff}}_{\Phi}[\Psi]
+ \mathcal{L}^{\mathrm{eff}}_{A}[\Psi]
+ \mathcal{L}^{\mathrm{eff}}_{g}[\Psi]
+ \Delta \mathcal{L}_{\mathrm{EFT}}.
\]

Minimalna realizacja tej klasy moze byc zapisana jako:

\[
\mathcal{L}_{A1} =
\frac12 G_{AB}(\Psi,\Phi) D_\mu \Psi^A D^\mu \Psi^B
+ \frac12 \partial_\mu \Phi\, \partial^\mu \Phi
- V_{\mathrm{eff}}(\Psi,\Phi)
- \frac14 Z_{IJ}(\Psi,\Phi) F^I_{\mu\nu} F^{J\mu\nu}
+ \frac{M^2_{\mathrm{eff}}(\Psi,\Phi)}{2} R
- \Lambda_{\mathrm{eff}}(\Psi,\Phi)
+ \Delta \mathcal{L}_{\mathrm{EFT}}.
\]

Gdzie:
- `Psi^A` pozostaje jedynym ontologicznie fundamentalnym polem,
- `Phi` jest traktowane jako efektywna warstwa porzadku zwiazana z `Psi`,
- `A_mu^I` i `g_{mu nu}` sa dopuszczone jako przyszle warstwy emergentne lub efektywne,
- `D_mu` jest zwykla pochodna w wariancie bez aktywnego gauge sector albo pochodna kowariantna po aktywacji `A_mu^I`,
- `F^I_{mu nu}` jest krzywizna polaczenia,
- `Delta L_EFT` reprezentuje wylacznie jawnie oznaczone wyzsze operatory tlumione skala `Lambda`.

## Twarde restrykcje na A1

1. lokalnosc,
2. Lorentz-covariance,
3. brak niejawnych nielokalnych kernel terms w samym ansatzu bazowym,
4. bazowo najwyzej drugi rzad pochodnych,
5. kazdy wyzszy operator musi byc jawnie wpisany jako EFT correction,
6. `Psi` pozostaje jedynym polem traktowanym jako ontologicznie fundamentalne,
7. `Phi`, `A_mu^I`, `g_{mu nu}` nie moga byc na tym etapie sprzedawane jako wspolfundamentalne byty,
8. zaden wybor grupy gauge nie moze byc jeszcze sprzedawany jako wyprowadzony,
9. brak roszczenia, ze fermiony Diraca sa juz uzyskane.

## Obowiazki recenzenckie odlozone poza A1

- istnienie rzeczywistego tla supersolitonowego jako rozwiazania E-L,
- identyfikacja fizycznego kernela fluktuacji,
- kontrola zero modes / gauge fixing,
- RG emergence,
- spinor/gamma emergence lub spin-bundle extension,
- anomaly closure,
- GR-limit equivalence.

## Kryterium zaliczenia A1

`A1` jest zaliczone wtedy i tylko wtedy, gdy:
- klasa dzialan jest jawnie zapisana,
- zalozenia sa jawnie oznaczone jako `physical` albo `technical`,
- ontologiczna rola `Psi` jako jedynego fundamentu jest jawnie odrozniona od warstw efektywnych,
- nie ma falszywego claimu unikalnosci albo full closure,
- dokument daje czysty punkt wejscia do `A2` i `A3`.
