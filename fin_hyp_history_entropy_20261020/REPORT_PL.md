# Zbadana hipoteza źródła: FIN-HYP-HIST-ME-001

Data: 2026-10-20. **Nowa hipoteza za zgodą autora, nie strict FIN i nie
nowo odkryta ogólna metoda maksymalizacji entropii.**

## Co rzeczywiście powstało

Zdefiniowano i zbadano konkretny postulat: spośród stacjonarnych historii
na zadanym nośniku, o zadanej mierze stanów i dozwolonych zmianach,
wybrać prawo o największej entropii Shannona na wewnętrzny krok.
W zamrożonym modelu jeden nieoznaczony pusty krok liczy się raz.

**Postulat rzeczywiście wybiera jedno prawo.** Nie założono osobno jego
markowskości ani odwracalności: wynikają warunkowo z maksimum i symetrii
grafu. Udowodniono istnienie, dodatniość, jednoznaczność, pełne prawo
historii i dokładny certyfikat optymalności przez KL.

To jest warunkowy model źródła dynamiki, a nie tylko odwrócenie wcześniej
założonego Q. Nadal jednak jego nośnik, miara, wsparcie, wewnętrzne kroki
i sama zasada maksimum są jawnymi przesłankami, nie wyprowadzeniem z
dotychczasowej ontologii nadsolitonu.

## Ważny wynik przeciwny oczekiwaniu

**Zamrożony postulat nie odtwarza rodzimego HB dla niejednorodnego N2 FIN,
nawet po zmianie zegara.** Wybiera inne proporcje następnych stanów.
W dokładnym przykładzie czterech etykiet jedno prawdopodobieństwo wynosi
1/5 zamiast 15/71, różnica -4/355. Ogólny dowód obejmuje również odziedziczony
dodatni, niejednorodny rozkład q12/A7, bez dopasowywania cyfr A7.

Zachowano ten niekorzystny wynik. Nie zmieniono liczby pustych kroków,
nośnika ani funkcjonału, żeby odzyskać HB. Nie oznacza to empirycznego
obalenia hipotezy: nie dostarczono pomiarów. Odrzucono jej konkretny most
do istniejącego modelu HB, nie samą możliwość istnienia innej fizyki.

## Dlaczego zwykłe hasło „maksymalna entropia” nie wystarcza

| Dokładnie określona wersja | Wynik dla dwóch kopii |
|---|---|
| Zamrożony nośnik oznaczonych stanów, jeden pusty krok | Unikalne prawo, na ogół nie HB nawet modulo zegar |
| Maksymalizacja bezpośrednio na samych różnicach etykiet | Inne zadanie; daje iid p i gałąź HB |
| Dwa odróżniane puste kanały prób | Nieprzyjęta kontrola; odtwarza HB tylko w specjalnym N2 |

Wyprowadzono dokładnie, skąd bierze się różnica: przy zmianie różnicy
istnieją dwa mikroskopowe sposoby, a przy pustym kroku w zamrożonym modelu
jeden. Po redukcji pozostaje dodatkowy człon aktywności razy log 2.
Pominięcie go zmienia postulat. Dodatkowy wymierny przykład ośmiu stanów
pokazuje, że nawet odzyskanie HB po dodaniu oznaczonych pustych kanałów
w N2 nie daje ogólnego twierdzenia dla wielu kopii.

## Postulat ma sprawdzalną treść

Nowe prawo wyklucza wcześniejszy kontrmodel przerzedzenia przejść odległości
trzy. Przewiduje podwójny iloraz szans równy 1, podczas gdy kontrmodel
przy zadanym g=3 daje 1/100. Jest to rzeczywiste ograniczenie swobody
aktywności **po dodaniu nowego postulatu**. Trzeba równocześnie zachować
niezgodność z pierwotnym HB — nie wolno raportować tylko korzystnego testu.

## Kontrole i granice

Weryfikator przechodzi **48 ograniczonych kontroli dokładnych**: ułamki,
symboliczne logarytmy, źródła 16- i 144-stanowe, przeplatanie z ilorazem,
zmiana przepływu o tych samych marginesach, kontrola pamięci i test ośmiu
stanów. Nie są to 48 dowodów. Nie ma numerycznej symulacji Wszechświata,
dużej enumeracji liczności ani pomiaru fizycznego.

Przeszukano 5286 jawnie zindeksowanych plików pod kątem wcześniejszych
odpowiedników; znaleziono i przeczytano wcześniejsze granice MaxCal oraz
zależności od priory. Znane ograniczenie nie zostało przemianowane na nową
odkrytą zasadę. Nowość lokalna to dopuszczona teraz konkretna hipoteza
i jej opłacony, niekorzystny test FIN. Literatura i awaria Firecrawl mają
osobny zapis w [LITERATURE_AND_REUSE.md](LITERATURE_AND_REUSE.md).

**Strict F01 pozostaje otwarte.** Istnieje już jawny warunkowy model wyboru
historii, lecz nie jest on niezależnym źródłem założonego HB ani prawem
natury. Oryginalny FIN, AGENTS, stare pakiety i plany zachowano bez zmian;
nie wykonano commita ani eksperymentu. Dalsze badanie opisuje
[NEXT_DECISION_AND_AUDIT_ORDER.md](NEXT_DECISION_AND_AUDIT_ORDER.md).
