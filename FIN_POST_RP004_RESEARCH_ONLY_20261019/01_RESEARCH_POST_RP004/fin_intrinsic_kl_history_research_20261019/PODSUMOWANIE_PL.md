# FIN — podsumowanie badań z 19.10.2026 (dla laika)

**Pytanie:** Czy nadsoliton jako ewolucja skorelowanych procesów musi się samodzielnie zmieniać według prawa wynikającego z informacji?

**Wynik:** Przebadano jeden konkretny wcześniejszy pomysł FIN (GEN-654): minimalizowanie sumarycznej różnicy informacyjnej między kolejnymi stanami. Przy matematycznie poprawnej normalizacji daje to energię kinetyczną informacji (metryka Fishera–Rao). Jeżeli znany jest tylko stan początkowy, najmniejszy koszt ma brak zmian. Jeśli znamy również stan końcowy, otrzymujemy wyznaczoną geometrycznie drogę między rozkładami, ale warunek końcowy trzeba dostarczyć dodatkowo. Taka droga nie wybiera kierunku czasu.

Sprawdzono też dwa możliwe sposoby uratowania pomysłu:

1. **Czy naturalna droga najmniejszego kosztu odtwarza dyfuzję ścisłego operatora FIN?** Nie: dla żadnego niejednostajnego stanu początkowego trajektoria cieplna nie minimalizuje samego działania Fishera pomiędzy własnymi końcami.
2. **Czy gradient entropii Fishera odtwarza znaną dyfuzję?** Nie: przy równowadze daje jednakowe tempo zaniku wszystkich perturbacji, podczas gdy FIN ma sześć różnych dodatnich wartości spektralnych. Sama zmiana jednostki czasu nie wystarcza.

**Pozytywny element:** Jeżeli wprowadzimy **określoną, dostosowaną do grafu mobilność z logarytmiczną średnią**, dyfuzja FIN staje się dokładnym przepływem gradientowym entropii. Ale sama ta mobilność nie jest jeszcze wyprowadzona z ontologii; dobiera się ją tak, aby uzyskać znane równanie. To matematycznie pokazuje, jakiego typu konkretnego obiektu brakuje.

**Test liczbowy:** dla jednego jawnego początkowego rozkładu: działanie rzeczywistej trajektorii cieplnej ≈ 0,243273, natomiast koszt najkrótszej drogi łączącej te same końce ≈ 0,186642. Różnica jest dodatnia. Dla wszystkich dopuszczalnych rozkładów nierówność ma niezależny dowód analityczny.

**Znaczenie dla ToE:** To warunkowy wynik matematyczny, a nie potwierdzenie fundamentalności FIN. Ściśle odrzuca niektóre *sposoby wyprowadzania* dynamiki z informacji, nie intuicję o informacyjnym nadsolitonie. Priorytet F01 (niekołowe źródło ewolucji) nadal jest otwarty.

**Metodologia:** Źródła AGENTS i GEN rzeczywiście odczytano; wcześniej znane przeszkody z teorii aktywności i wyboru źródła zachowano jako poprzedniki, nie nowe odkrycia. Matematyka została zapisana po angielsku w `PAPER_INTRINSIC_DYNAMICS_FIN.md`; `verify.py` to kontrola ograniczona, nie formalizacja dowodu.
