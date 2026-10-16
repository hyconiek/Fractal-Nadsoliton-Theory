# Wyniki nowych badań: rzeczywista lokalna domena FIN

Ustalono nowe twierdzenie matematyczne. W rzeczywistym obszarze metastabilnym MT każdy z dwunastu udziałów populacji jest co najmniej `103/1250` (z zachowaniem ściślejszej nierówności). Dlatego dla **każdego legalnego przejścia pomiędzy licznościami w domenie** natywny heat-bath ma tempo co najmniej `a_0 N`, gdzie `a_0=(103/15000)e^(-24)≈2,59·10^-13`. Liczba jest bardzo konserwatywnym dolnym oszacowaniem, nie przewidywaną obserwowaną intensywnością.

To usuwa lukę między dwoma wcześniejszymi sposobami sformułowania zadania: ograniczeniem względnym przewodności przez `K` razy przewodności HB i pierwotnym ograniczeniem sumy intensywności przez `C N`. **W tej konkretnej lokalnej domenie** ograniczoność jednego warunku dla wszystkich N jest równoważna ograniczoności drugiego. Poprzednie przykłady rozbieżności dotyczą innych grafów bez dodatniej dolnej granicy na referencyjne krawędzie.

Jest także dokładny trzeci, wygodniejszy warunek. Dla każdego zbioru całych orbit symetrii liczymy sumę rzeczywistego źródła dryfu parzystości, a w mianowniku sumę wag Gibbsa mnożoną przez liczbę dozwolonych skierowanych wyjść z tych stanów. Maksimum tych ilorazów to `Q_N`. Udowodniono:

`Q_N <= C_N^vert <= K_N^* <= Q_N/a_0`.

W rezultacie wystarczy **i jest konieczne** wykazać `sup_N Q_N<∞`, by rozstrzygnąć istnienie jednolitego limitu intensywności na dotychczasowym grafie RB. Ocena `Q_N` dla wszystkich ogromnych rzeczywistych grafów MT pozostaje nierozstrzygnięta. Nie skonstruowano nowych generatorów ani nie znaleziono przeciwprzykładu dla IBA.

Pakiet zawiera dokładny dowód, małe wymierne kontrole (niebędące modelem q12), zakres nierozstrzygniętego lematu, rejestr roszczeń i kod weryfikacyjny. Pierwotny model i AGENTS nie są zmieniane; status fizyczny ToE pozostaje bez zmian.
