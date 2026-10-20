# FIN — cztery zdarzenia wystarczą, aby rozpoznać różnicę niewidoczną w powrotach

Kontynuacja poprzedniej kampanii z dwoma różnymi generatorami na 12 stanach o identycznych prawdopodobieństwach powrotu dla dowolnego czasu i liczby skoków.

**Nowy wynik:** obserwator, który nie zna nazw stanów i tylko sprawdza, czy kolejne stany są równe, odróżni te dwa modele dopiero przy czterech kolejnych stanach. Wzorzec `A→B→A→B` ma prawdopodobieństwa różniące się dokładnie o `11/34629528` (~3,18·10^-7). Całe rozkłady pięciu możliwych wzorców równości dla czwórki stanów różnią się w odległości całkowitej wariacji `11/17314764` (~6,35·10^-7). Wzorce najwyżej trzech stanów są identyczne w obu modelach. Nie potrzeba sekund; numerujemy rzeczywiste zmiany, ale sam proces Markowa jest wcześniej założony.

**Drugi wynik:** jeśli znamy częstości naprzemiennych historii ABAB... do długości 12 stanów, z dokładnych tożsamości Newtona odzyskujemy wszystkie 11 wartości prawdopodobieństw wyjścia z danego stanu (wraz z krotnościami). Są one jednak nieuporządkowane. Nie wiemy, które przypadają na konkretne odległości w pierścieniu. Dla dwóch jawnych modeli można mieć te same częstości naprzemiennych historii dla wszystkich długości, a inną częstość powrotu po trzech skokach.

**Znaczenie:** historia skorelowanych zdarzeń może zawierać informacje o organizacji procesu niewidoczne w prostych korelacjach; nie jest to jednak wyprowadzenie prawa dynamiki z pierwotnego nadsolitonu. F01/F03/F06 pozostają otwarte. Wszystkie wyniki dotyczą zadanej klasy 12-stanowych modeli Markowa, nie natywnego modelu N2, nie przestrzeni fizycznej ani obserwacji eksperymentalnych.
