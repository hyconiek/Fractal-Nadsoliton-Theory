# RB i RBCUT: co przyjęto i co rzeczywiście wyprowadzono dalej

Data: 2026-10-16. Wspólny przegląd obu wskazanych pakietów. Oryginalne
pliki, archiwa oraz wcześniejsze zmiany w repo pozostają nienaruszone.

Najważniejszy postęp nie polega na kolejnym przepisaniu problemu przepływu.
Ostatni etap RBCUT dowodzi, że w **rzeczywistym lokalnym obszarze FIN**
natywne szybkości pojedynczych przejść mają dolne ograniczenie a0 N,
niezależne od N. Dlatego ograniczenie względem natywnego przepływu i
pierwotny warunek całkowitej szybkości CN są równoważne co do istnienia
jednolitej stałej. Starszy kontrprzykład ogólnografowy tego nie obala:
jego szybkości referencyjne degenerują, inaczej niż w lokalnym FIN.

To nadal **nie jest dowód istnienia potrzebnej rodziny kinetyk**.
Nie ograniczono jeszcze wszystkich przekrojów rzeczywistego grafu liczności.
Uniwersalne uśrednianie nawiasu IBA bez dodatkowych warunków pozostaje
nierozstrzygnięte; twierdzenia KR pod NC/IBA zachowują ważność.

## Poprawki zamiast cichego przyjęcia

- W czterowęzłowym przykładzie całkowita szybkość w środkowym węźle wynosi
  2(1+1/N), nie 1+1/N. Sens jakościowy przykładu pozostaje poprawny.
- Opis znaku potencjału dualnego wymaga odwrócenia słownego kierunku;
  jego poprawny wzór nie wymaga zmiany.
- Rzeczywiste dane Gibbsa nie dają automatycznie przepływu wymiernego.
- Przejścia zachowujące parzystość mogą przecinać fazy. W konstrukcji
  o rzadkim przecieku wolno zachować tylko te zgodne z jej podporą.

Wszystkie korekty są zapisane w osobnym audycie, nie w archiwalnych źródłach.

## Własna praca wykraczająca poza sprawdzenie paczek

Wyprowadziłem ograniczenia dla trzech klas przekrojów: globalnych prefiksów
parzystości, zbiorów niezmienniczych względem nieparzystego obrotu etykiet
oraz przekrojów z dostatecznie odsłoniętym dozwolonym brzegiem. Pozostaje
ściślej określona klasa nieregularnych przekrojów izolujących interfejs.
Jej istotą jest **nadmiar przepływu zabronionego na zewnątrz po odjęciu
całego natywnego przepływu do środka**, a nie sama mała masa powłoki.

Osobno dowiodłem, że proste usunięcie niedozwolonych przejść ma niezanikający
błąd dryfu w ważonej skali krytycznej: dolne ograniczenie wynosi M/2>0.
Jego nieskalowany błąd może jednocześnie dążyć do zera. To konkretna
przeszkoda dla pozornego domknięcia przez „małą poprawkę na powłoce”,
nie dowód niemożliwości każdej naprawy.

## Następne zlecenie

[SC-001–SC-004](../FIN_SIGNED_INTERFACE_ORBIT_CUT_NEXT_CHAIN_20261016.md):
udowodnić jedną nierówność kompensacji dla pozostałych rzeczywistych
przekrojów albo skonstruować prawdziwy ciąg przekrojów łamiących ograniczenie.
Dopiero dodatni wynik uruchamia budowę szybkości i wniosek o IBA. Jeśli
nierówność pozostaje otwarta, oddać precyzyjny brakujący lemat i zatrzymać
tę próbę; nie zamieniać kolejnego równoważnego zapisu w pozorny postęp.

Sprawdzono 47 wpisów manifestu RB, 75 wpisów pakietu zbiorczego i 71 plików
w czterech oryginalnych ZIP-ach. Pięć skryptów naukowych odtwarza wyniki
bajtowo; niezależny weryfikator ma 22 ograniczone kontrole. Jeden oryginalny
sprawdzacz paczki zatrzymuje się na historycznej ścieżce /mnt/data; jego
manifest i dostępne wejścia sprawdzono osobno, bez ogłaszania jego PASS.
Kontrole skończone nie zastępują dowodu dla całej rodziny N.

Zysk dla fizyki to precyzyjniejsze rozpoznanie, które prawa kolektywne
wynikają ze statyki i dryfu, a które zależą od swobodnej kinetyki. Nie
otrzymujemy stąd źródła A7, przestrzeni, jednostek, zachowanej energii,
fizycznego zegara, QM/GR ani ToE. Pakiet jest włączony do warstwy przeglądu
i rejestru AGENTS w drzewie roboczym; ten przegląd nie wykonuje commita.
