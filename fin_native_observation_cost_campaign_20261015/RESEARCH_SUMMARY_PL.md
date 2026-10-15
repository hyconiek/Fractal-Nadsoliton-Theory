# NOC-001–NOC-005 — wyjaśnienie po polsku

**Nowy wynik (warunkowy):** dla dwóch wcześniej ustalonych matematycznych wariantów FIN można odróżnić je z prawdopodobieństwem błędu dążącym do zera, obserwując coraz mniejszy procent elementów, ale rejestrując ich zachowanie odpowiednio długo. Nie jest to jednak test oparty wyłącznie na znakach: wykorzystuje średnią z `cos(W)`, więc zakłada dokładną skalę amplitudy.

Z góry ustalono `k_N=floor(N^(15/16))` obserwowanych elementów i `M_N=floor(N^(1/4))` odczytów. Ich udział `alpha_N=k_N/N` maleje. Bez dodawania nowych sił i szumu wyprowadzono błąd klasyfikacji ograniczony trzema malejącymi potęgami `N^-1/8, N^-3/16, N^-13/32`, składnikiem brzegowym oraz błędem wcześniejszej lokalnej relaksacji i sprzężenia. Ważne jest wykorzystanie korektora Poissona i martyngału, a nie fałszywe założenie, że kolejne odczyty są niezależne.

**Ograniczenie:** wynik dziedziczy warunki matematyczne MT/KO/DKO; nie podaje praktycznej liczby elementów, czasu fizycznego ani rozwiązania problemu, skąd przyroda miałaby wziąć operator A7. Kolejne możliwe badanie dotyczyłoby *koniecznego* kosztu obserwacji, a nie dalszej weryfikacji tej samej wystarczającej konstrukcji.
