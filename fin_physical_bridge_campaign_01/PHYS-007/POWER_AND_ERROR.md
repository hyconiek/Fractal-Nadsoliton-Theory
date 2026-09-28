# PHYS-007 — power and error budget

Pierwszy test zamrożono jako **N=2 conditional pair-difference histogram**, `g=3`, plus `g=0` negative control. To test statyczny: anchor `i=0`, próbka `j` ma dokładną predykcję `softmax[(g/2)A[0,j]]`. Nie używa bariery 327, protokołu 333 ani task335.

Primary finite alternatives: full Potts, flat P7, PERT_34_10P, PERT_35_10P. Minimalny idealny TV w tej klasie wynosi 0.011787. Przed walidacją każdy model+aparatura musi wejść w zamrożoną całkowitą kopertę kalibracyjną TV <= 0.003; po odjęciu dwóch kopert pozostaje dodatnia separacja >= 0.005787. Dwie bliższe sondy (5% perturbacja k3/k6 i 2% leakage k1/k2) są tylko sensitivity checks, bo ich kontrast nie przechodzi tego budżetu.

Przy idealnie znanych rozkładach najgorszy primary Chernoff equal-prior bound <=5% wymaga około 13,383 niezależnych prób. To nie jest uniwersalna kontrola type-I. Skrajnie konserwatywny 12-bin union/Hoeffding bound dla empirycznego promienia TV odpowiadającego połowie robust gap wynosi ~13,272,367 prób; pokazuje koszt pełnej distribution-free ostrożności. Realny plan ma użyć prerejestrowanego multinomial likelihood testu/parametric calibration po uzyskaniu rzeczywistej macierzy błędu odczytu.

Finite-g approximation error dla N=2 pair law = 0 wewnątrz modelu. Niepewność `g` ma być kalibrowana niezależnie; cel ±0.2% przesuwa FIN prediction o maks. TV≈0.000923 i musi mieścić się w całej kopercie 0.003.

Jeśli próbki pochodzą z jednej trajektorii, należy użyć ESS z autokorelacji; liczby powyżej zakładają niezależne reset/anchor cycles. Koszt pieniężny pozostaje OPEN, bo brak istniejącej zatwierdzonej platformy i throughput.
