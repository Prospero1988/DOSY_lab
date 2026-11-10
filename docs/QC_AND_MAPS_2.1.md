# DOSY Lab 2.1 · Kalibracja QC i nieprzypisane sygnały na mapach ILT

## Strychnina w CDCl₃

Odtworzyłem 20 ROI z projektu użytkownika. Porównanie przed i po zmianie wykorzystuje ten sam JDF, granice ROI, 16 warstw, maski fitów, fazę i korekcję bazy zapisaną w projekcie. Referencją w pliku jest warstwa 2. Dla porównania widoczności przyjąłem podłogę 0% i wzmocnienie 1×. Maksymalna różnica między 20 wartościami D przed i po zmianie wynosi dokładnie zero.

| Wynik | Przed | Po |
| --- | ---: | ---: |
| QC `UWAGA` | 20 | 2 |
| QC `INFO` | 0 | 18 |
| Mapa oryginalna, wyrównana i Shared-D, każda osobno | 20 | 20 |
| Mapa Global ILT | 4 | 18 |
| Mapa Adaptive Local ILT | 6 | 16 |

Przed zmianą 2 z 6 ROI na mapie AILT, S2 i S11, pochodziły z niezależnego dopasowania monoeksponencjalnego. Pozostałe 4 pochodziły z rekonstrukcji ILT. Po zmianie wszystkie 16 ROI na mapie AILT pochodzi z rekonstrukcji. Dla Global ILT i AILT liczba widocznych sygnałów była równa liczbie ROI z dodatnimi danymi na mapie przy ustawieniach powyżej. W tym projekcie maskowanie wspólnych binów ani podłoga wyświetlania nie odpowiadały za utratę sygnałów.

[Wyniki przed zmianą](../benchmarks/qc_maps_before.json) i [wyniki po zmianie](../benchmarks/qc_maps_after.json) podają dla każdego ROI granice ppm, D, R², RMS, poziom QC, wartości i progi diagnostyk kształtu, fazy i bazy, S/N, poziom szumu fazy i porównania D. Zawierają też przejście od obecności widma i całki przez fit mono, udział w optymalizacji, akceptację profilu, macierz mapy aż do wyświetlenia sygnału.

## Co wywoływało ostrzeżenia

Przed zmianą wystarczało wykrycie trendu z b, nawet jeśli różnica była zbyt mała, żeby zmienić D w istotnym stopniu. Dawne testy wykryły trend bazy w 17 ROI, zmianę kształtu w 16 i zmianę fazy w 3. Po zmianie liczby wykrytych trendów wynoszą odpowiednio 17, 18 i 3. Testy pozostają aktywne. S/N w projekcie był wystarczający do diagnostyki wszystkich 16 warstw każdego ROI. Główną przyczyną nadmiaru ostrzeżeń były rzeczywiste, lecz małe odchylenia i zbyt surowe przypisanie im statusu, a nie szum na końcu zaniku.

- **Kształt i faza.** Zespolony profil porównuję z ważonym wzorcem rzędu 1. Szum obu kwadratur jest szacowany osobno dla warstw na podstawie cichszych bloków widma po usunięciu lokalnego trendu. Do porównania włączam warstwy z konserwatywnie oszacowanym S/N ≥ 5 i wymagam co najmniej ośmiu różnych b. Warstw z niewystarczającą informacją nie usuwam z dopasowania zaniku. W diagnostyce fazy odejmuję wspólne przesunięcie i podaję względną fazę oraz jej skalę szumową. Losowe zmiany znaków całych warstw zachowują korelację wzdłuż ppm. Próg dla wykrycia trendu to 0,05 podzielone przez liczbę testowanych efektów. To reguła przesiewowa oparta na 511 realizacjach, przybliżonej symetrii i niezależności błędów między warstwami.
- **Wpływ na D.** Dla kształtu porównuję zaniki dwóch części wzorca o zbliżonej energii. Dla fazy porównuję całki ze znakiem z modułem współczynnika projekcji. Dla bazy porównuję D z korektą i bez niej na tych samych granicach, b, masce i wagach. Raportuję też stosunek RMS korekty do RMS całek przed korektą. Diagnostyka wpływu stosuje skalę przesiewową 1,96 × (SE₁ + SE₂), bez założenia niezależności obu fitów. Nie jest to przedział ufności dla różnicy. Jeśli efekt nie przekracza tej skali albo brak danych do oceny, status to INFO. Istotna w tym porównaniu różnica daje UWAGA. Niedostępny lub nieważny model ma osobny status, a główny fit nadal korzysta z oryginalnych całek ze znakiem.
- **Widok i eksport.** W tabeli najedź na QC, aby zobaczyć osobno diagnozę dopasowania i pomiaru. Uwagi diagnostyczne trafiają też do eksportu wyników.

## Gdzie znikały sygnały

W danych wejściowych, całkach i fitach mono zachowano wszystkie ROI. Większość zniknięć występowała po optymalizacji, kiedy sam alarm kształtu lub fazy wykluczał przypisanie profilu do pojedynczego D. Teraz taki efekt blokuje pojedyncze przypisanie, gdy wpływa na D. Resztę widmową oceniam po odjęciu oczekiwanej energii szumu od sumy kwadratów. Dla profilu przypiętego do warstwy referencyjnej uwzględniam również propagację szumu tej warstwy. Zachowałem próg względnej reszty 15% dla nadwyżki ponad szum. Przy bardzo słabej referencji, nieznanym szumie lub silnie skorelowanych błędach modelu to kryterium ma ograniczenia. S/N referencji poniżej 5 blokuje rekonstrukcję, kiedy można go oszacować.

Sam dobry zanik całego przedziału ILT nie wystarcza, by narysować sygnał przy D jego maksimum. Dla pojedynczej narysowanej składowej sprawdzam oddzielnie zanik dla tej wartości D. Zachowałem wcześniejszą tolerancję RMS względem indywidualnego dopasowania, większą z 1,5 × RMS i RMS + 0,005, w znormalizowanych jednostkach. Przy wielu składowych udziały poniżej 5% mogą uczestniczyć w dopasowaniu profili jako składniki pomocnicze, lecz ich masy nie przenoszę na inne D. Gdy po pominięciu dodatniego udziału zostaje jedno maksimum, ROI pozostaje bez przypisania.

| ROI | Ustalenie |
| --- | --- |
| S2, 7,186–7,294 ppm | D mono = 24,270 × 10⁻¹⁰ m²/s przy R² = 0,999977. Odczyty kształtu i fazy zmieniają D o około 1,05 i 1,06 × 10⁻¹⁰ m²/s, wobec progów 0,259 i 0,144 × 10⁻¹⁰ m²/s. Przy wspólnym maksimum Global ILT, 21,278 × 10⁻¹⁰ m²/s, względny RMS wynosi 2,75% przy limicie 0,62%, mimo że dla całego przedziału wynosi 0,35%. Zgodnie z przypisaniem użytkownika w S2 nakładają się sygnały strychniny i chloroformu. Podane D opisuje zanik całej całki S2. Ani dopasowanie monoeksponencjalne, ani ta rekonstrukcja ILT nie dostarcza osobnych D obu składników. |
| S11, okolice 0 ppm | D mono = 30,216 × 10⁻¹⁰ m²/s. Wrażliwość D na kształt i odejmowanie bazy przekracza odpowiednie progi. Wspólne maksimum nie odtwarza zaniku. AILT nie znalazł stabilnego niezależnego rozdzielenia. |
| S1 i S17 | Global ILT przypisuje profile. W AILT zanik przy maksimum 21,482 × 10⁻¹⁰ m²/s przekracza limit RMS odpowiednio 0,991% vs 0,936% i 0,862% vs 0,836%. Pozostają nieprzypisane w AILT przy tych ustawieniach. |

Domyślna nakładka pokazuje nieprzypisane ROI w bursztynie na zmierzonym widmie 1D. Opcjonalne fioletowe kreski na mapie wskazują pozycję mono D. Kreski nie mają intensywności rekonstrukcji i nie zasilają projekcji D. Z listy można odczytać powód, liczby i użyte progi. Legenda oraz PNG rozróżniają rekonstrukcję, nakładkę oraz pozycje mono D.

## Kontrole syntetyczne

80 realizacji na warunek, ziarno 2102, zmienny szum między warstwami i korelacja szumu wzdłuż ppm. `Alarm przed` oznacza dawny surowy alarm. `UWAGA po` oznacza nową klasę ryzyka dla D. Odrzucenie dotyczy przypisania pojedynczego profilu w tej kontroli, nie utraty sygnału z widma.

| Przypadek | Alarm przed / 80 | UWAGA po / 80 | Odrzucenie mapy przed → po |
| --- | ---: | ---: | ---: |
| Pojedynczy profil, mały szum | 5 | 0 | 5 → 0 |
| Pojedynczy profil, szum na końcu zaniku | 2 | 0 | 80 → 0 |
| Stała faza 25° | 5 | 0 | 5 → 0 |
| Faza zależna od b, zakres 55° | 80 | 74 | 80 → 74 |
| Nakładanie dwóch profili o różnym D | 80 | 79 | 80 → 80 |

Przy zaszumionym końcu zaniku wszystkie 80 realizacji dało INFO z powodu niedostatecznej informacji do diagnostyki. Przy zmiennej fazie trend wykryto w 80/80. Sześć realizacji pozostaje na poziomie INFO i przechodzi próbę mapy, mimo wprowadzonego zaburzenia. Przy nakładaniu jeden przypadek kończy się jako INFO, lecz próba mapy odrzuca wszystkie pojedyncze przypisania. Wyniki określają zachowanie dla tych danych syntetycznych. Potrzeba niezależnych pomiarów, by wyznaczyć częstość pomyłek dla innych próbek, sprzętu i zakresów b. Kształt odzyskany przez ILT zależy od zakresu b, siatki i regularyzacji. Zgodność zaniku i profilu nie dowodzi tożsamości składnika chemicznego. [Dane symulacji](../benchmarks/qc_maps_synthetic.json)

## Odtworzenie

Z katalogu repozytorium, w PowerShellu z zainstalowanymi zależnościami aplikacji.

```powershell
$env:OPENBLAS_NUM_THREADS = '1'
& '.venv/Scripts/python.exe' benchmarks/qc_maps.py 'projects/Strychnina_CDCl3_bpp_ste_led_diffusion-1-1.proj'
& '.venv/Scripts/python.exe' benchmarks/qc_synthetic.py
```

Wynik `before` odtworzono z kodu `aedf2b8` za pomocą `qc_maps.py --before <projekt.proj>`. Projekt ma SHA-256 `1031695f62961dd3aedf5194004e6536019c0ab2e0a0cb4cdc0201114cb18121`, a źródłowy JDF `d65b9384caddd4439fc916a4331bbd7226ac83eeb4d5de1d4242c4566d525b55`. Benchmark bez argumentu używa automatycznie wykrytych ROI wzorcowego JDF, więc ten tryb nie odtwarza zapisanych ustawień użytkownika. Testy `test_scientific_remediation.py` odtwarzają kontrprzykłady z audytu `42895f1`, a `test_qc_maps.py` obejmuje rozróżnienie szumu od zmiany profilu, błędne maksimum wspólnego przedziału i pominiętą masę. Testy GUI sprawdzają oddzielną prezentację i eksport nieprzypisanych sygnałów.

Na gałęzi 2.1 przeszło 98 testów, w tym 61 testów obliczeń i 37 testów GUI. Ponownie uruchomiłem przypadki z audytu dotyczące metadanych b, processingu, map AILT i odtwarzalności eksportu. Składnik o udziale 2,91% pozostał nierozdzielony, a oba profile z rozstrzygalnego przypadku kontrolnego zostały odtworzone z maksymalnym błędem 4,45 × 10⁻¹⁶ w jednostkach tego testu. Wspólne biny nakładających się ROI nie dostały intensywności w rekonstrukcji. Całki i dane w eksporcie odtworzyły się dokładnie w trzech przykładowych JDF. W osobnym teście GUI otworzyłem zapisany projekt strychniny. Tabela pokazała 18 INFO i 2 UWAGA. W AILT 16 ROI znalazło się w rekonstrukcji, a cztery pozostałe w nakładce.
