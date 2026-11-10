"""In-app help and native window flags for ILT dialogs."""
from PySide6 import QtCore, QtWidgets

from i18n import current_language, tr


ILT_WINDOW_FLAGS = (QtCore.Qt.WindowType.Window | QtCore.Qt.WindowType.WindowTitleHint
                    | QtCore.Qt.WindowType.WindowSystemMenuHint | QtCore.Qt.WindowType.WindowMaximizeButtonHint
                    | QtCore.Qt.WindowType.WindowCloseButtonHint)

ILT_HELP = {
    'pl': '''
# ILT pojedynczej całki

1. W głównym oknie przetwórz widmo, dodaj całkę i wykonaj dopasowanie. Kliknij **Analizuj ILT** w wierszu sygnału. Potrzeba co najmniej sześciu różnych aktywnych wartości b.
2. Sprawdź zanik i punkty wyłączone z dopasowania. ILT używa bieżącej migawki całek ze znakiem, listy b, maski i wag z dopasowania. Otwórz okno ponownie po zmianie tych danych.
3. Zacznij od ustawień domyślnych. Zmieniaj po jednym parametrze i sprawdzaj jednocześnie rozkład D oraz reszty. Po wpisaniu wartości naciśnij Enter lub przejdź do innego pola.

## Parametry

**Auto λ / stabilność** skanuje zakres λ bez zmiany bieżącego rozwiązania. Kliknij próbę, aby obejrzeć podgląd. Zastosowanie sugestii wymaga potwierdzenia w oknie skanu.

- **λ — siła regularyzacji.** Suwak jest logarytmiczny, a pole przyjmuje dokładną wartość od `1e-12` do `100`. Domyślnie `1e-5`. Mniejsze λ pozwala ściślej dopasować dane, w tym szum, i może ujawnić wąskie, niestabilne maksima. Większe λ mocniej tłumi takie zmiany i może poszerzyć lub zlać maksima.
- **Kara.** Rząd zerowy ogranicza kwadrat gęstości. Druga pochodna (domyślna) ogranicza krzywiznę, więc preferuje płynniejsze rozkłady. Ta sama wartość λ może dawać inny efekt po zmianie kary, zakresu lub typu siatki.
- **D min / D max.** Wpisz dodatnie granice w jednostkach `10⁻¹⁰ m²/s`, z D min poniżej D max. Domyślny zakres obejmuje od D mono/20 do 20 × D mono, z bezpiecznym zakresem zastępczym, gdy D mono wynosi zero. Jeśli maksimum wypada przy brzegu, sprawdź szerszy zakres i wrażliwość na λ.
- **Punkty siatki / Siatka D.** Dostępne są 50, 100, 200 i 400 punktów oraz siatka liniowa lub logarytmiczna. Domyślnie 100 punktów logarytmicznych. Gęstsza siatka dokładniej próbuje rozkład. Rozdzielczość eksperymentalna zależy od danych.
- **Sprawdź wrażliwość ILT** porównuje rozkład dla mniejszej i większej λ oraz drugiej współrzędnej siatki. Raportuje maksima i przedział środkowych 90% masy. To rozrzut modelu, a nie przedział ufności.
- **Logarytmiczna oś D** zmienia tylko sposób rysowania wykresu. **Przywróć ustawienia** wraca do wartości początkowych. **Przelicz ILT** uruchamia obliczenie od razu. Zmiany parametrów przeliczają się też automatycznie po krótkiej przerwie.

## Wykresy i interpretacja

- **A — zanik.** Punkty pomiarowe, zapisana monoeksponenta i rekonstrukcja ILT mają wspólną skalę Iₛ, równą największej bezwzględnej całce wśród aktywnych punktów. Krzyżyki wskazują punkty wyłączone.
- **B — rozkład.** Nieujemne wagi opisują wkład przedziałów D do sygnału. Na wykresie są pokazane jako gęstość o polu jeden w znormalizowanej współrzędnej siatki, liniowej w D lub w log(D). Linia przerywana oznacza D mono w zakresie siatki. W podsumowaniu podano położenia i względne wysokości maksimów.
- **C — reszty.** Porównaj pomiar minus przewidywanie i RMS aktywnych punktów. Sprawdź systematyczne odchylenia oraz komunikaty o zbieżności.

Przeciągnij prostokąt, aby przybliżyć wykres. **Przesuwaj** włącza przesuwanie. Dwuklik lub **Cały zakres** resetuje widok. **Pełny ekran** powiększa wybrany wykres. Użyj **Zapisz PNG** albo **Zapisz wszystkie wykresy PNG**, aby wyeksportować widok.

ILT dopasowuje sumę zaników exp(−bD) z nieujemnymi udziałami i karą Tichonowa. Więcej maksimów ani mniejszy RMS nie wystarczają do rozpoznania dodatkowego związku. Szerokość piku ILT zależy od regularyzacji. Wyświetlane SE dotyczy wyłącznie dopasowania monoeksponencjalnego.
''',
    'en': '''
# Single-integral ILT

1. Process the spectrum, define an integral and fit it in the main window. Click **Analyze ILT** in its row. At least six distinct included b values are required.
2. Check the attenuation and the excluded points. ILT uses a snapshot of the signed integrals, b values, mask and fit weights. Reopen the window after changing those data.
3. Start with the defaults. Vary one setting at a time and check both the D distribution and the residuals. Press Enter or leave the field to apply a typed value.

## Parameters

**Auto λ / stability** scans a λ range without changing the active result. Click a trial to preview it and apply a suggestion explicitly.

- **λ — regularization strength.** Use the log slider or type an exact value from `1e-12` to `100`. The default is `1e-5`. A smaller λ can fit noise more closely and produce narrow, unstable maxima. A larger λ can broaden or merge maxima.
- **Penalty.** Zeroth order penalizes squared density. The default second derivative penalizes curvature. The effect of the same λ can change with the penalty, range and spacing.
- **D min / D max.** Enter positive bounds in `10⁻¹⁰ m²/s`, with D min below D max. The default spans D mono/20 to 20 × D mono, with a fallback when D mono is zero. Check a wider range and sensitivity to λ if a maximum reaches a boundary.
- **Grid points / D grid.** Choose 50, 100, 200 or 400 points and linear or logarithmic spacing. The default is 100 log-spaced points. A denser grid samples the solution more finely. Experimental resolution depends on the data.
- **Check ILT sensitivity** compares smaller and larger λ and the alternative grid coordinate. It reports maxima and the interval containing the central 90% of the mass. That interval describes the fitted distribution.
- **Logarithmic D axis** changes the plot display only. **Restore defaults** resets the settings. **Recalculate ILT** starts a solve immediately. Parameter changes also trigger a solve after a short delay.

## Plots and interpretation

- **A — attenuation.** Measurements, the cached monoexponential fit and the ILT reconstruction share the scale Iₛ, the largest absolute included integral. Crosses mark excluded points.
- **B — distribution.** Nonnegative cell weights describe contributions at different D values. The plotted curve has unit area in the normalized grid coordinate, linear in D or log(D). The dashed line marks the monoexponential D within the range. The summary lists maxima and their relative heights.
- **C — residuals.** Compare measurement minus prediction and the RMS over included points. Inspect systematic deviations and the solver status.

Drag a rectangle to zoom. Use **Pan** to move the view and double-click or **Full range** to reset. **Full screen** enlarges one plot. Use **Save PNG** or **Save all plots as PNG** to export the view.

ILT fits a sum of exp(−bD) decays with nonnegative contributions and a Tikhonov penalty. Extra maxima or a lower RMS alone do not establish extra species. Peak width depends on regularization. The reported SE belongs to the monoexponential fit.
''',
}


GLOBAL_ILT_HELP = {
    'pl': '''
# Global ILT dla wielu całek

**Auto λ / stabilność** porównuje rozwiązania przy stałym μ. Zmiana μ wymaga ponownego skanu.

1. Przetwórz widmo, wyznacz całki i poczekaj na ich dopasowanie w głównym oknie.
2. Zaznacz co najmniej dwa wiersze tabeli przez Ctrl lub Shift. Każdy ROI musi mieć poprawne dopasowanie i co najmniej sześć różnych aktywnych wartości b. Kliknij **Global ILT** obok **Dopasuj / odśwież**.
3. Zacznij od ustawień domyślnych. Porównuj rozkłady, reszty i udziały. Zmieniaj λ i μ osobno, żeby sprawdzić, które maksima utrzymują się po zmianie ustawień.

Analiza używa migawki bieżących całek ze znakiem, zweryfikowanej listy b oraz masek poszczególnych ROI. Po zmianie danych w głównym oknie otwórz analizę ponownie.

## λ i μ

- **λ — regularyzacja Tichonowa** (domyślnie `1e-6`). Karze dużą sumę kwadratów gęstości w każdym ROI. Przy mniejszym λ rozkłady mogą być węższe i bardziej podatne na szum. Przy większym λ udziały zwykle rozkładają się szerzej, a bliskie maksima mogą się zlać.
- **μ — wspólna rzadkość** (domyślnie `1e-4`). Karze normę wektora udziałów wszystkich ROI dla każdego położenia D. Większe μ sprzyja skupieniu udziałów na mniejszej liczbie wspólnych położeń. Pojedynczy ROI może mieć zerowy udział w danej składowej. Zbyt duże μ może usunąć słabą składową lub przesunąć maksimum. Wpisz `0`, aby wyłączyć sprzężenie między ROI.

Suwaki są logarytmiczne. Dokładną wartość możesz wpisać w pole i zatwierdzić Enterem lub zmianą pola. λ jest w zakresie `1e-12`–`100`, a μ od `0` do `100`. Kara μ uwzględnia liczbę ROI, lecz jej praktyczny efekt zależy też od doboru sygnałów, szumu, siatki i zakresu D.

## Siatka, wagi i solver

- **D min / D max** podaj w `10⁻¹⁰ m²/s`. Obie wartości muszą być dodatnie. Zakres początkowy sięga od 1/20 najmniejszego do 20-krotności największego dodatniego D mono. Przy maksimum na brzegu sprawdź szerszy zakres oraz wrażliwość na λ i μ.
- **Punkty siatki / Siatka D.** Wybierz 50, 100, 200 lub 400 punktów i siatkę liniową lub logarytmiczną. Domyślnie 100 punktów logarytmicznych. **Logarytmiczna oś D** zmienia wyłącznie sposób rysowania.
- **Wagi z SE całek.** Jeśli wszystkie aktywne punkty danego ROI mają poprawne SE, używane są odwrotności tych SE. W pozostałych ROI punkty mają równe wagi. Każda kolumna wag jest osobno normalizowana. Silniejsze sygnały są skalowane względem własnej największej bezwzględnej całki, a przy rekonstrukcji amplituda wraca do jednostek pomiaru.
- **Limit iteracji / Tolerancja.** Domyślnie 10 000 i `2e-6`. Tolerancja dotyczy znormalizowanej reszty stacjonarności. Zmniejszenie jej lub zwiększenie limitu może wydłużyć obliczenia. Po osiągnięciu limitu wynik jest oznaczony jako przybliżony. Panel podaje wartość funkcji celu, resztę stacjonarności i liczbę iteracji.
- **Porównaj niezależne ILT** dopasowuje każdy ROI osobno na tej samej siatce, z tym samym λ, maską i wagami. **Przywróć ustawienia** przywraca wartości początkowe. **Przelicz ILT** uruchamia obliczenie od razu. Zmiany parametrów uruchamiają je też po krótkiej przerwie.

## Jak czytać wyniki

- **A — rozkład wspólny.** Wyświetla gęstość z normy udziałów ROI przy każdym D. Zaznaczone maksima to kandydaci do interpretacji.
- **B — rozkłady ROI.** Menu **Pokaż rozkłady ROI** przełącza widoczność poszczególnych krzywych. Domyślnie każda ma pole jeden. **Wspólna skala amplitudy ROI** zachowuje proporcje odtworzonych amplitud. Kropkowane linie wskazują położenia wspólnych maksimów.
- **C — rekonstrukcja.** Wybierz ROI z listy. Wykresy pokazują punkty pomiarowe, zapisaną monoeksponentę, Global ILT, opcjonalne niezależne ILT oraz ich reszty. Wszystkie mają wspólną skalę dla danego ROI. RMS liczy się z aktywnych punktów.
- **D — udziały.** Każda komórka sumuje wagi w przedziale D wyznaczonym między sąsiednimi minimami rozkładu wspólnego. Kolumny domyślnie pokazują procent całego odtworzonego sygnału ROI. Można przełączyć na oryginalne jednostki całek dla przewidywanego b = 0. Najedź kursorem, by odczytać wartość.
- **E — wyniki.** Tabela zestawia ppm, D mono, dominujące D, udziały w kolejności rosnącego D, RMS i status obliczeń.

Na wykresach możesz zaznaczyć prostokąt, przesuwać, zresetować zakres dwuklikiem albo użyć pełnego ekranu. **Zapisz PNG** zapisuje pojedynczy wykres, a **Zapisz wszystkie wykresy PNG** tworzy wspólną figurę z aktualnie wybranym ROI w panelu C.

Wspólne D i zbieżność solvera nie wystarczają do przypisania cząsteczki. Udziały opisują całki NMR, a ich interpretacja zależy od jakości danych, zakresu b, nakładania sygnałów i ewentualnej wymiany. To analiza wybranych ROI. Nie odtwarza widma 2D w rozdzielczości ppm.
''',
    'en': '''
# Global ILT for selected integrals

**Auto λ / stability** compares solutions at fixed μ. Run a new scan after changing μ.

1. Process the spectrum, define the integrals and complete their fits.
2. Ctrl-click or Shift-click at least two rows. Each ROI needs a valid fit and six distinct included b values. Click **Global ILT** next to **Fit / refresh**.
3. Start with the defaults. Inspect distributions, residuals and contributions. Vary λ and μ separately to check which maxima persist.

The window uses a snapshot of signed integrals, verified b values and ROI-specific masks. Reopen it after editing the main data.

## λ and μ

- **λ — Tikhonov regularization** (default `1e-6`). It penalizes the squared density within each ROI. Lower λ can yield sharper, noise-sensitive maxima. Higher λ can broaden distributions and merge nearby maxima.
- **μ — joint sparsity** (default `1e-4`). It penalizes the vector norm across ROIs at each D. Larger μ favors fewer shared positions and allows a given ROI to contribute zero at any position. Excessive μ can suppress weak populations or shift maxima. Type `0` to turn off coupling.

The sliders are logarithmic. Type an exact value and press Enter or leave the field to apply it. λ ranges from `1e-12` to `100` and μ from `0` to `100`. The μ penalty accounts for the number of ROIs. Its practical effect also depends on the selected signals, noise, grid and D range.

## Grid, weighting and solver

- **D min / D max** are positive bounds in `10⁻¹⁰ m²/s`. The initial range spans 1/20 of the smallest to 20 times the largest positive monoexponential D. Check the range and sensitivity to λ and μ if a maximum reaches an edge.
- **Grid points / D grid.** Choose 50, 100, 200 or 400 cells and linear or logarithmic spacing. The default is 100 log cells. **Logarithmic D axis** changes the display only.
- **Weight by integral SE.** ROIs with a complete set of valid included uncertainties use inverse SE. The rest use equal point weights. Each weight column is normalized separately. Intensities are scaled by the largest absolute included integral within each ROI and converted back for reconstruction.
- **Iteration limit / Tolerance.** The defaults are 10,000 and `2e-6`. Tolerance applies to the normalized stationarity residual. Tightening it or raising the iteration limit can take longer. The panel reports objective value, stationarity residual and iterations. Iteration-limited results are labeled approximate.
- **Compare independent ILT** fits each ROI separately using the same grid, λ, masks and weights. **Restore defaults** resets the settings. **Recalculate ILT** solves immediately. Edits also trigger a solve after a short delay.

## Read the views

- **A — joint distribution.** The plotted density comes from the norm across ROI contributions at each D. The markers show candidate maxima.
- **B — ROI distributions.** Use **Show ROI distributions** to control visibility. Each curve has unit area by default. **Common ROI amplitude scale** retains the ratios of reconstructed amplitudes. Dotted lines mark shared maxima.
- **C — reconstruction.** Choose an ROI to see observations, the cached monoexponential fit, Global ILT, optional independent ILT and their residuals. They use a shared scale within the ROI. RMS uses the included points.
- **D — contributions.** Each cell sums fitted weights in an interval divided by minima of the joint distribution. Each ROI column defaults to percentages of its fitted total. Original integral units refer to the fitted prediction at b = 0. Hover over a cell to read it.
- **E — results.** The table lists ppm, monoexponential D, dominant D, contributions by increasing D, RMS and solver status.

Use rectangle zoom, pan, full range (or double-click) and full screen on the plots. **Save PNG** saves one plot. **Save all plots as PNG** combines them with the current ROI in panel C.

Shared D and solver convergence cannot establish molecular identity. Contributions describe NMR integrals and depend on data quality, b coverage, signal overlap and exchange. This analysis covers selected ROIs and does not reconstruct a ppm-resolved 2D spectrum.
''',
}


def add_ilt_info(layout, parent, title, manual):
    button = QtWidgets.QPushButton('Info')
    button.setObjectName('ilt_info_button')
    button.setToolTip(tr('Otwórz instrukcję ILT i opis parametrów.'))
    button.clicked.connect(lambda: show_ilt_info(parent, title, manual))
    layout.addWidget(button)


def show_ilt_info(parent, title, manual):
    dialog = parent.findChild(QtWidgets.QDialog, 'ilt_info_window')
    if dialog is not None and dialog.isVisible():
        dialog.findChild(QtWidgets.QTextBrowser).setMarkdown(manual[current_language()])
        dialog.raise_()
        dialog.activateWindow()
        return
    dialog = QtWidgets.QDialog(parent, ILT_WINDOW_FLAGS)
    dialog.setObjectName('ilt_info_window')
    dialog.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose)
    dialog.setWindowTitle(tr(title))
    dialog.resize(820, 680)
    layout = QtWidgets.QVBoxLayout(dialog)
    browser = QtWidgets.QTextBrowser()
    browser.document().setDocumentMargin(14)
    browser.setMarkdown(manual[current_language()])
    layout.addWidget(browser, 1)
    close = QtWidgets.QPushButton(tr('Zamknij'))
    close.clicked.connect(dialog.close)
    layout.addWidget(close, alignment=QtCore.Qt.AlignmentFlag.AlignRight)
    dialog.show()
