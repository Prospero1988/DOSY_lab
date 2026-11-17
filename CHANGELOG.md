# Changelog

## 2.1.4 · 2026-11-17

- Highlight blocking errors in red and warnings in amber. Show both ppm intervals when a baseline mask intersects an integral.
- Show a baseline extrapolation warning only when an integral extends beyond the fitted range and baseline subtraction is applied. Name the affected integrals in the warning.
- Use **INFO** for measurement diagnostics and reserve **WARNING** for issues with the decay fit. Mark unavailable current results with a dash in the main QC table.
- Mark draggable panel dividers with 1 px lines that brighten on hover.

## 2.1.3 · 2026-11-17

- Add a Processing preference to show or hide the blue baseline regions independently of the fitted baseline preview.
- Put **Select baseline** and **Clear baseline** in Processing and remove their main-toolbar buttons.

## 2.1.2 · 2026-11-17

- Add optional per-layer solvent-tail baseline correction for presaturated proton DOSY with a recorded irradiation position in ppm.
- Guide background selection in an amber Processing panel. Show or hide the orange wing guides independently of the blue mask and the yellow dashed fitted baseline.
- Apply manually entered ILT regularization values with **Enter**, including λ in individual ILT and λ and μ in Global ILT.
- Show a confirmation message and close the regularization scan when no clear λ recommendation is available.

## 2.1.1 · 2026-11-13

- Add main-window shortcuts **M** to edit, **D** to split and **R** to delete the selected integral.
- Add a dedicated **Set pivot** mode to the phase correction dialog.
- Confirm exit and offer to save the current project before closing.

## 2.1.0 · 2026-11-12

- Withhold conflicting automatic JEOL b values and support explicitly sourced manual b lists.
- Stabilize bounded monoexponential fits for shifted and late b ranges, with conditional uncertainty diagnostics.
- Separate fit and measurement QC. Add noise-aware complex-profile diagnostics and baseline sensitivity checks.
- Add Exploratory, QC-filtered and Mono-fit reference policies for ILT maps, with unresolved-ROI overlays and policy-labelled PNG exports.
- Add opt-in λ scans, noise-based discrepancy checks, conditional fixed-μ scans, peak stability diagnostics and per-group AILT recommendations.
- Organize the public user guide, scientific documentation, tests, reproducible benchmarks and the historical 2.0.0 audit.

## 2.0.0 · 2026-11-10

- Added exploratory Individual ILT with adjustable D grid and penalties, decay and residual comparison, solver diagnostics, and PNG export.
- Added exploratory Global ILT for multiple selected ROIs with joint sparsity, ROI distributions, a component contribution matrix and independent-fit comparison.
- Expanded the 2D map with display alignment, shared monoexponential fitting, Global ILT profiles and Adaptive Local ILT. The original monoexponential map remains the default.
- Added a linked optional diffusion projection, guide clearing with Delete, and Z/P/F navigation for the map and phasing window.
- Added reusable integral definitions with `.ints` import and export. Templates contain ppm bounds, IDs and groups. Fits are recomputed on the current measurement.
- Added in-app ILT and Global ILT guides, and enabled native title-bar close buttons on both analysis windows.
- Preserved the conventional fitting method, raw JEOL JDF support and the existing project and processing file schemas.

## 1.0.0

- Classical JEOL DOSY processing, conventional diffusion fits, D statistics, a native-D map and analysis exports.
- English and Polish interface.
