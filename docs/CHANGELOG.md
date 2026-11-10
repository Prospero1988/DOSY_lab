# Changelog

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
