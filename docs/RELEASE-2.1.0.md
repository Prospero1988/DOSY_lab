# DOSY Lab 2.1.0 release notes

DOSY Lab analyzes classical JEOL diffusion NMR from raw JDF acquisitions. Version 2.1.0 strengthens validation of the diffusion-weighting metadata and gives users more control over how exploratory ILT results are shown and assessed.

- **Diffusion weighting and monoexponential fits.** The reader checks independent sources of b when available. Conflicting inputs leave b unresolved until a sourced manual list is supplied. The fit uses the lowest included b as its numerical amplitude reference and retains the original signed integrals and point selections.
- **Scientific QC.** The fit table separates fit validity from evidence of phase, complex-profile or baseline sensitivity. Hover over the QC status to see the measurements and decision thresholds. SE and intervals are conditional on the fit, model and processing.
- **ILT map policies.** Individual, Global and Adaptive Local ILT offer **Exploratory**, **QC-filtered** and **Mono-fit reference** views. The exploratory map can show a reconstructed profile with a QC warning. The 1D overlay identifies unresolved measured signals. Exported PNGs identify the policy used.
- **Auto λ / stability.** Scan the regularization range, inspect the L-curve or noise discrepancy, and track maxima. Apply suggestions explicitly, or retain manual settings. Global ILT and grouped AILT keep μ fixed during the scan. AILT can apply a separate λ to each compatible group.

A reproducible maximum is evidence of stability under the settings tested. Resolution depends on gradient coverage, noise and overlap. The automatic choice of λ can be ambiguous and cannot establish a count of physical populations.

Run **DOSY Lab.cmd** on Windows with Python 3.11 or newer. The same `.proj`, `.list` and `.ints` formats are supported. Existing projects verify their source JDF by checksum and can find a relocated file with the same name and content in the workspace.

[Installation and quick start](../README.md#installation-on-windows) · [Complete manual](INSTRUKCJA.md) · [Scientific methods](SCIENTIFIC-METHODS.md) · [QC validation](QC_AND_MAPS_2.1.md) · [ILT map comparison](ILT_MAP_RESTORATION_2.1.md) · [Regularization benchmarks](ILT_REGULARIZATION_2.1.md) · [Audit history](../audit/README.md) · [Changelog](../CHANGELOG.md)
