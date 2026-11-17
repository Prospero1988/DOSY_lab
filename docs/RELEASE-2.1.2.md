# DOSY Lab 2.1.2 release notes

## Presaturation baseline correction

- **Solvent tail (presaturation)** fits a dispersive tail together with the polynomial baseline on manually selected peak-free regions. When applied, the correction estimates an independent baseline for every gradient layer. The standard polynomial baseline remains the default.
- **Show presaturation wings** highlights where to select background around the irradiation position. **Hide presaturation wings** clears the orange guides. **Show fitted baseline** displays the model as a yellow dashed line.
- The guide includes a synthetic before-and-after example and instructions for selecting the background. [Presaturation and residual solvent](INSTRUKCJA.md#presaturation-and-residual-solvent)

## ILT controls

- Press **Enter** to apply a typed λ in Individual ILT or a typed λ or μ in Global ILT.
- When a scan yields no clear λ recommendation, acknowledge the message to close the scan and return to the active ILT analysis.

[User guide](INSTRUKCJA.md) · [Changelog](../CHANGELOG.md)
