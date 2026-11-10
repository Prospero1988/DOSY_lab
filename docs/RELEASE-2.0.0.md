# DOSY Lab 2.0.0 release notes

The conventional DOSY workflow carries over from 1.x. Open a raw JEOL JDF or an existing project, process the spectrum, choose the integration ranges and fit the signed attenuation. A saved `.ints` template lets you use the same ranges, IDs and groups on another referenced measurement.

- **Individual ILT** opens from a fitted ROI. Compare the measured attenuation, cached monoexponential prediction and regularized reconstruction. Explore how the grid and penalty affect the distribution.
- **Global ILT** jointly analyzes selected ROIs. Its shared-position penalty allows each signal to contribute independently. Use the reconstruction and contribution matrix to see which ROIs support each candidate maximum.
- **2D DOSY map** starts at the independent fitted D for each ROI. Compare it with display alignment, shared monoexponential fits, Global ILT profiles or Adaptive Local ILT. The map now has an optional linked diffusion projection. Use **Z** to zoom, **P** to pan, **F** to reset and **Delete** to clear the guides.

All ILT analyses are exploratory. Similar diffusion estimates and persistent maxima can guide further checks, but do not identify molecules or establish how many populations are present. Map modes use selected ROIs and apply different assumptions about sharing and grouping. They leave the conventional fit table intact.

Extract the source archive and run **DOSY Lab.cmd** on Windows with Python 3.11 or later. Existing `.proj` and `.list` files keep their formats. See the [installation guide](README.md#installation-on-windows), [user manual](INSTRUKCJA.md), [scientific methods](SCIENTIFIC-METHODS.md) and [changelog](CHANGELOG.md).
