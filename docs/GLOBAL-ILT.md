# Experimental Global ILT for selected ROIs

Open **DOSY Lab.cmd** from the DOSY Lab 2.1.0 directory. Load and process a measurement, define the integrals and wait for their conventional fits. In the signal table, Ctrl-click or Shift-click at least two fitted rows, then click **Global ILT** next to **Fit / refresh**. The button shows the selected row count. Each included ROI needs six distinct active b values. The window uses a snapshot of the current integrals, masks and fits. Reopen it to capture main-window edits.

For example, select existing ROIs around 7.22, 4.15, 2.84 and 1.32 ppm when those signals occur in the sample. The module reads the integration regions already present in the table. The four coordinates are an example and do not define integration bounds or default ROIs.

## Views and controls

Open **Info** for the in-app instructions and a separate explanation of $\lambda$ and $\mu$. Use the title-bar **×** or **Esc** to close the analysis window.

- **A, Joint distribution** shows the unit-area density formed from the row norms of the joint solution. Markers identify candidate maxima in $\mathrm{m^2/s}$. The distribution summarizes the selected curves mathematically and cannot be interpreted as a molecular concentration distribution.
- **B, ROI distributions** overlays each ROI's solution. The visibility menu affects this view only. The default shows each nonzero distribution at unit area. **Common ROI amplitude scale** restores original fitted amplitudes and divides all densities by the largest fitted area. Dotted lines mark shared candidate positions. **Logarithmic D axis** changes the view and is independent of grid spacing.
- **C, Reconstruction** shows the original signed data and residuals for the chosen ROI. Monoexponential predictions come from the cached fit. Global ILT uses the same b, included points and positive intensity scale. **Compare independent ILT** adds separately solved curves using the same grid, $\lambda$, masks and weighting as Global ILT, with $\mu$ set to zero for the comparison. Curves use absolute b. The displayed y-axis scale is the largest absolute included integral for the selected ROI.
- **D, Component contributions** shows the sum of fitted cell masses in each component interval, defaulting to percent of the total fitted mass for each ROI. Use the checkbox for original integral units at $b=0$. The fitted extrapolated amplitudes depend on spectral intensity, number of contributing nuclei, ROI boundaries and acquisition conditions. They are not calibrated concentrations.
- **E, Results** lists the ROI and its midpoint ppm, conventional D, the largest interval contribution's D, all interval shares in increasing D order, included-point RMS values and solver status. RMS is expressed relative to the ROI's intensity scale.

Use rectangle zoom, **Pan**, double-click or **Full range** to reset, and **Full screen** for each plot. Hover to read coordinates or a heatmap cell value. PNG export uses DOSY Lab's light styling, 1600-pixel panel width and 300 dpi metadata. **Save all plots as PNG** stacks A, B, C (attenuation and residuals for the current ROI), and D. Both the plotted ROI overlay and the exports respect visibility settings.

Enter exact $\lambda$ and $\mu$ values or use the logarithmic sliders. Type `0` for $\mu$ to remove coupling. Set the positive D bounds in $`10^{-10}\ \mathrm{m^2/s}`$ and choose 50, 100, 200 or 400 linear or logarithmic cells. The default is 100 log cells, $\lambda=10^{-6}$, $\mu=10^{-4}$, 10,000 iterations and tolerance $2\times 10^{-6}$. The default range extends from 1/20 of the smallest positive conventional D to 20 times the largest. If no positive conventional D is available, it uses the inverse b span as its center. The center bounds for the controls are $10^{-15}$–$`10^{-5}\ \mathrm{m^2/s}`$. Adjust the range if the relevant attenuations require it. Changes start a background solve after a short delay. Obsolete runs are cancelled cooperatively and their results are discarded. A failed calculation leaves the last result labeled as previous, with export disabled until the current calculation finishes.

## Joint objective and scaling

Let $Y$ have b layers as rows and selected ROIs as columns. For a positive D grid, $K_{ij}=\exp(-b_iD_j)$ using DOSY Lab's existing b values in $\mathrm{s/m^2}$ and D in $\mathrm{m^2/s}$. Let $s_m$ be the largest absolute included integral in ROI $m$, and $Z_{im}=Y_{im}/s_m$. The solver estimates nonnegative cell masses $F_{jm}$ in normalized units from

```math
\min_{F\ge0}\|W\odot(KF-Z)\|_F^2+\lambda P\|F\|_F^2 + \mu\sqrt M\sum_j \|F_{j,:}\|_2.
```

Here $P$ is the number of D cells and $M$ the number of ROIs. The row penalty encourages shared grid positions and allows an individual ROI to contribute zero at any position. Its $\sqrt M$ multiplier balances the sum of the per-ROI losses when identical columns are replicated. Different sets of ROIs can still require a different $\mu$. The raw prediction is $(KF)_{im}s_m$. Excluded points have weight zero and remain visible. Negative measured integrals keep their sign. A positive signal and valid input are checked per ROI.

Without noise weighting, included points have equal weights within each ROI. **Weight by integral SE** uses inverse SE where all included uncertainties for that ROI are valid. ROIs without a complete uncertainty vector use equal point weights. Each column's weights have unit Euclidean norm. Thus they change the relative trust among its points and preserve equal total weight across ROI columns. This does not estimate correlations or impose a cross-ROI SNR ranking. Rescaling an ROI and its uncertainties together leaves the normalized problem invariant.

The cells have width $1/P$ in a normalized coordinate $u\in[0,1]$ linear in D or log D. $P F$ is a density with respect to $u$. The zeroth-order Tikhonov term approximates the integrated squared density and the group sum approximates its integrated row norm, so changing P does not by itself multiply the penalty strength. The A plot shows $P G/\sum_jG_j$ with $`G_j=\|F_{j,:}\|_2`$. Local maxima require at least 5% prominence relative to the largest row norm. Minima between neighboring maxima divide the grid into intervals. The contribution matrix sums the original-unit cell masses in those intervals. Changing the regularization can change both the number of intervals and their boundaries.

## Optimization and limits

The solver uses accelerated proximal gradient with adaptive restart. It takes a fixed step from the largest eigenvalue of the weighted kernel Gram matrix across ROIs and the Tikhonov term. The proximal operator projects onto the positive orthant and then applies $L_2$ shrinkage to each diffusion row. Every 20 iterations it checks the Frobenius norm of the proximal gradient mapping divided by the norm of the smooth gradient at zero. The result records this stationarity residual, the objective, the iteration count and whether the tolerance was met. An iterate stopped at the limit is labeled approximate. A small stationarity residual establishes numerical convergence of the selected objective. It does not establish physical resolution. [FISTA method](https://doi.org/10.1137/080716542)

Weak signal, sparse sampling of early or late b, correlated noise, imperfect phase or baseline, overlap and exchange kinetics can change the inferred maxima. A shared D is insufficient to identify a molecule. Exchange and model mismatch may produce apparent diffusion behavior. The window flags edge maxima and strong suppression at the high-D end. Global ILT here uses integrated ROIs. It does not reconstruct a chemical-shift-resolved DOSY 2D spectrum. [Simultaneous Laplace inversion for DOSY](https://doi.org/10.1016/j.jmr.2017.03.004)

## Regularization guidance

**Auto λ / stability** previews a parameter scan without changing the active solution. Global ILT keeps μ fixed during the scan. Select **Apply recommendation** to use a suggested value. [Selection methods and validation](ILT_REGULARIZATION_2.1.md)

## Reproducible validation

Run `.\.venv\Scripts\python.exe tests\run.py` from the application directory. The synthetic fixtures in `tests/test_global_ilt.py` use 65 b values, including zero and dense sampling at low b up to a maximum of $`35\times 10^9\ \mathrm{s/m^2}`$. The reported maxima use a 100-cell log grid from `0.3` to `40` in $`10^{-10}\ \mathrm{m^2/s}`$ and a solver tolerance of $10^{-6}$. The fixtures test amplitude and noise differences, two shared populations, partial ROI membership, close components, parameter sensitivity, rescaling, signed/masked data, cancellation and the iteration limit. At $\mu=0$ the joint solver reproduces the separable NNLS predictions to within $1.5\times 10^{-4}$ in the tested signal units. The Qt tests use existing integral rows and verify selection, comparison, plot export, cancellation and preservation of the cached monoexponential fits.

In the partly shared noiseless example the input D values `1.5, 7, 24` produced maxima `1.47, 7.04, 25.13`, in units of $`10^{-10}\ \mathrm{m^2/s}`$. Under additive noise with SD `0.004`, the tested sequence of $(\lambda,\mu)$ settings $(10^{-8},10^{-6})$, $(10^{-6},10^{-4})$, $(10^{-4},0.003)$ and $(10^{-3},0.02)$ produced 2, 3, 2 and 1 maxima, respectively. The fastest component is therefore parameter-sensitive in this fixture. At noise SD `0.003`, the tested pair `4.0, 4.5` produced one maximum near `4.32`. In four noisy single-D realizations, matched independent and Global ILT both localized the dominant D near the true value. The largest difference in original-unit residual RMS between methods was below $4\times 10^{-4}$. These fixtures do not establish general superiority or resolution limits for experimental data.
