# DOSY Lab 2.0.0 scientific methods

## Measurements and diffusion weighting

Each gradient layer contains a complex direct FID. The application apodizes and optionally zero-fills the FID, Fourier transforms it, then applies the common phase correction and any supported digital filter correction. Integration uses the real, signed spectrum after optional baseline subtraction. The baseline is estimated separately on each layer when subtraction is applied. Integral endpoints are interpolated and the area is calculated by the trapezoidal rule on the ppm axis. An ROI (region of interest) therefore yields a signed intensity sequence $I(b_i)$ in arbitrary intensity × ppm. Display rescaling of integrals leaves the fit inputs unchanged. [Processing and fit controls](INSTRUKCJA.md)

The diffusion weighting $b_i$ is in $\mathrm{s/m^2}$ and the translational diffusion coefficient $D$ is in $\mathrm{m^2/s}$. For an ideal rectangular gradient pair with duration $\delta$, separation $\Delta$, amplitude $g$ and gyromagnetic ratio $\gamma$, the diffusion weighting is

```math
b=(\gamma g\delta)^2(\Delta-\delta/3)
```

The actual weighting depends on the gradient waveform and timing. DOSY Lab uses an explicit JDF b axis, a supported and verified BPP-STE-LED or BPP-DSTE-LED encoding model, or a manually supplied list in acquisition order. This calculation relies on the acquisition metadata and calibration. An incorrect b scale biases every inferred D in the same experiment. Repeated b values remain separate observations. A conflict among recorded timings or between comparable explicit and modeled b values disables automatic diffusion fitting. The FIDs remain available and a documented manual b list can be supplied.

## Conventional monoexponential fitting

The Stejskal–Tanner attenuation model for one freely diffusing contribution is

```math
I(b)=I_0\exp(-bD)
```

$I_0$ is the extrapolated integral at zero weighting. A fitted D may be apparent when the ROI contains overlapping signals, exchange or other departures from the model.

The fit uses bounded nonlinear least squares with $I_0\ge0$ and $D\ge0$. It uses at least three distinct included b values. The fitter uses the smallest active b as its amplitude reference and scales b by the active span. It divides intensities by the largest absolute included intensity. A profiled amplitude supplies the initial estimate before bounded nonlinear least squares. The amplitude at zero weighting is extrapolated when representable. Optional positive point standard errors $\sigma_i$ weight the residuals by $1/\sigma_i$. Negative observations keep their sign. Excluded observations remain in the project. A fit's $\mathrm{SE}(D)$ comes from the local parameter covariance and is conditional on the model, data and weights. Across selected ROIs, SD describes the dispersion of fitted D values and SEM is descriptive. Errors may be correlated across ROIs in one experiment.

## Measurement-quality diagnostics

Each active layer has a robust spectral noise estimate from the quieter half of detrended blocks in both quadratures. The complex ROI diagnostic uses a weighted rank-one template and requires projected S/N of at least five in at least eight distinct active b values. Layers below that criterion are excluded only from this diagnostic. A common phase offset is removed from the projected phase traces. Relative phase and its estimated noise scale are reported separately. Random whole-layer sign changes preserve within-layer correlation across ppm. This b-trend screen assumes roughly symmetric, independent layer errors and estimates its template from the same data. It can miss nonmonotonic changes or overlapping species with indistinguishable attenuation.

QC distinguishes model validity from measurement diagnostics. A detected trend is informational unless an alternative reading changes D beyond a screening scale of $1.96(\mathrm{SE}_1+\mathrm{SE}_2)$. Shape is checked by comparing two portions of the complex profile, phase by comparing the signed integrals with template projection magnitudes. For baseline subtraction, the comparison uses identical ROI bounds, b values, fit mask and weights with and without the correction. The fit detail also reports baseline correction RMS relative to the uncorrected integrals. These sensitivity comparisons are heuristic and their screening scale is not a confidence interval for the difference. The original signed fit inputs and fitted values are preserved. [Benchmark results](QC_AND_MAPS_2.1.md)

## Individual inverse Laplace transform

For a continuous distribution, the diffusion model is $`I(b)=\int a(D)e^{-bD}\,dD`$. DOSY Lab discretizes it using the recorded $b_i$ and candidate diffusion grid centers $D_j$. The kernel is $K_{ij}=\exp(-b_i D_j)$. The absolute recorded b values enter the exponent, even when none is zero. Finite and noisy measurements constrain only a limited number of directions in this kernel. Very different distributions can explain the same attenuation within the noise, so the inversion needs regularization.

For $z_i=I(b_i)/s$ with $s$ the largest absolute included integral, the solver finds the nonnegative cell masses $f$ by augmented nonnegative least squares.

```math
\min_{f\ge0}\ \|W(Kf-z)\|_2^2+\lambda\|Lf\|_2^2.
```

The weights have unit Euclidean norm across included points. In the unweighted case this yields mean squared residuals in normalized units. In weighted mode it preserves the relative inverse-SE weights of the conventional fit. The reconstruction is $sKf$. Each ROI needs six distinct included b values.

The $P$ equal cells are defined in a dimensionless coordinate $u\in[0,1]$, linear in D or log D. The zeroth-order penalty uses $`L=\sqrt{P}\,\mathrm{Id}`$. The smoothness penalty uses $L=\sqrt P(0.05P)^2\Delta^2$, where $\Delta^2$ is the centered second-difference operator on the cell masses. The plotted unit-area density in $u$ is $Pf_j/\sum_k f_k$. Its height is distinct from cell mass and depends on the chosen coordinate. Increasing $\lambda$ can suppress or merge maxima. [Settings and diagnostics](ILT.md) · [Diffusion NMR inversion and regularization](https://pubmed.ncbi.nlm.nih.gov/21683632/)

## Global ILT

The joint solver uses $M$ selected ROI attenuation columns. For each ROI $m$, $s_m$ is the largest absolute included integral and $Z_{im}=I_m(b_i)/s_m$. It estimates the nonnegative cell-mass matrix $F$ using the objective

```math
\min_{F\ge0}\ \|W\odot(KF-Z)\|_F^2+\lambda P\|F\|_F^2
+\mu\sqrt M \sum_j\|F_{j,:}\|_2 .
```

Included points have unit-norm weights per ROI, optionally using a complete valid inverse-SE vector. Excluded points have weight zero. The ridge term with $\lambda$ penalizes squared density in the normalized coordinate. The row penalty with $\mu$ favors a common set of D cells while allowing zero contributions in any individual ROI. The $\sqrt M$ factor scales the penalty relative to replicated data columns. The practical effect of $\mu$ can still change with ROI selection. The normalized problem is invariant to a common positive rescaling of the measured intensities and SE values within one ROI. Multiplication by $s_m$ restores its original signal units.

The accelerated proximal gradient solver projects onto nonnegative masses and applies L2 shrinkage row by row. The step uses the largest weighted Gram-matrix eigenvalue plus the regularization term, and adaptive restart limits momentum. A normalized proximal-gradient mapping tests stationarity every 20 iterations. The reported objective, iteration count and stationarity measure describe the numerical solve. [Global ILT controls and validation](GLOBAL-ILT.md) · [Accelerated proximal gradient method](https://doi.org/10.1137/080716542)

The global summary is formed from $`G_j=\|F_{j,:}\|_2`$ and plotted as $PG_j/\sum_kG_k$. Candidate maxima use a prominence threshold relative to the largest maximum, and grid edges are flagged. Minima between selected maxima define component intervals. The contribution matrix sums the cell masses in each interval, normalized within each ROI for percent values or multiplied by $s_m$ for original integral units at $b=0$. Both are model-derived signal contributions. Molecular concentration would require independent calibration and signal assignment.

## Adaptive Local ILT and the 2D display

Adaptive Local ILT uses the cached conventional fits to propose groups. It screens for positive D, suitable active points and fit quality, then bounds the symmetric distance between every pair of candidate members by the chosen tolerance (5% by default). For two estimates the distance is $`200|D_1-D_2|/(D_1+D_2)`$ percent. A shared monoexponential fit screens out candidates whose normalized residual RMS rises above both 1.5 times their independent RMS and that RMS plus 0.005. A joint local inversion requires at least six distinct active b values, a positive reference integral and an estimated reference S/N of at least five. When a reference SE is unavailable, residual RMS supplies a screening proxy.

Within a group each ROI is still independently normalized. The default amplitude-balance threshold is a 10-fold difference from the strongest reference integral. Above this ratio a weaker ROI can receive less influence according to its estimated S/N. The bounded factor multiplies its contribution to the residual objective. Each group uses its own joint inversion with a 100-cell logarithmic grid, bounds from the smallest conventional D divided by 20 to the largest multiplied by 20, default settings $\lambda=10^{-6}$ and $\mu=10^{-4}$.

The exploratory policy offers an independent ILT for ungrouped ROIs with valid attenuation inputs. For the QC-filtered policy, an eligible ungrouped ROI can be split using nearby $\lambda$ settings, consistent interior maxima and improved residual RMS. Accepted multicomponent profiles are estimated by nonnegative fits to the measured spectra at each ppm bin using decay bases fixed by the recovered component intervals. The map screens for ill-conditioned bases, reconstruction residuals and inconsistent amplitudes. Unresolved ROIs remain visible in the reference-spectrum overlay by default and contribute no reconstructed ILT intensity. Bins shared by overlapping ROI definitions are masked. [Grouping, diagnostics and reconstruction criteria](DOSY-ALIGNMENT.md)

The original map places each positive real reference profile at its conventional D. Display alignment can move an existing profile to a representative group D. Shared monoexponential fitting estimates one D per compatible group. Individual, Global and Adaptive Local ILT display reconstructed component profiles within the selected ROIs. In the exploratory policy, nonnegative fitted profiles are capped by the ILT predicted area of each component and by the measured positive profile at each ppm bin. QC diagnostics remain visible regardless of the selected policy. The Gaussian D-axis display width is adjustable. The optional D projection sums positive displayed map intensity across ppm before applying visual gain and threshold.

## Interpreting results

| Output | Meaning |
| --- | --- |
| Experimental attenuation | Signed integrals measured across gradient layers after processing. |
| Conventional and shared D | Point estimates conditional on a monoexponential model and the selected b values. |
| ILT distributions | Regularized estimates conditional on the D grid, penalties, weights and ROI selection. |
| Aligned and reconstructed maps | ROI-based visualizations that also depend on grouping, reconstruction gates and display settings. |

Maxima can change in number, position and width when regularization or the D grid changes. Peak width does not provide a confidence interval. Similar D values alone cannot establish molecular identity, and a local maximum is insufficient evidence for a separate population. Global and adaptive analyses add assumptions about shared positions or the compatibility of groups. Cleaner bands can result from those assumptions or from visual settings. Inspect signed attenuation, line shapes, residuals, b coverage and sensitivity before interpreting the map. Overlap, S/N, gradient sampling, exchange kinetics and model selection all affect the inference. Validation fixtures demonstrate behavior for their chosen conditions. They do not establish a general resolution limit or experimentally verify molecular assignments. [Effect of S/N on DOSY resolution](https://mr.copernicus.org/articles/2/733/2021/)
