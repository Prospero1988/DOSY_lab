# Experimental inverse Laplace analysis

Launch **DOSY Lab.cmd** from the DOSY Lab 2.1.0 directory. Open a JDF or project and complete the usual processing, ROI integration and monoexponential fit. Click **Analyze ILT** in the **ILT** column for the signal you want to explore. The window displays a snapshot of that accepted fit. Reopen it after changing the processing, ROI, fit mask or weights in the main window.

## Controls and plots

Open **Info** in the analysis window for the in-app guide and an explanation of $\lambda$ and the two penalties. The title-bar **×** and **Esc** both close the window.

| Control | Effect |
| --- | --- |
| **λ** | Set the regularization strength using the logarithmic slider or enter an exact value from `1e-12` to `100`. Press Enter to apply typed values. Larger values usually yield smoother or broader solutions and may merge features. |
| **D min / D max** | Set positive grid endpoints in units of $`10^{-10}\ \mathrm{m^2/s}`$. A broad range can introduce regions the acquisition barely constrains. |
| **Grid points / D grid** | Choose 50, 100, 200 or 400 cells, equally spaced in linear D or log D. More cells improve numerical sampling of the distribution. Physical resolution depends on the attenuation data. |
| **Penalty** | Choose zeroth-order Tikhonov or second-derivative smoothness. The latter is the starting setting, with $\lambda=10^{-5}$, 100 log-spaced cells and a range centered geometrically on the monoexponential D from $D/20$ to $20D$. If that D is zero, the initial range is based on the inverse active b span. The center is bounded to $10^{-15}$–$`10^{-5}\ \mathrm{m^2/s}`$ for the controls. |
| **Check ILT sensitivity** | Compare the current result with λ divided and multiplied by ten when in range, and with the other grid coordinate at fixed bounds and cell count. The panel lists maxima, normalized residual RMS and the central 90% mass interval. None of the trials replaces the active result. |
| **Logarithmic D axis** | Toggle the distribution plot axis without solving again or changing the density convention. |
| **Restore defaults / Recalculate ILT** | Restore the starting settings or solve immediately using the displayed inputs. Other fit changes also trigger a solve after a short pause. |

- **A** shows the signed measurements, cached monoexponential curve and ILT reconstruction. Crosses identify points excluded from both fits. Both curves and measured points use the same positive intensity scale.
- **B** shows unit-area density in the normalized grid coordinate. The dashed reference marks the monoexponential D when it lies within the grid range. The summary lists candidate maxima in $\mathrm{m^2/s}$ with their heights relative to the largest maximum.
- **C** shows measured minus predicted residuals for each model. RMS uses included points. With weighted fitting, both weighted RMS values are also shown using the same relative weights.

Drag a rectangle to zoom, enable **Pan** to move the view, and double-click or use **Full range** to reset. Hover over a plot for axis coordinates. Each plot has a fullscreen button. PNG export uses the current view, DOSY Lab's light figure style and a 1600-pixel plot width. The single-panel and combined figures carry 300 dpi metadata. The three-panel export stacks A, B and C. Pending or failed computations can leave the last valid plots visible with a status message. Export is available once the current settings have produced a solution.

## Model and scaling

For the existing DOSY values $b_i$ in $\mathrm{s/m^2}$ and grid centers $D_j$ in $\mathrm{m^2/s}$, the model is $I(b_i) = \sum_j a_j e^{-b_iD_j}$ with $a_j\ge0$. It uses the stored, baseline-processed signed integrals and the same included points as the monoexponential fit. Negative observations contribute signed residuals. At least six distinct included b values are required. When the accepted monoexponential fit was weighted, the inverse-sigma weights retain their relative ratios. The ILT uses the original absolute b values, including when the smallest b is greater than zero.

The intensity scale $I_s$ is the largest absolute included integral, as in the existing fitter. SciPy NNLS solves the augmented system for $f_j=a_j/I_s$, minimizing the mean squared normalized residual plus $`\lambda\|Lf\|^2`$. In weighted mode the residual multipliers are rescaled to unit RMS before taking the mean, so changing the overall scale of the uncertainties does not change the ILT penalty balance. The exponential kernel uses the reciprocal scales $b\times10^{-10}$ and $D/10^{-10}$ to keep the exponent dimensionless.

The grid consists of $N$ equal cells in a dimensionless coordinate $u\in[0,1]$, linear in D or log D. Each weight $f_j$ is a contribution per cell, and $\rho_j=N f_j$ approximates density per unit $u$. The zeroth-order penalty approximates the integral of $\rho^2$. The second-order penalty uses centered second differences in $u$, approximating the integral of $(0.05^2\rho''(u))^2$. This fixed normalized length keeps the curvature penalty on a practical $\lambda$ scale. The plotted density is $\rho_j/\sum f_j$, which integrates to one under the cell convention. Peak heights therefore depend on the grid coordinate and are distinct from integrated population fractions. $\lambda$ is a tuning parameter and its practical effect can change with the grid range, spacing and penalty.

The effective kernel rank counts singular values above $10^{-6}$ of the largest singular value for the weighted data matrix. This threshold is diagnostic and does not estimate the number of diffusion components. Candidate local maxima require at least 5% prominence relative to the highest value, and the dominant bin is reported even at an edge. The software flags a dominant edge peak and strong suppression at the high-D end when present.

## Regularization guidance

**Auto λ / stability** previews a parameter scan without changing the active solution. Select **Apply recommendation** to use a suggested value. [Selection methods and validation](ILT_REGULARIZATION_2.1.md)

## Interpretation and validation

The inverse problem is ill-conditioned. A smaller RMS is expected to be possible with a flexible model and is insufficient evidence for multiple populations. Vary $\lambda$, the grid range and the grid spacing to check which features persist. Inspect signal quality and the b coverage before interpreting a maximum. Close components can be indistinguishable at realistic noise levels. A maximum can reflect regularization, a grid boundary, overlapping NMR signals or processing errors. The width of a regularized peak is not an uncertainty estimate. The displayed SE belongs to the cached monoexponential fit. Automatic $\lambda$ selection and ILT confidence intervals are outside this experimental implementation. See [inversion and regularization for diffusion NMR](https://pubmed.ncbi.nlm.nih.gov/21683632/) and [the effect of S/N on DOSY resolution](https://mr.copernicus.org/articles/2/733/2021/).

Run `.\.venv\Scripts\python.exe tests\run.py` from the application directory after the launcher installs its environment. The numerical tests cover a monoexponential, a well-separated biexponential at nonzero minimum b, noisy data under different penalties, unresolved close components, low intensities, excluded and signed points, and failure paths. The offscreen Qt tests check the table action, the cached data, obsolete results, independent grid and axis controls, fullscreen and PNG exports. These synthetic fixtures establish behavior for their tested b sampling and noise, with no general resolution claim.
