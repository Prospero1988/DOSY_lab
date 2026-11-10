"""Generate audit figures. Requires Matplotlib in addition to the release dependencies."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from matplotlib.ticker import MultipleLocator
import numpy as np

from run_audit import B, DUNIT, ILTParameters, decay, local_fixture
from adaptive_ilt import adaptive_local_ilt
from dosy_map import roi_bands
from ilt import solve_ilt

OUT = Path(__file__).resolve().parent / "results"
R = json.loads((OUT/"benchmark.json").read_text(encoding="utf-8"))
D_LABEL = r"$D\ [10^{-10}\ \mathrm{m^2/s}]$"
B_LABEL = r"$b\ [10^9\ \mathrm{s/m^2}]$"
plt.rcParams.update({"font.size":10, "axes.spines.top":False, "axes.spines.right":False,
                     "axes.labelcolor":"#0f172a", "text.color":"#0f172a",
                     "figure.facecolor":"white", "axes.facecolor":"white"})


def save(fig, name):
    fig.savefig(OUT/f"{name}.png",dpi=180,bbox_inches="tight")
    plt.close(fig)


def fit_and_inversion():
    fig, axes = plt.subplots(2,2,figsize=(10.2,7.2),layout="constrained")
    colors = ["#2563a6","#c14b12","#734aa1"]
    for ax, key, title in (
        (axes[0,0],"solvent_tail","A  Solvent tail biases an analyte ROI"),
        (axes[0,1],"gradient_phase_55deg","B  Layer-dependent phase biases D"),
    ):
        row = R["processing"][key]
        x = np.asarray(row["b9"])
        y = np.asarray(row["y"])
        ax.plot(x,y/y[0],"o",ms=3.4,color=colors[0],label="Signed ROI integral")
        d = row["fit"]["D10"]
        r2 = row["fit"]["R2"]
        ax.plot(x,np.asarray(row["prediction"])/y[0],color=colors[1],
                label=f"Mono fit  D={d:.3g}, R²={r2:.5f}")
        ax.plot(x,np.exp(-.38*x),"--",lw=1,color="#172f4d",label="Analyte truth  D=3.8")
        ax.set(title=title,xlabel=B_LABEL,ylabel="Integral / first integral")
        ax.legend(frameon=False,fontsize=8.5)
    for ax, title, ds, amps in (
        (axes[1,0],"C  One true D, regularization-dependent width",[4.],[1.]),
        (axes[1,1],"D  Two true Ds, regularization-dependent count",[2.,16.],[.6,.4]),
    ):
        for d in ds:
            ax.axvline(d,color="#172f4d",ls="--",lw=1.1,label="True D" if d==ds[0] else None)
        for lam,color in zip((1e-8,1e-5,1e-2),colors):
            p = ILTParameters(.25*DUNIT,40*DUNIT,points=100,regularization=lam)
            z = solve_ilt(B*1e9,decay(B,ds,amps),p)
            ax.plot(z.diffusion/DUNIT,z.density,color=color,label=f"λ={lam:g}")
        ax.set(title=title,xlabel=D_LABEL,ylabel="Density per normalized log D",xscale="log",xlim=(.3,39))
        ax.legend(frameon=False,fontsize=8.5)
    save(fig,"fit_and_ilt")


def map_comparison():
    specification = [([3.],[1.],.003,[0]),([3.],[.008],.0009,[0]),
        ([3.25],[1],.016,[0]),([3.25],[.6],.012,[0]),([9.],[.45],.004,[0]),
        ([9.],[.2],.005,[0]),([3.,9.],[.22,.15],0.,[-.032,.032]),([4.],[0.],.001,[0])]
    dataset,rows = local_fixture(specification)
    ailt = adaptive_local_ilt(dataset,rows,15.,10.,1e-6)
    original,_ = roi_bands(dataset.ppm,dataset.spectrum[0],rows[:4])
    grouped,_ = roi_bands(dataset.ppm,dataset.spectrum[0],rows[:4],profiles=ailt.profiles)
    dgrid = np.linspace(2.45,3.8,461)
    ppm = np.linspace(.60,4.40,1300)
    fwhm = .13

    def gaussian(center):
        return np.exp(-4*np.log(2)*((dgrid-center)/fwhm)**2)

    truth = sum(gaussian(ds[0])[:,None]*a[0]*np.exp(-.5*((ppm-i)/.036)**2)[None,:]
                for i,(ds,a,_,_) in enumerate(specification[:4],1))
    def raster(bands):
        result = np.zeros((len(dgrid),len(ppm)))
        for band in bands:
            trace = np.interp(ppm,band.ppm,np.maximum(band.intensity,0.),left=0,right=0)
            result += gaussian(band.diffusion/DUNIT)[:,None]*trace
        return result
    fig, axes = plt.subplots(1,3,figsize=(12.1,3.5),layout="constrained",sharey=True)
    for ax,title,values in zip(axes,("Known D and profiles","Native D of each ROI","Adaptive Local ILT, tolerance 15%"),
                               (truth,raster(original),raster(grouped))):
        im = ax.imshow(values,extent=(ppm[0],ppm[-1],dgrid[0],dgrid[-1]),origin="lower",
                       aspect="auto",cmap="Blues",norm=PowerNorm(.43,vmin=.004,vmax=1))
        ax.axhline(3.,ls="--",lw=.75,color="#c14b12",alpha=.6)
        ax.axhline(3.25,ls="--",lw=.75,color="#c14b12",alpha=.6)
        ax.set(title=title,xlabel=r"$^1$H chemical shift [ppm]",xlim=(4.4,.6))
        ax.xaxis.set_major_locator(MultipleLocator(1.))
    axes[0].set_ylabel(D_LABEL)
    fig.colorbar(im,ax=axes,label="Reference intensity (same scale)",shrink=.8)
    fig.suptitle("Synthetic map detail with common display width 0.13 (audit zoom)")
    save(fig,"map_fidelity")


if __name__=="__main__":
    fit_and_inversion()
    map_comparison()
