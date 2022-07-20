"""Plotting functions"""

import numpy as np
import matplotlib.pyplot as plt

from piff_syst_plots.funcs import compute_res, bin_res_by_color
from piff_syst_plots.core import Plotter


class PlotStarsPerCCD(Plotter):
    """Rough histogram of nstars/CCD

    The 'factor' parameter accounts the fact that we are making this plot with
    the reserve stars, which are typically 20% of the total,  Thuse the default value of 5.
    """

    default_config = dict(
        nbins=70,
        factor=5.0,
        figsize=(12, 8),
    )

    def _make_plot(self, cat):
        fig = plt.figure(figsize=self.config["figsize"])
        _, nstars = np.unique(
            np.stack((cat["EXPNUM"], cat["CCDNUM"]), axis=1), axis=0, return_counts=True
        )
        plt.hist(nstars * self.config["factor"], bins=self.config["nbins"])
        plt.xlabel("Nstars / CCD")
        return fig


class PlotFootprint(Plotter):
    """Make a hexbin plot of the catalog footprint

    This isn't a pretty plot, it is just for sanity checking
    """

    default_config = dict(
        figsize=(15, 12),
    )

    def _make_plot(self, cat):
        fig, ax = plt.subplots(1, 1, figsize=self.config["figsize"])

        catra = cat["RA"]
        catra[catra > 180.0] -= 360.0
        im = ax.hexbin(catra, cat["DEC"], bins="log", mincnt=1)
        ax.set_xlabel("RA")
        ax.set_ylabel("DEC")
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("# stars")
        ax.invert_xaxis()
        # ax.set_ylim(cat['dec'].min(), cat['dec'].max())
        return fig


class PlotFluxByBand(Plotter):
    """Makes a single plot with all the fluxes"""

    default_config = dict(
        zeropt=30.0,
        figsize=(16, 22),
        bands="grizY",
        nbins=50,
        xmin=14,
        xmax=24,
    )

    def _make_plot(self, cat):
        fig, ax = plt.subplots(1, 1, figsize=self.config["figsize"])
        bins = np.linspace(
            self.config["nbins"], self.config["xmin"], self.config["xmax"]
        )
        for band in self.config["bands"]:
            flux = cat["FLUX"][cat["BAND"] == band]
            zeropt = 30.0
            mag = zeropt - 2.5 * np.log10(flux)
            print(band, min(mag), max(mag))
            ax.hist(mag, bins=bins, histtype="step", label=band, lw=3)
        ax.set_xlabel("mag")
        ax.legend()
        return fig


class PlotColorsByBandMulti(Plotter):
    """plot color distribution per band on multiple plots"""

    default_config = dict(
        figsize=(16, 22),
        bands="grizY",
    )

    def _make_plot(self, cat):
        bands = self.config["bands"]
        fig, axs = plt.subplots(len(bands), 1, figsize=self.config["figsize"])

        for i, band in enumerate(bands):
            if band == "z":
                color = "IZ_COLOR"
            else:
                color = "GI_COLOR"
            chist = cat[color][cat["BAND"] == band]
            ax = axs[i]
            ax.hist(chist, bins=50)
            ax.set_xlabel(color)
        return fig


class PlotColorsByBandSingle(Plotter):
    """plot color distribution per band on single plot"""

    default_config = dict(
        figsize=(12, 8),
        bands="grizY",
    )

    def _make_plots(self, cat):

        bands = self.config["bands"]
        fig, ax = plt.subplots(1, 1, figsize=self.config["figsize"])
        bins = np.linspace(0, 3.5, 50)
        for band in bands:
            if band == "z":
                color = "IZ_COLOR"
            else:
                color = "GI_COLOR"
            chist = cat[color][cat["BAND"] == band]
            ax.hist(chist, bins=bins, histtype="step", label=band, lw=3)
        ax.legend()
        ax.set_xlim(0, 3.5)
        ax.set_xlabel("g-i or i-z")
        return fig


class PlotColorVMagByBand(Plotter):
    """plot color v mag distribution per band on a grid
    of 2 X nband plots"""

    default_config = dict(
        figsize=(16, 24),
        bands="grizY",
    )

    def _make_plots(self, cat):

        bands = self.config["bands"]
        fig, axs = plt.subplots(len(bands), 2, figsize=self.config["figsize"])
        cmap = plt.get_cmap("jet")
        for i, band in enumerate(bands):
            ax = axs[i][0]
            # plot g-r_color vs g_mag 2d histogram with log density
            if band == "z":
                color = "IZ_COLOR"
            else:
                color = "GI_COLOR"
            chist = cat[color][cat["BAND"] == band]
            flux = cat["FLUX"][cat["BAND"] == band]
            zeropt = 30.0
            mag = zeropt - 2.5 * np.log10(flux)
            im = ax.hexbin(chist, mag, bins="log", mincnt=1, cmap=cmap)
            ax.set_ylim(14, 30)
            ax.set_ylabel("{} mag".format(band))
            ax.set_xlabel(color)
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label("# stars")

            # plot color histogram in mag slices
            ax = axs[i][1]
            nbins = 7
            cbins = np.linspace(15.0, 21.0, nbins)
            for j in range(len(cbins) - 1):
                k = j + 1
                mask = np.logical_and(
                    np.greater_equal(mag, cbins[i]), np.less(mag, cbins[k])
                )
                color_slice = chist[mask]
                loedge = "{0:3.1f}".format(cbins[i])
                hiedge = "{0:3.1f}".format(cbins[k])
            ax.hist(
                color_slice,
                density=True,
                label="{} < {} mag < {}".format(loedge, band, hiedge),
                histtype="step",
                lw=2,
                color=cmap(k / nbins),
            )
            ax.legend(ncol=1, fontsize=11)
            ax.set_xlabel(color)
            ax.set_ylabel("Normalized count")
        return fig


class PlotSizeAndEllipticityByBand(Plotter):
    """Plot the size and the ellipiticy on
    a 2 X nband grid of plots"""

    default_config = dict(
        figsize=(16, 24),
        bands="grizY",
    )

    def _make_plots(self, cat):
        bands = self.config["bands"]
        fig, axs = plt.subplots(len(bands), 2, figsize=self.config["figsize"])
        for i, band in enumerate(bands):
            ax = axs[i][0]
            # plot size histogram
            size = cat["T_DATA"][cat["BAND"] == band]
            ax.hist(size, bins=np.linspace(0, 2.5, 50))
            ax.set_xlim(0, 2.5)
            ax.set_xlabel("T [arcsec^2]")
            ax.set_ylabel(band)

            # plot e1 and e2 histograms
            ax = axs[i][1]
            e1 = cat["G1_DATA"][cat["BAND"] == band]
            e2 = cat["G2_DATA"][cat["BAND"] == band]
            ax.hist(
                e1, bins=np.linspace(-0.3, 0.3, 50), histtype="step", label="e1", lw=2
            )
            ax.hist(
                e2, bins=np.linspace(-0.3, 0.3, 50), histtype="step", label="e2", lw=2
            )
            ax.set_xlim(-0.3, 0.3)
            ax.legend()
            ax.set_xlabel("Ellipticity")
        return fig


class PlotSizeByBand(Plotter):
    """Plot the distributions of sizes for each band
    on single plot"""

    default_config = dict(
        figsize=(16, 24),
        bands="grizY",
        logscale=False,
    )

    def _make_plots(self, cat):
        bands = self.config["bands"]
        fig, axs = plt.subplots(1, 1, figsize=self.config["figsize"])
        for band in bands:
            ax = axs
            # plot size histogram
            size = cat["T_DATA"][cat["BAND"] == band]
            print(
                "{} mean, median, std: {}".format(
                    band, (np.mean(size), np.median(size), np.std(size))
                )
            )
            ax.hist(size, np.linspace(0, 2.5, 50), label=band, histtype="step", lw=3)
        ax.set_xlim(0, 2.5)
        ax.legend()
        ax.set_xlabel("T [arcsec^2]")
        return fig


class PlotSeeingByBand(Plotter):
    """Plot the distributions of seeing for each band
    on single plot"""

    default_config = dict(
        figsize=(16, 24),
        bands="grizY",
    )

    def _make_plots(self, cat):

        bands = self.config["bands"]
        fig, axs = plt.subplots(1, 1, figsize=self.config["figsize"])
        for band in bands:
            ax = axs
            # plot size histogram
            size = cat["PSF_FWHM"][cat["BAND"] == band]
            ax.hist(size, np.linspace(0, 2.5, 50), label=band, histtype="step", lw=3)
        ax.set_xlim(0, 2.5)
        ax.legend()
        ax.set_xlabel("Exposure median FWHM [arcsec]")
        return fig


class PlotTversusFWHMByBand(Plotter):
    """Plot T v. FWHM by band on a single plot"""

    default_config = dict(
        figsize=(20, 18),
        bands="griz",
    )

    def _make_plot(self, cat):
        bands = self.config["bands"]
        fig, axs = plt.subplots(2, 2, figsize=self.config["figsize"])
        for i, band in enumerate(bands):
            # <T> vs (color, mag)
            ax = axs[int(i / 2)][i % 2]
            # if band == 'z':
            #     color = 'IZ_COLOR'
            #     cmax = 0.7
            # else:
            #     color = 'GI_COLOR'
            #     cmax = 3.5
            data = cat[cat["BAND"] == band]
            tsize = data["T_DATA"]
            fwhm = data["PSF_FWHM"]

            cmap = plt.get_cmap("jet")
            im = ax.hexbin(
                tsize,
                fwhm,
                mincnt=1,
                bins="log",
                cmap=cmap,
                extent=(0, 3.5, 0.5, 3.0),
                vmin=1.0,
                vmax=1e6,
            )
            ax.plot(np.linspace(0, 3.5, 2), np.linspace(0.5, 3.0, 2), c="k", alpha=0.3)
        # ax.set_aspect('equal')
        ax.set_xlim(0, 3.5)
        ax.set_ylim(0.5, 3.0)
        ax.set_xlabel(r"$T [arcsec^2]$")
        ax.set_ylabel(r"$FWHM [arcsec]$")
        _ = plt.colorbar(im, ax=ax)
        # cbar.set_label(r'$\left<T_{\mathrm{PSF}}\right>$')
        return fig


def add_profile_by_color_to_ax(
    m, dT, ax, color, cmin=0.0, cmax=3.5, tmin=0.0, tmax=1.0  # min_mused=None
):
    """Add a profile plot to an axis"""
    mag_bins = np.linspace(cmin, cmax, 30)
    print("col_bins = ", mag_bins)

    index = np.digitize(m, mag_bins)
    bin_dT = [dT[index == i].mean() for i in range(1, len(mag_bins))]
    print("bin_dT = ", bin_dT)
    bin_dT_err = [
        np.sqrt(dT[index == i].var() / len(dT[index == i]))
        for i in range(1, len(mag_bins))
    ]
    print("bin_dT_err = ", bin_dT_err)

    # Fix up nans
    for i in range(1, len(mag_bins)):
        if i not in index:
            bin_dT[i - 1] = 0.0
            bin_dT_err[i - 1] = 0.0

    ax.set_ylim(tmin, tmax)
    # ax.plot([cmin,cmax], [0,0], color='black')
    # ax.plot([min_mused,min_mused],[-1,1], color='Grey')
    # ax.fill( [min_mag,min_mag,min_mused,min_mused], [-1,1,1,-1], fill=True, color='Grey',alpha=0.3)
    _ = ax.errorbar(mag_bins[:-1], bin_dT, yerr=bin_dT_err, color="b", fmt="o")
    # ax.axhline(y=0.003, linewidth=4, color='grey')
    # ax.legend([t_line], [r'$\delta T$'])
    ax.set_ylabel(r"$\left<T_{data}\right> \quad({\rm arcsec}^2)$")
    ax.set_xlabel(color)
    plt.tight_layout()


class PlotSizeVColorByBands(Plotter):
    """Plot the size v. color in a
    2 X nband grid of plots"""

    default_config = dict(
        figsize=(16, 24),
        bands="griz",
    )

    def _make_plots(self, cat):
        bands = self.config["bands"]
        fig, axs = plt.subplots(4, 2, figsize=self.config["figsize"])
        for i, band in enumerate(bands):
            ax = axs[i][0]
            # plot size vs color 2d hist with log density
            if band == "z":
                color = "IZ_COLOR"
            else:
                color = "GI_COLOR"
            # cat = cat[cat['T_DATA'] < 2.5]
            chist = cat[color][cat["BAND"] == band]
            size = cat["T_DATA"][cat["BAND"] == band]
            cmap = plt.get_cmap("jet")
            im = ax.hexbin(chist, size, bins="log", mincnt=1, cmap=cmap)
            ax.set_ylim(0.0, 1.6)
            ax.set_ylabel("{} band: T [arcsec^2]".format(band))
            ax.set_xlabel(color)
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label("# stars")

            # plot size vs color profile
            ax = axs[i][1]
            if band == "z":
                cmax = 0.7
            else:
                cmax = 3.5
            add_profile_by_color_to_ax(chist, size, ax, color, cmax=cmax)
        return fig


def plot_residual_on_axis(
    band,
    color,
    min_edge,
    max_edge,
    bins,
    bin_dT,
    bin_dT_err,
    bin_dTfrac,
    bin_dTfrac_err,
    bin_de1,
    bin_de1_err,
    bin_de2,
    bin_de2_err,
    cutlim=None,
    Tlims=(-0.0075, 0.0075),
    Tfraclims=(-0.02, 0.02),
    elims=(-6.0e-4, 6.0e-4),
    label=None,
    axes=None,
):
    """Draw a residual plot on an axis"""
    ax = axes[0]
    ax.set_title("{} band".format(band))
    ax.set_ylim(Tlims)
    ax.plot([min_edge, max_edge], [0, 0], color="black")
    if cutlim is not None:
        ax.plot([cutlim, cutlim], [-1, 1], color="Grey")
        ax.fill(
            [min_edge, min_edge, cutlim, cutlim],
            [-1, 1, 1, -1],
            fill=True,
            color="Grey",
            alpha=0.3,
        )
    _ = ax.errorbar(bins[:-1], bin_dT, yerr=bin_dT_err, fmt="o", label=label)
    # ax.axhline(y=0.003, linewidth=4, color='grey')
    # ax.legend([t_line], [r'$\delta T$'])
    ax.set_ylabel(r"$(T_{\rm PSF} - T_{\rm model}) \quad({\rm arcsec}^2)$")

    _ = bin_de2, bin_de2_err

    ax = axes[1]
    ax.set_ylim(Tfraclims)
    ax.plot([min_edge, max_edge], [0, 0], color="black")
    if cutlim is not None:
        ax.plot([cutlim, cutlim], [-1, 1], color="Grey")
        ax.fill(
            [min_edge, min_edge, cutlim, cutlim],
            [-1, 1, 1, -1],
            fill=True,
            color="Grey",
            alpha=0.3,
        )
    _ = ax.errorbar(bins[:-1], bin_dTfrac, yerr=bin_dTfrac_err, label=label, fmt="o")
    # ax.legend([t_line], [r'$\delta T$'])
    ax.set_ylabel(r"$(T_{\rm PSF} - T_{\rm model})/ T_{\rm PSF}$")

    ax = axes[2]
    ax.set_ylim(elims)
    ax.plot([min_edge, max_edge], [0, 0], color="black")
    if cutlim is not None:
        ax.plot([cutlim, cutlim], [-1, 1], color="Grey")
        ax.fill(
            [min_edge, min_edge, cutlim, cutlim],
            [-1, 1, 1, -1],
            fill=True,
            color="Grey",
            alpha=0.3,
        )
    _ = ax.errorbar(bins[:-1], bin_de1, yerr=bin_de1_err, label=label, fmt="o")
    # e2_line = ax.errorbar(bins[:-1], bin_de2, yerr=bin_de2_err, label=label, fmt='o')
    # ax.axhline(y=0.0002, linewidth=4, color='grey')
    # ax.axhline(y=-0.0002, linewidth=4, color='grey')
    # ax.legend()
    # ax.legend([e1_line, e2_line], [r'$e_1$', r'$e_2$'])
    ax.set_ylabel(r"$e_{\rm PSF} - e_{\rm model}$")
    ax.set_xlim(min_edge, max_edge)
    ax.set_xlabel(color)
    plt.tight_layout()


class PlotResidualsByBand(Plotter):
    """Plot the residuals by band on an
    3 x nBand grid"""

    default_config = dict(
        figsize=(20, 12),
        bands="griz",
    )

    def _make_plots(self, cat):
        bands = self.config["bands"]
        fig, axes = plt.subplots(
            3, 4, figsize=self.config["figsize"], sharey="row", sharex="col"
        )
        for band in bands:
            data = cat[cat["BAND"] == band]
            fracsizeres, sizeres, e1res, e2res = compute_res(data)
            if band == "z":
                color = "IZ_COLOR"
            else:
                color = "GI_COLOR"
            (
                min_edge,
                max_edge,
                bins,
                bin_dT,
                bin_dT_err,
                bin_dTfrac,
                bin_dTfrac_err,
                bin_de1,
                bin_de1_err,
                bin_de2,
                bin_de2_err,
            ) = bin_res_by_color(data[color], sizeres, fracsizeres, e1res, e2res)
            plot_residual_on_axis(
                band,
                color,
                min_edge,
                max_edge,
                bins,
                bin_dT,
                bin_dT_err,
                bin_dTfrac,
                bin_dTfrac_err,
                bin_de1,
                bin_de1_err,
                bin_de2,
                bin_de2_err,
                Tlims=None,
                Tfraclims=None,
                elims=None,
                axes=axes,
            )
        return fig


def plot_bin_by_mag_on_axis(
    band,
    m,
    dT,
    dTfrac,
    de1,
    de2,
    mmin=15.0,
    mmax=21.0,
    min_mused=16.5,
    label=None,
    axes=None,
    xlabel="Magnitude",
):
    """Plot something"""
    min_mag = mmin
    max_mag = mmax
    print(band, min_mag, max_mag)
    mag_bins = np.linspace(min_mag, max_mag, 30)
    # print('col_bins = ',mag_bins)

    index = np.digitize(m, mag_bins)
    # print('len(index) = ',len(index))
    bin_de1 = [de1[index == i].mean() for i in range(1, len(mag_bins))]
    # print('bin_de1 = ',bin_de1)
    bin_de2 = [de2[index == i].mean() for i in range(1, len(mag_bins))]
    # print('bin_de2 = ',bin_de2)
    bin_dT = [dT[index == i].mean() for i in range(1, len(mag_bins))]
    # print('bin_dT = ',bin_dT)
    bin_dTfrac = [dTfrac[index == i].mean() for i in range(1, len(mag_bins))]
    # print('bin_dTfrac = ',bin_dTfrac)
    bin_de1_err = [
        np.sqrt(de1[index == i].var() / len(de1[index == i]))
        for i in range(1, len(mag_bins))
    ]
    # print('bin_de1_err = ',bin_de1_err)
    bin_de2_err = [
        np.sqrt(de2[index == i].var() / len(de2[index == i]))
        for i in range(1, len(mag_bins))
    ]
    # print('bin_de2_err = ',bin_de2_err)
    bin_dT_err = [
        np.sqrt(dT[index == i].var() / len(dT[index == i]))
        for i in range(1, len(mag_bins))
    ]
    # print('bin_dT_err = ',bin_dT_err)
    bin_dTfrac_err = [
        np.sqrt(dTfrac[index == i].var() / len(dTfrac[index == i]))
        for i in range(1, len(mag_bins))
    ]
    # print('bin_dTfrac_err = ',bin_dTfrac_err)

    # Fix up nans
    for i in range(1, len(mag_bins)):
        if i not in index:
            bin_de1[i - 1] = 0.0
            bin_de2[i - 1] = 0.0
            bin_dT[i - 1] = 0.0
            bin_dTfrac[i - 1] = 0.0
            bin_de1_err[i - 1] = 0.0
            bin_de2_err[i - 1] = 0.0
            bin_dT_err[i - 1] = 0.0
            bin_dTfrac_err[i - 1] = 0.0

    ax = axes[0]
    ax.set_title("{} band".format(band))
    ax.set_ylim(-0.005, 0.005)
    ax.plot([min_mag, max_mag], [0, 0], color="black")
    ax.plot([min_mused, min_mused], [-1, 1], color="Grey")
    ax.fill(
        [min_mag, min_mag, min_mused, min_mused],
        [-1, 1, 1, -1],
        fill=True,
        color="Grey",
        alpha=0.1,
    )
    # t_line = ax.errorbar(mag_bins[:-1], bin_dT, yerr=bin_dT_err, color='blue', fmt='o')
    _ = ax.errorbar(mag_bins[:-1], bin_dT, yerr=bin_dT_err, label=label, fmt="o")
    # ax.axhline(y=0.003, linewidth=4, color='grey')
    # ax.legend([t_line], [r'$\delta T$'])
    ax.set_ylabel(r"$(T_{\rm PSF} - T_{\rm model}) \quad({\rm arcsec}^2)$")

    ax = axes[1]
    ax.set_ylim(-0.02, 0.02)
    ax.plot([min_mag, max_mag], [0, 0], color="black")
    ax.plot([min_mused, min_mused], [-1, 1], color="Grey")
    ax.fill(
        [min_mag, min_mag, min_mused, min_mused],
        [-1, 1, 1, -1],
        fill=True,
        color="Grey",
        alpha=0.1,
    )
    # t_line = ax.errorbar(mag_bins[:-1], bin_dTfrac, yerr=bin_dTfrac_err, color='blue', fmt='o')
    _ = ax.errorbar(
        mag_bins[:-1], bin_dTfrac, yerr=bin_dTfrac_err, label=label, fmt="o"
    )
    # ax.legend([t_line], [r'$\delta T$'])
    ax.legend(ncol=2, columnspacing=0.5)
    ax.set_ylabel(r"$(T_{\rm PSF} - T_{\rm model})/ T_{\rm PSF}$")

    ax = axes[2]
    ax.set_ylim(-6.0e-4, 6.0e-4)
    ax.plot([min_mag, max_mag], [0, 0], color="black")
    ax.plot([min_mused, min_mused], [-1, 1], color="Grey")
    ax.fill(
        [min_mag, min_mag, min_mused, min_mused],
        [-1, 1, 1, -1],
        fill=True,
        color="Grey",
        alpha=0.1,
    )
    # e1_line = ax.errorbar(mag_bins[:-1], bin_de1, yerr=bin_de1_err, color='red', mfc='white', fmt='o')
    # e2_line = ax.errorbar(mag_bins[:-1], bin_de2, yerr=bin_de2_err, color='blue', fmt='o')
    e1_line = ax.errorbar(
        mag_bins[:-1], bin_de1, yerr=bin_de1_err, mfc="white", fmt="o", label=label
    )
    e2_line = ax.errorbar(
        mag_bins[:-1], bin_de2, yerr=bin_de2_err, fmt="o", label=label
    )
    # ax.axhline(y=0.0002, linewidth=4, color='grey')
    # ax.axhline(y=-0.0002, linewidth=4, color='grey')
    # ax.legend(ncol=2, columnspacing=0.5)
    ax.legend([e1_line, e2_line], [r"$e_1$", r"$e_2$"], loc="lower left")
    ax.set_ylabel(r"$e_{\rm PSF} - e_{\rm model}$")
    ax.set_xlim(min_mag, max_mag)
    ax.set_xlabel(xlabel)
    plt.tight_layout()


class PlotProfileVMagByBands(Plotter):
    """Make profile plots by magnitude in a
    3 X nband grid"""

    default_config = dict(
        figsize=(24, 12),
        bands="griz",
    )

    def _make_plot(self, cat):

        bands = self.config["bands"]
        fig, axes = plt.subplots(
            3, 4, figsize=self.config["figsize"], sharey="row", sharex="col"
        )
        # tlims = [(0,0.6), (0.6,1.0), (1.0, 1.5), (1.5, 2.5)]
        for i, band in enumerate(bands):
            data = cat[cat["BAND"] == band]
            fracsizeres, sizeres, e1res, e2res = compute_res(data)
            zeropt = 30.0
            mag = zeropt - 2.5 * np.log10(data["FLUX"])
            print("Total in band: ", len(data))

            plot_bin_by_mag_on_axis(
                band, mag, sizeres, fracsizeres, e1res, e2res, axes=axes[:, i]
            )
        return fig


class PlotPsfByBand(Plotter):
    """Plot the PSF FWHM by band in a
    3 X nband grid"""

    default_config = dict(
        figsize=(24, 12),
        bands="griz",
    )

    def _make_plots(self, cat):
        bands = self.config["bands"]
        fig, axes = plt.subplots(
            3, 4, figsize=self.config["figsize"], sharey="row", sharex="col"
        )

        # tlims = [(0, 0.6), (0.6, 1.0), (1.0, 1.5), (1.5, 2.5)]
        fwhmlims = [(0.5, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 3.0)]
        for i, band in enumerate(bands):
            data = cat[cat["BAND"] == band]
            fracsizeres, sizeres, e1res, e2res = compute_res(data)
            zeropt = 30.0
            mag = zeropt - 2.5 * np.log10(data["FLUX"])
            print("Total in band: ", len(data))
            for lims in fwhmlims:
                t_cut = np.logical_and(
                    data["PSF_FWHM"] > lims[0], data["PSF_FWHM"] < lims[1]
                )
            mag_cut = mag[t_cut]
            fracsizeres_cut = fracsizeres[t_cut]
            sizeres_cut = sizeres[t_cut]
            e1res_cut = e1res[t_cut]
            e2res_cut = e2res[t_cut]

            print(
                band,
                "{} - {}asec: ".format(lims[0], lims[1]),
                len(mag_cut),
                "frac: ",
                len(mag_cut) / len(data),
            )
            plot_bin_by_mag_on_axis(
                band,
                mag_cut,
                sizeres_cut,
                fracsizeres_cut,
                e1res_cut,
                e2res_cut,
                axes=axes[:, i],
                label=str(lims),
            )
        return fig
