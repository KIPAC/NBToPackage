"""Utility functions for piff_syst_plots"""

import numpy as np


def compute_res(d):
    """Compute the shape residuals"""
    de1 = d["G1_DATA"] - d["G1_MODEL"]
    de2 = d["G2_DATA"] - d["G2_MODEL"]
    dt = d["T_DATA"] - d["T_MODEL"]
    dtfrac = dt / d["T_DATA"]
    # print("mean de = ", np.mean(de1), np.mean(de2))
    # print("mean dt = ", np.mean(dt))

    return dtfrac, dt, de1, de2


def bin_res_by_color(m, dT, dTfrac, de1, de2, min_edge=0.0, max_edge=3.5):
    """Compute profile of residuals binned by color"""
    bins = np.linspace(min_edge, max_edge, 30)
    # print("col_bins = ", bins)

    index = np.digitize(m, bins)
    # print('len(index) = ',len(index))
    bin_de1 = [de1[index == i].mean() for i in range(1, len(bins))]
    # print('bin_de1 = ',bin_de1)
    bin_de2 = [de2[index == i].mean() for i in range(1, len(bins))]
    # print('bin_de2 = ',bin_de2)
    bin_dT = [dT[index == i].mean() for i in range(1, len(bins))]
    # print('bin_dT = ',bin_dT)
    bin_dTfrac = [dTfrac[index == i].mean() for i in range(1, len(bins))]
    # print('bin_dTfrac = ',bin_dTfrac)
    bin_de1_err = [
        np.sqrt(de1[index == i].var() / len(de1[index == i]))
        for i in range(1, len(bins))
    ]
    # print('bin_de1_err = ',bin_de1_err)
    bin_de2_err = [
        np.sqrt(de2[index == i].var() / len(de2[index == i]))
        for i in range(1, len(bins))
    ]
    # print('bin_de2_err = ',bin_de2_err)
    bin_dT_err = [
        np.sqrt(dT[index == i].var() / len(dT[index == i])) for i in range(1, len(bins))
    ]
    # print('bin_dT_err = ',bin_dT_err)
    bin_dTfrac_err = [
        np.sqrt(dTfrac[index == i].var() / len(dTfrac[index == i]))
        for i in range(1, len(bins))
    ]
    # print('bin_dTfrac_err = ',bin_dTfrac_err)

    # Fix up nans
    for i in range(1, len(bins)):
        if i not in index:
            bin_de1[i - 1] = 0.0
            bin_de2[i - 1] = 0.0
            bin_dT[i - 1] = 0.0
            bin_dTfrac[i - 1] = 0.0
            bin_de1_err[i - 1] = 0.0
            bin_de2_err[i - 1] = 0.0
            bin_dT_err[i - 1] = 0.0
            bin_dTfrac_err[i - 1] = 0.0

    return (
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
    )
