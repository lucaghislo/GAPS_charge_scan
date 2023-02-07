import numpy as np
import pandas as pd
import scipy
from os import *
from scipy.optimize import minimize
import matplotlib.pyplot as plt

raw_pedestal_tau4_filepath = r"C:\Users\ghisl\Google Drive UniBG\UniBG\CORSI\PhD\GAPS\SSL_Berkeley\L4R0M012345_automated_tests\MODULE_400\IT_L4R0M0_Gigi_m27.2C\data\Pedestals_tau4.dat"


def single_norm(x, *args):
    (m1, s1, k1) = args
    return k1 * scipy.stats.norm.pdf(x, loc=m1, scale=s1)


def get_pedestal_peak(raw_data_ch, pedestal_ch, ch):

    binwidth = 1

    raw_pedestals = pd.read_csv(
        raw_pedestal_tau4_filepath, sep="\t", comment="#", header=None
    )

    raw_pedestals_ch = raw_pedestals[
        (raw_pedestals.loc[:, 3] == ch)
        & ((raw_pedestals.loc[:, 2] == 0) | (raw_pedestals.loc[:, 2] == 1))
    ]
    print(raw_pedestals_ch)
    print(raw_pedestals_ch.loc[:, 4])

    raw_pedestals_ch = raw_pedestals_ch.to_numpy()

    # Extract mean pedestal from pedestal raw data
    (counts, bins, patches) = plt.hist(
        raw_pedestals_ch[
            (raw_pedestals_ch >= pedestal_ch * 0.9)
            & (raw_pedestals_ch <= pedestal_ch * 1.1)
        ],
        bins=np.arange(pedestal_ch * 0.9, pedestal_ch * 1.1, 1),
    )

    bins = bins[1::]
    x = bins
    y = counts

    print(x)
    print(y)

    params = [pedestal_ch, 1, 1]
    fitted_params, _ = scipy.optimize.curve_fit(
        single_norm, x, y, p0=params, maxfev=20000
    )

    fitted_pedestal = fitted_params

    # Extract profile from histogram
    (counts, bins, patches) = plt.hist(
        raw_data_ch[
            (raw_data_ch >= pedestal_ch * 0.9) & (raw_data_ch <= pedestal_ch * 1.1)
        ],
        bins=np.arange(pedestal_ch * 0.9, pedestal_ch * 1.1, 1),
    )

    bins = bins[1::]
    x = bins
    y = counts

    params = [pedestal_ch, 1, 1]
    fitted_params, _ = scipy.optimize.curve_fit(
        single_norm, x, y, p0=params, maxfev=20000
    )

    return fitted_params[0], fitted_params[1], fitted_pedestal[0], fitted_pedestal[1]
