import os.path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy


def single_norm(x, *args):
    (m1, s1, k1) = args
    return k1 * scipy.stats.norm.pdf(x, loc=m1, scale=s1)


def get_pedestal_auto(filepath_ped, ch):
    # Het raw pedestal data from automated test
    raw_pedestals = pd.read_csv(filepath_ped, sep="\t", comment="#", header=None)
    raw_pedestals_ch = raw_pedestals[raw_pedestals.loc[:, 3] == ch]
    raw_pedestals_ch = raw_pedestals_ch.loc[:, 4].to_numpy()
    mean_ped = np.mean(raw_pedestals_ch)

    # Extract mean pedestal from pedestal raw data
    (counts, bins, patches) = plt.hist(
        raw_pedestals_ch[
            (raw_pedestals_ch >= mean_ped * 0.9) & (raw_pedestals_ch <= mean_ped * 1.1)
        ],
        bins=np.arange(mean_ped * 0.9, mean_ped * 1.1, 1),
    )

    bins = bins[1::]
    x = bins
    y = counts

    params = [mean_ped, 1, 1]
    fitted_params, _ = scipy.optimize.curve_fit(
        single_norm, x, y, p0=params, maxfev=20000
    )

    ped_auto = fitted_params[0]

    return ped_auto
