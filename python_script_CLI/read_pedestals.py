import os.path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy


def single_norm(x, *args):
    (m1, s1, k1) = args
    return k1 * scipy.stats.norm.pdf(x, loc=m1, scale=s1)


def get_pedestal_auto(filepath_ped, ch):
    # Get raw pedestal data from automated test
    raw_pedestals = pd.read_csv(filepath_ped, sep="\t", comment="#", header=None)
    raw_pedestals_ch = raw_pedestals[raw_pedestals.loc[:, 3] == ch]
    raw_pedestals_ch = raw_pedestals_ch.loc[:, 4].to_numpy()

    mean_ped = np.mean(raw_pedestals_ch)
    ped_auto = mean_ped

    return ped_auto
