import matplotlib.pyplot as plt
import numpy as np
import os.path
import scipy as sp
from scipy.stats import norm

from erf_function import *


# INPUT PARAMETERS
# x: DAC_inj converted in keV
# ch_data: channel trigger profile
def compute_ERF(x, ch_data):

    ch_data = [ch_i / 100 for ch_i in ch_data]

    # Interpolate erf function
    popt, pcov = sp.optimize.curve_fit(erf_function, x, ch_data)

    # Get mean: threshold
    mu = popt[0]

    # Get standard deviation: dispersion
    sigma = popt[1] * 2.35

    return mu, sigma


def compute_ERF_thrscan(dac_range, events):
    """
    It takes a list of DAC values and a list of event counts, and returns the mean and standard
    deviation of the DAC values where the event counts are non-zero

    :param dac_range: the DAC values that were scanned
    :param events: list of events per DAC value
    :return: The mean and standard deviation of the dac_range_lim array.
    """
    dac_range_lim = []
    events_lim = []
    max_counter = 0
    reached = True
    for i in range(0, len(dac_range)):

        if not reached and events[i] == 0:
            dac_range_lim = []
            events_lim = []

        if events[i] > 0 and max_counter < 1:
            if events[i] == 100:
                max_counter = max_counter + 1

            if reached:
                events_lim.append(0)
                dac_range_lim.append(dac_range[i - 1])

            reached = False
            dac_range_lim.append(dac_range[i])
            events_lim.append(events[i])

    mu = np.mean(dac_range_lim)
    sigma = np.std(dac_range_lim) * 2.35

    return mu, sigma
