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
    # Interpolate erf function
    popt, pcov = sp.optimize.curve_fit(erf_function, x, ch_data)

    # Get mean: threshold
    mu = popt[0]

    # Get standard deviation: dispersion
    sigma = popt[1] * 2.35

    return mu, sigma
