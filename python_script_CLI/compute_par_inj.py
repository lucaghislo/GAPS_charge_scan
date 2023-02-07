import os as os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit

from plot_config import *
from erf_function import *
from read_transfer_function import get_fdt, read_transfer_function
from read_pedestals import get_pedestal_auto
from calculate_xray_gain import get_linear_gain_realfdt


def linear_model(x, m, q):
    return m * x + q


def get_parasitic_injection(pedestal_filepath, fdt_filepath, ch, pt):
    # Read channel pedestal given peaking time
    pedestal_ch = get_pedestal_auto(pedestal_filepath, ch)

    # Calculate transfer function linear gain and estimate pedestal
    (fdt_gain, fdt_pedestal) = get_linear_gain_realfdt(fdt_filepath, ch, pt, 200)

    # Determine parasitic injection
    ch_par_inj = abs(pedestal_ch - fdt_pedestal)

    print(
        "Ch. "
        + str(ch)
        + "\tGain: "
        + str(np.round(fdt_gain, 5))
        + "\tPedestal fdt:"
        + str(np.round(fdt_pedestal, 5))
        + "\tPedestal auto: "
        + str(np.round(pedestal_ch, 5))
        + "\tPar. inj.: "
        + str(np.round(ch_par_inj, 5))
    )

    return ch_par_inj
