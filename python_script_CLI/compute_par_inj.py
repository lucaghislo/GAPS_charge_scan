import os as os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit

from plot_config import *
from erf_function import *
from read_transfer_function import get_fdt, read_transfer_function
from read_pedestals import get_pedestal, read_pedestals
from calculate_xray_gain import get_linear_gain_realfdt


def linear_model(x, m, q):
    return m * x + q


# TODO
# Compute parasitic injection from pedestal difference with transfer function interpolation


def get_parasitic_injection(pedestal_filepath, fdt_filepath, ch, pt):
    # Read channel pedestal given peaking time
    pedestal_ch = get_pedestal(read_pedestals(pedestal_filepath), ch, pt)

    # Read channel transfer function given peaking time
    (fdt_cal_v, fdt_ch_out) = get_fdt(read_transfer_function(fdt_filepath), ch, pt)

    # Calculate transfer function linear gain and estimate pedestal

    (fdt_gain, fdt_pedestal) = get_linear_gain_realfdt(fdt_filepath, ch, pt, 200)

    ch_par_inj = abs(pedestal_ch - fdt_pedestal)

    return ch_par_inj
