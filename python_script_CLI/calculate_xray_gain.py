from scipy.optimize import curve_fit
from read_transfer_function import *
import numpy as np

conv_factor = 0.841


def linear_model(x, m, q):
    return m * x + q


# Linear interpolation with ADU as output
def get_linear_gain_realfdt(filepath, ch, pt, max_dacinj, outpath=""):

    # Get fdt data for given ch and pt
    cal_v, out = get_fdt(read_transfer_function(filepath), ch, pt)
    cal_v_kev = [cal_v_i * conv_factor for cal_v_i in cal_v]

    max_index = np.where(np.array(cal_v) >= max_dacinj)[0][0]

    out_selected = out[1:max_index]
    cal_v_kev_selected = cal_v_kev[1:max_index]

    popt, pcov = curve_fit(linear_model, cal_v_kev_selected, out_selected)

    gain = popt[0]
    pedestal = abs(popt[1])

    return gain, pedestal
