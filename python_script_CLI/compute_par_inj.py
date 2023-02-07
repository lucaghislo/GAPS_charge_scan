import os as os
import numpy as np
import pandas as pd

from plot_config import *
from erf_function import *
from read_pedestals import get_pedestal_auto
from calculate_xray_gain import get_linear_gain_realfdt


def linear_model(x, m, q):
    return m * x + q


def get_parasitic_injection(
    pedestal_filepath, fdt_filepath, ch, pt, summary_data_filepath
):
    # Read channel pedestal given peaking time
    pedestal_ch = get_pedestal_auto(pedestal_filepath, ch)

    # Calculate transfer function linear gain and estimate pedestal
    (fdt_gain, fdt_pedestal) = get_linear_gain_realfdt(fdt_filepath, ch, pt, 200)

    # Determine parasitic injection
    ch_par_inj = abs(pedestal_ch - fdt_pedestal)

    with open(summary_data_filepath, "a") as fp:
        fp.write(
            str(ch)
            + "\t"
            + str(np.round(fdt_gain, 5))
            + "\t"
            + str(np.round(fdt_pedestal, 5))
            + "\t"
            + str(np.round(pedestal_ch, 5))
            + "\t"
            + str(np.round(ch_par_inj, 5))
            + "\n"
        )
    fp.close()

    print(
        "Ch. "
        + str(ch)
        + "\tGain: "
        + str(np.round(fdt_gain, 5))
        + "\tPedestal fdt: "
        + str(np.round(fdt_pedestal, 5))
        + "\tPedestal auto: "
        + str(np.round(pedestal_ch, 5))
        + "\tPar. inj.: "
        + str(np.round(ch_par_inj, 5))
    )

    return ch_par_inj
