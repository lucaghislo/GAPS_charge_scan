import os as os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

from plot_config import *
from error_function_calculator import compute_ERF, compute_ERF_thrscan
from erf_function import *
from charge_scan_noinj import charge_scan_noinj
from charge_scan import charge_scan
from threshold_scan import threshold_scan

# Configuration
# filename = "IT_L4R0M0_Gigi_charge_scan_THR_205_FTH_MX.dat"
# ch_min = 0
# ch_max = 31

while True:
    # Request user input
    filename = input("Charge or threshold scan filename in \input folder: ")
    ch_min = int(input("First channel: "))
    ch_max = int(input(" Last channel: "))

    # LaTex interpreter
    plt.rcParams.update({"text.usetex": True, "font.family": "serif"})

    # Folders
    input_folder = "input"
    output_folder = "output"

    # PLOT CONFIGURATION
    # Label size
    matplotlib.rcParams["axes.labelsize"] = 13
    # Tick label size
    matplotlib.rcParams["xtick.labelsize"] = 13
    matplotlib.rcParams["ytick.labelsize"] = 13
    # Figure size
    matplotlib.rcParams["figure.figsize"] = 6.4 * 1.5, 4.8 * 1.5
    # Legend font size
    matplotlib.rcParams["legend.fontsize"] = 10

    # Read data from file
    data = pd.read_csv(
        os.path.join(input_folder, filename), comment="#", sep="\t", header=None
    )
    data_bkp = data

    # Configuration
    conv_factor = 0.841
    channels = range(ch_min, ch_max + 1)
    filename = filename.replace(".dat", "")

    # Determine if charge scan or threshold scan
    n_events = data.iloc[0][2]
    threshold_col = data.iloc[:, 0]
    threshold_col = threshold_col.to_numpy()
    thr_unique = np.unique(threshold_col)
    charge_scan_flag = True

    if len(thr_unique) > 1:
        charge_scan_flag = False

    if charge_scan_flag:
        # Charge scan
        charge_scan(data, channels, conv_factor, output_folder, filename)

        # Charge scan with subtracted parasitic injection
        charge_scan_noinj(data, channels, conv_factor, output_folder, filename)

    else:
        # Threshold scan
        threshold_scan(data_bkp, channels, n_events, output_folder, filename)
