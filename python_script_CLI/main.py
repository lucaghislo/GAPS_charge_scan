import os as os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from plot_config import *
from erf_function import *
from charge_scan_noinj import charge_scan_noinj
from charge_scan import charge_scan
from threshold_scan import threshold_scan
from compute_par_inj import get_parasitic_injection


while True:
    # TODO check input correctness
    print("\n*** GAPS CHARGE SCAN TOOL v1.0 ***\n")

    # Request user input
    filename_chargescan = input("    Charge or threshold scan filepath: ")
    if filename_chargescan[0] == '"':
        filename_chargescan = filename_chargescan.replace('"', "")
    ch_min = int(input("                        First channel: "))
    ch_max = int(input("                         Last channel: "))
    excl_channels = input("  Excluded channels (comma separated): ")

    if excl_channels != "":
        excl_channels = excl_channels.split(",")
        excl_channels = [int(i) for i in excl_channels]
    else:
        excl_channels = []

    output_folder_filepath = input("               Output folder filepath: ")
    if output_folder_filepath[0] == '"':
        output_folder_filepath = output_folder_filepath.replace('"', "")

    # LaTex interpreter
    plt.rcParams.update({"text.usetex": True, "font.family": "serif"})

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
        filename_chargescan,
        comment="#",
        sep="\t",
        header=None,
    )
    data_bkp = data

    # Configuration
    conv_factor = 0.841
    channels = range(ch_min, ch_max + 1)
    channels = np.setdiff1d(channels, excl_channels)

    # Determine if charge scan or threshold scan
    n_events = data.iloc[0][2]
    threshold_col = data.iloc[:, 0]
    threshold_col = threshold_col.to_numpy()
    thr_unique = np.unique(threshold_col)
    charge_scan_flag = True

    if len(thr_unique) > 1:
        charge_scan_flag = False

    if charge_scan_flag:
        # Ask user to compensate parasitic injection or not
        comp_parinj_flag = input("Compensate parasitic injection? (y/n): ")

        # Charge scan with removal of parasitic injection
        if comp_parinj_flag == "y":
            # Get additional info when charge scan is selected and user wants to compensate parasitic injection
            filename_pedestal = input("         Pedestal from automated test: ")
            if filename_pedestal[0] == '"':
                filename_pedestal = filename_pedestal.replace('"', "")
            filename_fdt = input("Transfer function from automated test: ")
            if filename_fdt[0] == '"':
                filename_fdt = filename_fdt.replace('"', "")
            peaking_time = int(input("                Peaking time (0 to 7): "))

        # Charge scan without removal of parasitic injection
        # Always done
        (xmin, xmax) = charge_scan(data, channels, conv_factor, output_folder_filepath)

        # Charge scan with removal of parasitic injection
        if comp_parinj_flag == "y":
            # Get parasitic injection estimate to get proper charge scan
            allch_par_inj_estimate = []
            for ch in channels:
                ch_par_inj_estimate = get_parasitic_injection(
                    filename_pedestal, filename_fdt, ch, peaking_time
                )
                allch_par_inj_estimate.append(ch_par_inj_estimate)

            # Charge scan with subtracted parasitic injection
            charge_scan_noinj(
                data,
                channels,
                conv_factor,
                output_folder_filepath,
                xmin,
                xmax,
            )

    else:
        # Threshold scan
        threshold_scan(data_bkp, channels, n_events, output_folder_filepath)
