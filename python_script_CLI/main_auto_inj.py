import os as os
import numpy as np
import pandas as pd
from pathlib import Path as path
import glob
from os.path import join

from plot_config import *
from erf_function import *
from charge_scan_noinj import charge_scan_noinj
from charge_scan import charge_scan
from threshold_scan import threshold_scan
from compute_par_inj import get_parasitic_injection

# *** CONFIGURATION ***
# Layer
layer = 2

# Where charge scan raw data folders are located
root_filepath_base = (
    r"C:\Users\ghisl\Downloads\charge_scan_layers_computed\charge_scan_layer_2"
)

# Compensate parasitic injection?
comp_inj_flag = True
if comp_inj_flag:
    # Where pedestal and transfer function folders are located
    additional_data_filepath = r"C:\Users\ghisl\Google Drive UniBG\UniBG\CORSI\PhD\GAPS\SSL_Berkeley\layer2_module_tests"

# Do not edit
root_filepath = join(root_filepath_base, "*")
leaf_filepath = r"data\ChargeScan_fast.dat"
leaf_filepath_out = r"output"
all_folders = glob.glob(root_filepath)
all_folders = all_folders[1::]
alllayer_mask = join(root_filepath_base, "layer" + str(layer) + "_mask.txt")
if comp_inj_flag:
    alllayer_mask_inj = join(root_filepath_base, "layer" + str(layer) + "_mask_inj.txt")

for folder in all_folders:
    filename_chargescan = os.path.join(folder, leaf_filepath)
    output_folder_filepath = os.path.join(folder, leaf_filepath_out)

    # Change accordingly (by hand)
    row = int(folder[82:83])
    module = int(folder[84:85])

    # Read data from file
    data = pd.read_csv(
        filename_chargescan,
        comment="#",
        sep="\t",
        header=None,
    )

    # Configuration
    ch_min = 0
    ch_max = 31
    deactivate_thr = 30
    deactivate_enc = 10
    excl_channels = []

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

    # Charge scan without removal of parasitic injection
    (xmin, xmax) = charge_scan(
        data,
        channels,
        conv_factor,
        output_folder_filepath,
        deactivate_thr,
        deactivate_enc,
    )

    if comp_inj_flag:
        filename_test = join(
            join(
                join(
                    additional_data_filepath,
                    "MODULE" + str(layer) + str(row) + str(module) + "_fast",
                ),
                "*",
            ),
            "data",
        )

        leaf_filepath_data = glob.glob(filename_test)[0]

        # Pedestal
        filename_pedestal = join(leaf_filepath_data, "Pedestals_tau5.dat")

        # Transfer function
        filename_fdt = join(leaf_filepath_data, "TransferFunction_fast_tau5.dat")

        # Peaking time
        peaking_time = 5

        # Write estimated parameters to file
        summary_data_filepath = os.path.join(
            output_folder_filepath,
            "summary_inj_ch" + str(ch_min) + "-" + str(ch_max) + ".dat",
        )
        with open(summary_data_filepath, "w") as fp:
            fp.write("ch\ttf_gain\ttf_pedestal\tauto_pedestal\tpar_inj\n")
        fp.close()

        # Get parasitic injection estimate to get proper charge scan
        allch_par_inj_estimate = []
        for ch in channels:
            ch_par_inj_estimate = get_parasitic_injection(
                filename_pedestal,
                filename_fdt,
                ch,
                peaking_time,
                summary_data_filepath,
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
            allch_par_inj_estimate,
            deactivate_thr,
            deactivate_enc,
        )

    # Write activation mask to summary output file (no compensation)
    act_mask_path = join(output_folder_filepath, "ch0-31_activation_mask.txt")
    act_mask_file = open(act_mask_path, "r")
    act_mask = act_mask_file.read()
    with open(alllayer_mask, "a") as fp:
        fp.write(str(layer) + str(row) + str(module) + "\t" + str(act_mask) + "\n")
    fp.close()

    if comp_inj_flag:
        # Write activation mask to summary output file (compensation)
        act_mask_inj_path = join(
            output_folder_filepath, "ch0-31_activation_mask_inj.txt"
        )
        act_mask_file_inj = open(act_mask_inj_path, "r")
        act_mask_inj = act_mask_file_inj.read()
        with open(alllayer_mask_inj, "a") as fp:
            fp.write(
                str(layer) + str(row) + str(module) + "\t" + str(act_mask_inj) + "\n"
            )
        fp.close()
