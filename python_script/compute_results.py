import os as os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit

from plot_config import *
from error_function_calculator import compute_ERF, compute_ERF_thrscan
from erf_function import *


def linear_model(x, m, q):
    return m * x + q


# Configuration
min_ch = 0
max_ch = 31
pt = 5

thr_list = [
    205,
    206,
    207,
    208,
    # 209,
    # 210,
    # 211,
    # 212,
    # 213,
    # 214,
    # 215,
    # 216,
    # 217,
    # 218,
    # 219,
    # 220,
    # 221,
    # 222,
    # 223,
    # 224,
    # 225,
]
input_folder = "python_script\output"
channels = range(min_ch, max_ch + 1)

# Read data
# Fine thresholds
FTHR_thresholds = []
thr_index = 0
allch_thr_full = np.zeros(shape=(len(channels), len(thr_list)))
for thr in thr_list:
    # Get THR and ENC data
    filepath_root_thr = os.path.join(
        input_folder,
        "noFTHR_THR_"
        + str(thr)
        + "_pt"
        + str(pt)
        + "_ch_"
        + str(min_ch)
        + "-"
        + str(max_ch),
    )
    filepath_subfolder_thr = os.path.join(filepath_root_thr, "ENC_THR")
    filepath_thr = os.path.join(
        filepath_subfolder_thr, "ch" + str(min_ch) + "-" + str(max_ch) + "_THR_ENC.dat"
    )
    data_thr_enc_raw = pd.read_csv(filepath_thr, sep="\t", header=None)
    data_thr_enc_raw = data_thr_enc_raw.iloc[:, [1, 3]]

    # Channel threshold and ENC
    allch_thr = data_thr_enc_raw.iloc[:, 0]
    allch_enc = data_thr_enc_raw.iloc[:, 1]

    for i in range(0, len(allch_thr)):
        allch_thr_full[i, thr_index] = allch_thr[i]

    thr_index = thr_index + 1

    thr_mean = np.mean(allch_thr[allch_thr > 0])
    thr_sigma = np.std(allch_thr[allch_thr > 0])

    FTHR_thresholds.append(thr_mean)

    # Get channel dac_inj and events
    filepath_root_ch_0 = os.path.join(
        input_folder,
        "FTHR_THR_"
        + str(thr)
        + "_pt"
        + str(pt)
        + "_ch_"
        + str(min_ch)
        + "-"
        + str(max_ch),
    )
    filepath_root_ch = os.path.join(filepath_root_ch_0, "single_channels")
    filepath_subfolder_ch = os.path.join(filepath_root_ch, "data")

    plt.clf()
    for ch in channels:
        # Get ch data for threshold thr
        filepath_ch = os.path.join(
            filepath_subfolder_ch, "ch_" + str(ch) + "_THR_" + str(thr) + ".dat"
        )
        data_ch_raw = pd.read_csv(filepath_ch, sep="\t", header=None)
        dac_inj = data_ch_raw.iloc[:, 0]
        events = data_ch_raw.iloc[:, 1]
        plt.plot(dac_inj, events)
    # plt.show()

x = FTHR_thresholds[0 : len(FTHR_thresholds)]  # keV
y = thr_list[0 : len(thr_list)]  # DAC_thr code

plt.clf()
plt.plot(y, x)
plt.show()

# Acquire thr scan data for given channel
thr_scan_data = pd.read_csv(
    r"python_script\output\thrscan_150-255_FTHR_pt5_ch_0-31\ENC_THR\ch0-31_THR_ENC.dat",
    sep="\t",
    header=None,
)

thr_scan_values_allch = thr_scan_data.iloc[:, 1]
thr_scan_values_allch = thr_scan_values_allch.to_numpy()

paras_inj_allch = []

print("** Single channel analysis **")
for ch in channels:
    ch_data = allch_thr_full[ch, :]
    ch_data = [i / 0.841 for i in ch_data]
    popt, pcov = curve_fit(
        linear_model, ch_data[0 : len(ch_data)], thr_list[0 : len(thr_list)]
    )

    m = popt[0]
    q = popt[1]

    plt.plot(ch_data[0 : len(ch_data)], thr_list[0 : len(thr_list)])

    # Reconstruct linear model
    reconstructed_lin = []
    for i in range(0, len(ch_data)):
        x = ch_data[i]
        val = linear_model(x, m, q)
        reconstructed_lin.append(val)

    plt.plot(ch_data, reconstructed_lin)
    # plt.plot(thr_list, [0] * len(y))
    plt.plot([0] * len(range(200, 255)), range(200, 255))

    paras_inj_allch.append(abs(q - thr_scan_values_allch[ch]))

    print(str(ch) + ": " + str(m) + "\t" + str(q) + "\t" + str(paras_inj_allch[ch]))
    # plt.show()

plt.clf()
plt.plot(range(0, 32), paras_inj_allch, marker="o")
plt.xlabel("Channel")
plt.ylabel("Parasitic injection [DAC\_thr code]")
plt.title(r"\textbf{Estimated parasitic injection (no FTHR)}")
plt.savefig(
    r"python_script\output\parasitic_injection_estimation\parasitic_inj_noFTHR.pdf"
)
