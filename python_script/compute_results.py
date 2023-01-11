import os as os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

from plot_config import *
from error_function_calculator import compute_ERF, compute_ERF_thrscan
from erf_function import *

# Configuration
min_ch = 0
max_ch = 0  # 31
pt = 5

fthr_selector = ["noFTHR", "FTHR"]
thr_list = [
    205,
    206,
    207,
    208,
    209,
    210,
    211,
    212,
    213,
    214,
    215,
    216,
    217,
    218,
    219,
    220,
    221,
    222,
    223,
    224,
    225,
]
input_folder = "python_script\output"
channels = range(min_ch, max_ch + 1)

# Read data
# Fine thresholds
FTHR_thresholds = []
for thr in thr_list:
    # Get THR and ENC data
    filepath_root_thr = os.path.join(
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
    filepath_subfolder_thr = os.path.join(filepath_root_thr, "ENC_THR")
    filepath_thr = os.path.join(
        filepath_subfolder_thr, "ch" + str(min_ch) + "-" + str(max_ch) + "_THR_ENC.dat"
    )
    data_thr_enc_raw = pd.read_csv(filepath_thr, sep="\t", header=None)
    data_thr_enc_raw = data_thr_enc_raw.iloc[:, [1, 3]]

    # Channel threshold and ENC
    allch_thr = data_thr_enc_raw.iloc[:, 0]
    allch_enc = data_thr_enc_raw.iloc[:, 1]

    # keV -> DAC_inj code
    allch_thr = [i / 0.841 for i in allch_thr]
    allch_enc = [i / 0.841 for i in allch_enc]

    print(allch_thr)

    # print(thr)
    print(len(allch_thr[allch_thr - allch_enc > 0]))
    # print(max(allch_thr))
    # print("")

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

plt.clf()
plt.plot(thr_list, FTHR_thresholds)
plt.show()
