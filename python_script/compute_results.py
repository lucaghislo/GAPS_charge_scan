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
max_ch = 31
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
    # 222,
    # 223,
    # 224,
    # 225,
]
input_folder = "python_script\output"

# Read data
# Fine thresholds
FTHR_thresholds = []
for thr in thr_list:
    filepath_root = os.path.join(
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
    filepath_subfolder = os.path.join(filepath_root, "ENC_THR")
    filepath = os.path.join(
        filepath_subfolder, "ch" + str(min_ch) + "-" + str(max_ch) + "_THR_ENC.dat"
    )
    data_raw = pd.read_csv(filepath, sep="\t", header=None)
    data_raw = data_raw.iloc[:, [1, 3]]

    # Channel threshold and ENC
    allch_thr = data_raw.iloc[:, 0]
    allch_enc = data_raw.iloc[:, 1]

    # print(thr)
    print(len(allch_thr[allch_thr - allch_enc > 0]))
    # print(max(allch_thr))
    # print("")

    thr_mean = np.mean(allch_thr[allch_thr > 0])
    thr_sigma = np.std(allch_thr[allch_thr > 0])

    FTHR_thresholds.append(thr_mean)


plt.clf()
plt.plot(thr_list, FTHR_thresholds)
plt.show()
