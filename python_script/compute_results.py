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
# plt.show()

# Acquire thr scan data for given channel
thr_scan_data = pd.read_csv(
    r"python_script\output\thrscan_150-255_FTHR_pt5_ch_0-31\ENC_THR\ch0-31_THR_ENC.dat",
    sep="\t",
    header=None,
)

thr_scan_values_allch = thr_scan_data.iloc[:, 1]
thr_scan_values_allch = thr_scan_values_allch.to_numpy()

# Acquire pedestal data from automated test
pedestal_auto_raw = pd.read_csv(
    r"python_script\output\parasitic_injection_estimation\allch_pedestal_data.dat",
    sep="\t",
    header=None,
)

pedestal_auto = pedestal_auto_raw.iloc[:, 1]
pedestal_auto = pedestal_auto.to_numpy()

# Acquire estimated pedestal and gain from transfer function interpolation
pedestal_gain_fdt_raw = pd.read_csv(
    r"python_script\output\parasitic_injection_estimation\allchs_pt5_low_energy_gain_200_realfdt.dat",
    sep="\t",
    header=None,
)

gain_fdt = pedestal_gain_fdt_raw.iloc[:, 1]
gain_fdt = gain_fdt.to_numpy()
pedestal_fdt = pedestal_gain_fdt_raw.iloc[:, 2]
pedestal_fdt = pedestal_fdt.to_numpy()

estimated_coeff_filepath = (
    r"python_script\output\parasitic_injection_estimation\paras_inj_coeff_FTHR.dat"
)

est_coeff_handle = open(estimated_coeff_filepath, "w")
est_coeff_handle.write("ch\tmpar\tqpar\tqthr\tped\tpedinj\tchgain\n")
est_coeff_handle.close()
est_coeff_handle = open(estimated_coeff_filepath, "a")

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

    # print(
    #     str(ch)
    #     + ": "
    #     + str(m)
    #     + "\t"
    #     + str(q)
    #     + "\t"
    #     + str(paras_inj_allch[ch])
    #     + "\t"
    #     + str(pedestal_auto[ch])
    #     + "\t"
    #     + str()
    # )
    # plt.show()

    # Save data for given channel
    est_coeff_handle.write(
        str(ch)
        + "\t"
        + str(m)
        + "\t"
        + str(q)
        + "\t"
        + str(thr_scan_values_allch[ch])
        + "\t"
        + str(pedestal_auto[ch])
        + "\t"
        + str(pedestal_fdt[ch])
        + "\t"
        + str(gain_fdt[ch])
        + "\n"
    )

est_coeff_handle.close()

# Read coefficient file and compute parasitic injection comparison
data = pd.read_csv(
    r"python_script\output\parasitic_injection_estimation\paras_inj_coeff_FTHR.dat",
    sep="\t",
)

print(data)

channels = range(0, 32)
mpar = data["mpar"]
qpar = data["qpar"]
qthr = data["qthr"]
ped = data["ped"]
pedinj = data["pedinj"]
chgain = data["chgain"]

plt.clf()
plt.plot(range(0, 32), paras_inj_allch, marker="o")
plt.plot(range(0, 32), mpar * 40 + 8, marker="*")
plt.xlabel("Channel")
plt.ylabel("Parasitic injection [DAC\_thr code]")
plt.title(r"\textbf{Estimated parasitic injection (FTHR)}")

# Save plot to file
plt.savefig(
    r"python_script\output\parasitic_injection_estimation\parasitic_inj_FTHR.pdf"
)


# Parasitic injection estimate from charge + thr scan
par_inj_dacthr = [qthr[i] - qpar[i] for i in channels]
print(par_inj_dacthr)
print("")

parasitic_inj_pedestal = []
parasitic_inj_method = []

# Convert parasitic injection in DAC_inj
for ch in channels:
    # Iniezione parassita da differenza di piedistallo
    inj_pedestal = abs(ped[ch] - pedinj[ch]) * 0.841

    # Iniezione parassita
    par_inj_dacthr_ch = par_inj_dacthr[ch]  # [DAC_thr]

    # Thr da threshold scan
    qthr_ch = qthr[ch]  # [DAC_thr]

    # Parametri modello DAC_thr su Y e DAC_inj su X
    mpar_ch = mpar[ch]  # [DAC_thr/DAC_inj]
    qpar_ch = qpar[ch]  # [DAC_thr]

    # Parametri modello DAC_inj su Y e DAC_thr su X
    mpar1_ch = 1 / mpar_ch  # [DAC_inj/DAC_thr]
    qpar1_ch = qpar_ch / abs(mpar_ch)  # [DAC_inj]

    # Ordinata all'origine ottenuta da threshold scan
    qpar2_ch = qthr_ch / abs(mpar_ch)

    # print(str(mpar1_ch) + "\t" + str(qpar1_ch))

    inj_dacinj_temp = linear_model(par_inj_dacthr_ch, mpar1_ch, qpar1_ch)

    # Iniezione parassita
    inj_dacinj = abs(inj_dacinj_temp - qpar2_ch)  # [DAC_inj]

    # Inieizione parassita
    inj_kev = inj_dacinj * 0.841  # [keV]

    # Iniezione parassita OK
    inj_kev_ok = par_inj_dacthr_ch * abs(mpar1_ch) * 0.841  # [keV]

    parasitic_inj_pedestal.append(inj_pedestal)
    parasitic_inj_method.append(inj_kev_ok)

    print(
        str(ch)
        # + "\t"
        # + str(inj_kev)
        + "\t"
        + str(inj_pedestal)
        + "\t"
        + str(inj_kev_ok)
        + "\t"
        + str(inj_pedestal - inj_kev_ok)
    )

    # print(str(ch) + "\t" + str(par_inj_dacthr_ch * abs(mpar1_ch)))

plt.clf()
plt.plot(
    range(0, 32),
    parasitic_inj_pedestal,
    marker="o",
    label=r"From pedestal $\mu="
    + str(round(np.mean(parasitic_inj_pedestal), 2))
    + r"$ keV, $\sigma = "
    + str(round(np.std(parasitic_inj_pedestal), 2))
    + r"$ keV",
)
plt.plot(
    range(0, 32),
    parasitic_inj_method,
    marker="*",
    label=r"From charge/thr scan $\mu="
    + str(round(np.mean(parasitic_inj_method), 2))
    + r"$ keV, $\sigma = "
    + str(round(np.std(parasitic_inj_method), 2))
    + r"$ keV",
)
plt.title(r"\textbf{Parasitic injection estimate}")
plt.ylabel("Parasitic injection [keV]")
plt.xlabel("Channel")
plt.legend(loc="lower left")
plt.grid()
plt.savefig(
    r"python_script\output\parasitic_injection_estimation\parasitic_inj_FTHR_comparison.pdf"
)
