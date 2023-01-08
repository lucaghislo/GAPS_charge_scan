import os as os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

from plot_config import *
from error_function_calculator import compute_ERF
from erf_function import *

# Configuration
filename = "IT_L4R0M0_Gigi_charge_scan_THR_205_FTH_MX"
ch_min = 0
ch_max = 31

# TODO
# Request user input
# filename = input("Charge scan filepath: ")
# ch_min = int(input("First channel:"))
# ch_max = int(input(" Last channel:"))

# LaTex interpreter
plt.rcParams.update({"text.usetex": True, "font.family": "serif"})

# Folders
input_folder = "python_script\input"
output_folder = "python_script\output"

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

data = pd.read_csv(os.path.join(input_folder, filename + ".dat"), comment="#", sep="\t")
threshold = data.iloc[0][0]
n_events = data.iloc[0][2]
conv_factor = 0.841

channels = range(ch_min, ch_max + 1)
dac_range = np.unique(data.iloc[:, 1])

# All channels in the same plot
plt.clf()
ch_count = 0
for ch in channels:
    ch_data = data[data.iloc[:, 4] == ch]
    inj_range = ch_data.iloc[:, 1]
    events = ch_data.iloc[:, 3]
    inj_range = [inj_i * conv_factor for inj_i in inj_range]
    events = [ev_i / n_events * 100 for ev_i in events]
    (mu, sigma) = compute_ERF(inj_range, events)
    plt.plot(
        inj_range,
        events,
        label=str(ch)
        + " $\mu$ = "
        + str(round(mu, 2))
        + " keV\n $\sigma$ = "
        + str(round(sigma, 2))
        + " keV",
        linestyle="--" if ch_count >= len(channels) / 2 and len(channels) > 16 else "-",
    )
    ch_count = ch_count + 1

plt.title(
    r"\textbf{Charge scan (THR: " + str(threshold) + ", " + str(n_events) + " events)}"
)
plt.xlabel("Energy [keV]")
plt.ylabel("Probability [\%]")
num_columns = 1
if len(channels) > 16:
    num_columns = 2
plt.legend(
    title=r"\textbf{Channel}",
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    ncol=num_columns,
)
plt.grid()

output_folder_spec = os.path.join(output_folder, filename)
if not os.path.exists(output_folder_spec):
    os.mkdir(output_folder_spec)

allch_filename = os.path.join(
    output_folder_spec,
    "charge_scan_ch"
    + str(ch_min)
    + "-"
    + str(ch_max)
    + "_THR"
    + str(threshold)
    + ".pdf",
)
plt.savefig(allch_filename)
print("Saved: " + allch_filename + "\n")

parameters = np.zeros([1, 2])

# Save single channels
output_folder_spec_single = os.path.join(output_folder_spec, "single_channels")
if not os.path.exists(output_folder_spec_single):
    os.mkdir(output_folder_spec_single)

# Legend font size
matplotlib.rcParams["legend.fontsize"] = 13

for ch in channels:
    plt.clf()
    ch_data = data[data.iloc[:, 4] == ch]
    inj_range = ch_data.iloc[:, 1]
    events = ch_data.iloc[:, 3]
    inj_range = [inj_i * conv_factor for inj_i in inj_range]
    events = [ev_i / n_events * 100 for ev_i in events]
    (mu, sigma) = compute_ERF(inj_range, events)
    plt.plot(
        inj_range,
        events,
        label="$\mu$ = "
        + str(round(mu, 5))
        + " keV\n $\sigma$ = "
        + str(round(sigma, 5))
        + " keV",
    )
    plt.title(
        r"\textbf{Charge scan ch. "
        + str(ch)
        + " (THR: "
        + str(threshold)
        + ", "
        + str(n_events)
        + " events)}"
    )
    plt.xlabel("Energy [keV]")
    plt.ylabel("Probability [\%]")
    plt.grid()
    plt.legend(handlelength=0, handletextpad=0)
    plt.savefig(
        os.path.join(
            output_folder_spec_single,
            "charge_scan_ch" + str(ch) + "_THR_" + str(threshold) + ".pdf",
        )
    )

    parameters = np.vstack([parameters, [mu, sigma]])
    print("channel " + str(ch) + " -> mu: " + str(mu) + "\tsigma: " + str(sigma))
    print("Saved ch. " + str(ch))

parameters = parameters[1:, :]
ENC_THR_folder = os.path.join(output_folder_spec, "ENC_THR")

if not os.path.exists(ENC_THR_folder):
    os.mkdir(ENC_THR_folder)

# Write parameters to file
with open(
    os.path.join(
        ENC_THR_folder,
        "ch"
        + str(channels[0])
        + "-"
        + str(channels[len(channels) - 1])
        + "_THR_ENC.dat",
    ),
    "w",
) as filehandle:
    for i in range(0, len(channels)):
        print("Parameter: " + str(parameters[i, 0]))
        filehandle.write(
            "%d\t%f\t\t%f\n" % (channels[i], parameters[i, 0], parameters[i, 1])
        )

# Plot histogram of threshold data
plt.clf()
data = parameters[:, 0]
plot_data = [int(data_i) for data_i in data]
(n, bins, hist) = plt.hist(
    data,
)
plt.title(
    r"\textbf{Thresholds from charge scan}",
)
plt.xlabel("Threshold [keV]")
plt.ylabel("Count")

mu, std = norm.fit(data)

matplotlib.pyplot.text(
    min(bins),
    max(n),
    "$\mu$ = " + str(round(mu, 5)) + " keV\n $\sigma$ = " + str(round(std, 5)) + " keV",
    fontsize=13,
    verticalalignment="top",
    bbox=dict(facecolor="white", edgecolor="#cdcdcd", boxstyle="round,pad=0.35"),
)

plt.savefig(
    os.path.join(
        ENC_THR_folder,
        "ch"
        + str(channels[0])
        + "-"
        + str(channels[len(channels) - 1])
        + "_THR_hist.pdf",
    )
)

# Plot threshold derived from charge scan
plt.clf()
plt.plot(range(0, len(channels)), parameters[:, 0], marker="o")
plt.xlabel("Channel")
plt.ylabel("Threshold [keV]")
plt.title(
    r"\textbf{Thresholds from charge scan}",
)

plt.grid()

plt.savefig(
    os.path.join(
        ENC_THR_folder,
        "ch"
        + str(channels[0])
        + "-"
        + str(channels[len(channels) - 1])
        + "_THR_plot.pdf",
    )
)

# Plot ENC derived from charge scan
plt.clf()
plt.plot(range(0, len(channels)), parameters[:, 1], marker="o")
plt.xlabel("Channel")
plt.ylabel("ENC [keV]")
plt.title(
    r"\textbf{ENC from Charge Scan}",
)

plt.grid()

plt.savefig(
    os.path.join(
        ENC_THR_folder,
        "ch" + str(channels[0]) + "-" + str(channels[len(channels) - 1]) + "_ENC.pdf",
    )
)
