import os as os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

from plot_config import *
from error_function_calculator import compute_ERF
from erf_function import *

# Configuration
# filename = "IT_L4R0M0_Gigi_charge_scan_THR_205_FTH_MX.dat"
# ch_min = 0
# ch_max = 31

# Request user input
filename = input("Charge or threshold scan filename in \input folder: ")
ch_min = int(input("First channel: "))
ch_max = int(input(" Last channel: "))

print("\nWorking on it, be patient...\n")

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

if len(thr_unique > 1):
    charge_scan_flag = False
else:
    charge_scan_flag = True

if charge_scan_flag:
    # CHARGE SCAN
    # All channels in the same plot
    threshold = data.iloc[0][0]
    n_events = data.iloc[0][2]
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
            + " THR: "
            + str(round(mu, 2))
            + " keV\n ENC: "
            + str(round(sigma, 2))
            + " keV",
            linestyle="--"
            if ch_count >= len(channels) / 2 and len(channels) > 16
            else "-",
        )
        ch_count = ch_count + 1

    plt.title(
        r"\textbf{Charge scan (THR: "
        + str(threshold)
        + ", "
        + str(n_events)
        + " events)}"
    )
    plt.ylim((-5, 105))
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
        output_folder_spec, filename.replace(".dat", "") + ".pdf"
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
            label="THR: "
            + str(round(mu, 5))
            + " keV\n ENC: "
            + str(round(sigma, 5))
            + " keV",
        )
        plt.title(
            r"\textbf{Charge Scan ch. "
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
        plt.ylim((-5, 105))
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
        r"\textbf{Thresholds from Charge Scan}",
    )
    plt.xlabel("Threshold [keV]")
    plt.ylabel("Count")

    mu, std = norm.fit(data)

    matplotlib.pyplot.text(
        min(bins),
        max(n),
        "$\mu$ = "
        + str(round(mu, 5))
        + " keV\n $\sigma$ = "
        + str(round(std, 5))
        + " keV",
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
        r"\textbf{Thresholds from Charge Scan}",
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
            "ch"
            + str(channels[0])
            + "-"
            + str(channels[len(channels) - 1])
            + "_ENC.pdf",
        )
    )

else:
    # THRESHOLD SCAN
    # All channels in the same plot
    data = data_bkp
    plt.clf()
    ch_count = 0
    for ch in channels:
        min_thr = data.iloc[0][0]
        ch_data = data[data.iloc[:, 4] == ch]
        dac_range = ch_data.iloc[:, 0].to_numpy()
        events = ch_data.iloc[:, 3]
        events = [ev_i / n_events * 100 for ev_i in events]

        (mu, sigma) = compute_ERF(dac_range, events)

        print("DAC range")
        print(dac_range)

        print("\nEvents")
        print(events)

        plt.plot(
            dac_range,
            events,
            label=str(ch)
            + " THR: "
            + str(round(mu, 2))
            + " DAC\_thr code\n ENC: "
            + str(round(sigma, 2))
            + " DAC\_thr code",
            linestyle="--"
            if ch_count >= len(channels) / 2 and len(channels) > 16
            else "-",
        )
        ch_count = ch_count + 1

    plt.title(r"\textbf{Threshold Scan (" + str(n_events) + " events)}")
    plt.ylim((-5, 105))
    plt.xlabel("Discriminator Threshold [DAC\_thr code]")
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
        output_folder_spec, filename.replace(".dat", "") + ".pdf"
    )
    plt.savefig(allch_filename)
    print("Saved: " + allch_filename + "\n")

    # Print sigmoide
    plt.clf()
    plt.plot(dac_range, events, marker="o")
    plt.plot(dac_range, erf_function(dac_range, mu, sigma))
    plt.show()
    print(erf_function(dac_range, mu, sigma))

    print(len(dac_range))
    print(len(events))

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
        dac_range = ch_data.iloc[:, 0].to_numpy()
        events = ch_data.iloc[:, 3]
        events = [ev_i / n_events * 100 for ev_i in events]
        (mu, sigma) = compute_ERF(dac_range, events)
        plt.plot(
            dac_range,
            events,
            label="THR: "
            + str(round(mu, 5))
            + " DAC\_thr code\n ENC: "
            + str(round(sigma, 5))
            + " DAC\_thr code",
        )
        plt.title(
            r"\textbf{Threshold Scan ch. "
            + str(ch)
            + " ("
            + str(n_events)
            + " events)}"
        )
        plt.xlabel("Discriminator Threshold [DAC\_thr code]")
        plt.ylabel("Probability [\%]")
        plt.grid()
        plt.ylim((-5, 105))
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
        r"\textbf{Thresholds from Threshold Scan}",
    )
    plt.xlabel("Threshold [DAC\_thr code]")
    plt.ylabel("Count")

    mu, std = norm.fit(data)

    matplotlib.pyplot.text(
        min(bins),
        max(n),
        "$\mu$ = "
        + str(round(mu, 5))
        + " DAC\_thr code\n $\sigma$ = "
        + str(round(std, 5))
        + " DAC\_thr code",
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
    plt.ylabel("Threshold [DAC\_thr code]")
    plt.title(
        r"\textbf{Thresholds from Threshold Scan}",
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
    plt.ylabel("ENC [DAC\_thr code]")
    plt.title(
        r"\textbf{ENC from Threshold Scan}",
    )

    plt.grid()

    plt.savefig(
        os.path.join(
            ENC_THR_folder,
            "ch"
            + str(channels[0])
            + "-"
            + str(channels[len(channels) - 1])
            + "_ENC.pdf",
        )
    )
