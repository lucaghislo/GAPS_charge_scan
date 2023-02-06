import os as os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

from plot_config import *
from error_function_calculator import compute_ERF_thrscan
from erf_function import *

# THRESHOLD SCAN
def threshold_scan(data_bkp, channels, n_events, output_folder):
    # Legend font size
    matplotlib.rcParams["legend.fontsize"] = 8.5

    excl_channels = np.setdiff1d(range(0, 32), channels)

    print("\nTHRESHOLD SCAN\n")
    print("Working on it, be patient...\n")
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
        (mu, sigma) = compute_ERF_thrscan(dac_range, events)
        plt.plot(
            dac_range,
            events,
            label=str(ch) + " THR: " + str(round(mu, 2)) + " DAC\_thr code",
            linestyle="--"
            if ch_count >= len(channels) / 2 and len(channels) > 16
            else "-",
        )

        ch_count = ch_count + 1

    if len(excl_channels) != 0:
        plt.title(
            r"\textbf{Threshold Scan ("
            + str(n_events)
            + " events, excl. ch. "
            + str(excl_channels)[1:-1].replace("'", "")
            + ")}"
        )
    else:
        plt.title(r"\textbf{Threshold Scan (" + str(n_events) + " events" + ")}")
    plt.ylim((-5, 105))
    plt.xlabel("Discriminator Threshold [DAC\_thr code]")
    plt.ylabel("Probability [\%]")
    num_columns = 1
    if len(channels) > 16:
        num_columns = 1
    plt.legend(
        title=r"\textbf{Channel}",
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        ncol=num_columns,
    )
    plt.grid()

    output_folder_spec = output_folder
    if not os.path.exists(output_folder_spec):
        os.mkdir(output_folder_spec)

    allch_filename = os.path.join(
        output_folder_spec,
        "threshold_scan_ch"
        + str(channels[0])
        + "-"
        + str(channels[len(channels) - 1])
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
        dac_range = ch_data.iloc[:, 0].to_numpy()
        events = ch_data.iloc[:, 3]
        events = [ev_i / n_events * 100 for ev_i in events]
        (mu, sigma) = compute_ERF_thrscan(dac_range, events)
        plt.plot(
            dac_range,
            events,
            label="THR: "
            + str(round(mu, 5))
            + " DAC\_thr code\n ENC: "
            + str(round(sigma, 5))
            + " DAC\_thr code",
        )
        if len(excl_channels) != 0:
            plt.title(
                r"\textbf{Threshold Scan ch. "
                + str(ch)
                + " ("
                + str(n_events)
                + " events, excl. ch. "
                + str(excl_channels)[1:-1].replace("'", "")
                + ")}"
            )
        else:
            plt.title(
                r"\textbf{Threshold Scan ch. "
                + str(ch)
                + " ("
                + str(n_events)
                + " events"
                + ")}"
            )
        plt.xlabel("Discriminator Threshold [DAC\_thr code]")
        plt.ylabel("Probability [\%]")
        plt.grid()
        plt.ylim((-5, 105))
        plt.legend(handlelength=0, handletextpad=0)
        plt.savefig(
            os.path.join(
                output_folder_spec_single,
                "charge_scan_ch" + str(ch) + ".pdf",
            )
        )

        parameters = np.vstack([parameters, [mu, sigma]])
        print("* Channel " + str(ch) + " *")
        print(
            "   mu: "
            + str(mu)
            + " DAC_inj code\nsigma: "
            + str(sigma)
            + " DAC_inj code\n"
        )

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
    if len(excl_channels) != 0:
        plt.title(
            r"\textbf{Thresholds from Threshold Scan (excl. ch. "
            + str(excl_channels)[1:-1].replace("'", "")
            + ")}",
        )
    else:
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
    if len(excl_channels) != 0:
        plt.title(
            r"\textbf{Thresholds from Threshold Scan (excl. ch. "
            + str(excl_channels)[1:-1].replace("'", "")
            + ")}",
        )
    else:
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
    if len(excl_channels) != 0:
        plt.title(
            r"\textbf{ENC from Threshold Scan (excl. ch. "
            + str(excl_channels)[1:-1].replace("'", "")
            + ")}",
        )
    else:
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

    print("\n")
