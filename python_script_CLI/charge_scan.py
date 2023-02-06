import os as os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

from plot_config import *
from error_function_calculator import compute_ERF
from erf_function import *

# CHARGE SCAN
def charge_scan(data, channels, conv_factor, output_folder):
    print("\nCHARGE SCAN\n")
    print("Working on it, be patient...\n")

    excl_channels = np.setdiff1d(range(0, 32), channels)

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

    if len(excl_channels) != 0:
        plt.title(
            r"\textbf{Charge scan (THR: "
            + str(threshold)
            + ", "
            + str(n_events)
            + " events, excl. ch. "
            + str(excl_channels)[1:-1].replace("'", "")
            + ")}"
        )
    else:
        plt.title(
            r"\textbf{Charge scan (THR: "
            + str(threshold)
            + ", "
            + str(n_events)
            + " events"
            + ")}"
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
    axes = plt.gca()
    xmin, xmax = axes.get_xlim()

    output_folder_spec = output_folder
    if not os.path.exists(output_folder_spec):
        os.mkdir(output_folder_spec)

    allch_filename = os.path.join(
        output_folder_spec,
        "charge_scan_ch"
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

    output_folder_spec_single_plot = os.path.join(output_folder_spec_single, "plots")
    if not os.path.exists(output_folder_spec_single_plot):
        os.mkdir(output_folder_spec_single_plot)

    output_folder_spec_single_data = os.path.join(output_folder_spec_single, "data")
    if not os.path.exists(output_folder_spec_single_data):
        os.mkdir(output_folder_spec_single_data)

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
            + " events"
            + ")}"
        )
        plt.xlabel("Energy [keV]")
        plt.ylabel("Probability [\%]")
        plt.grid()
        plt.ylim((-5, 105))
        plt.legend(handlelength=0, handletextpad=0)
        plt.savefig(
            os.path.join(
                output_folder_spec_single_plot,
                "charge_scan_ch" + str(ch) + "_THR_" + str(threshold) + ".pdf",
            )
        )

        parameters = np.vstack([parameters, [mu, sigma]])
        print("* Channel " + str(ch) + " *")
        print("   mu: " + str(mu) + " keV\nsigma: " + str(sigma) + " keV\n")

        # Save data for given channel
        with open(
            os.path.join(
                output_folder_spec_single_data,
                "ch_" + str(ch) + "_THR_" + str(threshold) + ".dat",
            ),
            "w",
        ) as filehandle:
            for i in range(0, len(inj_range)):
                filehandle.write("   %f\t%f\n" % (inj_range[i], events[i]))

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
            r"\textbf{Thresholds from Charge Scan (excl. ch. "
            + str(excl_channels)[1:-1].replace("'", "")
            + ")}",
        )
    else:
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
    if len(excl_channels) != 0:
        plt.title(
            r"\textbf{Thresholds from Charge Scan (excl. ch. "
            + str(excl_channels)[1:-1].replace("'", "")
            + ")}",
        )
    else:
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
    if len(excl_channels) != 0:
        plt.title(
            r"\textbf{ENC from Charge Scan (excl. ch. "
            + str(excl_channels)[1:-1].replace("'", "")
            + ")}",
        )
    else:
        plt.title(r"\textbf{ENC from Charge Scan}")

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

    return xmin, xmax
