import os as os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

from plot_config import *
from error_function_calculator import compute_ERF
from erf_function import *

# CHARGE SCAN WITHOUT PARASITIC INJECTION ESTIMATED FROM PEDESTAL AND FDT INTERPOLATION
def charge_scan_noinj(
    data,
    channels,
    conv_factor,
    output_folder,
    xmin,
    xmax,
    par_inj_pedestal,
    deactivate_thr,
):

    # Legend font size
    matplotlib.rcParams["legend.fontsize"] = 10

    thr_limit_NaN = 300

    print("\nCHARGE SCAN")
    print("Charge scan without estimated parasitic injection\n")
    print("Working on it, be patient...\n")

    excl_channels = np.setdiff1d(range(0, 32), channels)

    # All channels in the same plot
    threshold = data.iloc[0][0]
    n_events = data.iloc[0][2]
    plt.clf()
    ch_count = 0
    lim_flags_nan = []
    lim_flags_thr = []
    for ch in channels:
        ch_data = data[data.iloc[:, 4] == ch]
        inj_range = ch_data.iloc[:, 1]
        events = ch_data.iloc[:, 3]
        inj_range = [inj_i * conv_factor - par_inj_pedestal[ch] for inj_i in inj_range]
        events = [ev_i / n_events * 100 for ev_i in events]
        (mu, sigma) = compute_ERF(inj_range, events)
        lim_flag_nan = False
        if mu > thr_limit_NaN:
            lim_flag_nan = True
        lim_flags_thr.append(mu < deactivate_thr or lim_flag_nan)
        lim_flags_nan.append(lim_flag_nan)
        plt.plot(
            inj_range,
            events,
            label=str(ch)
            + " THR: "
            + (str(round(mu, 2)) + " keV" if not lim_flag_nan else "nan")
            + " \n ENC: "
            + (str(round(sigma, 2)) + " keV" if not lim_flag_nan else "nan"),
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
            + r" events, parasitic injection removed, excl. ch. "
            + str(excl_channels)[1:-1].replace("'", "")
            + ")}"
        )
    else:
        plt.title(
            r"\textbf{Charge scan (THR: "
            + str(threshold)
            + ", "
            + str(n_events)
            + r" events, parasitic injection removed"
            + ")}"
        )
    plt.ylim((-5, 105))
    plt.xlabel("Energy [keV]")
    plt.ylabel("Probability [\%]")
    plt.xlim(xmin, xmax)
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

    output_folder_spec = output_folder
    if not os.path.exists(output_folder_spec):
        os.mkdir(output_folder_spec)

    allch_filename = os.path.join(
        output_folder_spec,
        "charge_scan_ch"
        + str(channels[0])
        + "-"
        + str(channels[len(channels) - 1])
        + "_inj.pdf",
    )
    plt.savefig(allch_filename)
    plt.savefig(allch_filename.replace(".pdf", ".png"))
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

    ch_counter = 0
    for ch in channels:
        plt.clf()
        ch_data = data[data.iloc[:, 4] == ch]
        inj_range = ch_data.iloc[:, 1]
        events = ch_data.iloc[:, 3]

        inj_range = [inj_i * conv_factor - par_inj_pedestal[ch] for inj_i in inj_range]
        events = [ev_i / n_events * 100 for ev_i in events]
        (mu, sigma) = compute_ERF(inj_range, events)
        plt.plot(
            inj_range,
            events,
            label="THR: "
            + (str(round(mu, 5)) + " keV" if not lim_flags_nan[ch_counter] else "nan")
            + "\n ENC: "
            + (
                str(round(sigma, 5)) + " keV"
                if not lim_flags_nan[ch_counter]
                else "nan"
            ),
        )
        if len(excl_channels) != 0:
            plt.title(
                r"\textbf{Charge Scan ch. "
                + str(ch)
                + " (THR: "
                + str(threshold)
                + ", "
                + str(n_events)
                + r" events, parasitic injection removed, excl. ch. "
                + str(excl_channels)[1:-1].replace("'", "")
                + ")}"
            )
        else:
            plt.title(
                r"\textbf{Charge Scan ch. "
                + str(ch)
                + " (THR: "
                + str(threshold)
                + ", "
                + str(n_events)
                + r" events, parasitic injection removed"
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
                "charge_scan_ch" + str(ch) + "_THR_" + str(threshold) + "_inj.pdf",
            )
        )
        plt.savefig(
            os.path.join(
                output_folder_spec_single_plot,
                "charge_scan_ch" + str(ch) + "_THR_" + str(threshold) + "_inj.png",
            )
        )

        parameters = np.vstack([parameters, [mu, sigma]])
        print("* Channel " + str(ch) + " *")
        print("   mu: " + str(mu) + " keV\nsigma: " + str(sigma) + " keV\n")

        # Save data for given channel
        with open(
            os.path.join(
                output_folder_spec_single_data,
                "ch_" + str(ch) + "_THR_" + str(threshold) + "_inj.dat",
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
            + "_THR_ENC_inj.dat",
        ),
        "w",
    ) as filehandle:
        filehandle.write("ch\tthr\tenc\n")
        for i in range(0, len(channels)):
            if not lim_flags_nan[i]:
                filehandle.write(
                    "%d\t%f\t\t%f\n" % (channels[i], parameters[i, 0], parameters[i, 1])
                )
            else:
                filehandle.write("%d\tnan\t\tnan\n" % (channels[i]))

    # Write activation mask to file
    with open(
        os.path.join(
            output_folder,
            "ch"
            + str(np.min(channels))
            + "-"
            + str(np.max(channels))
            + "_activation_mask_inj.txt",
        ),
        "w",
    ) as filehandle:
        allch_flag = ""
        for i in range(0, len(channels)):
            ch_flag = 0
            if not lim_flags_thr[i]:
                ch_flag = 1

            allch_flag = allch_flag + str(ch_flag)

        filehandle.write(
            "# Ch. " + str(np.min(channels)) + " to " + str(np.max(channels)) + "\n"
        )
        filehandle.write(allch_flag)
        filehandle.write(
            "\n\n# Ch. " + str(np.max(channels)) + " to " + str(np.min(channels)) + "\n"
        )
        filehandle.write(allch_flag[len(allch_flag) :: -1])

    filehandle.close()

    # Plot histogram of threshold data
    plt.clf()
    for i in range(0, len(parameters[:, 0])):
        if lim_flags_nan[i]:
            parameters[i, 0] = "nan"
            parameters[i, 1] = "nan"

    data = []
    plot_data = []
    for i in range(0, len(parameters[:, 0])):
        if not lim_flags_nan[i]:
            plot_data.append(int(parameters[i, 0]))
            data.append(parameters[i, 0])

    (n, bins, hist) = plt.hist(data)

    if len(excl_channels) != 0:
        plt.title(
            r"\textbf{Thresholds from Charge Scan (parasitic injection removed, excl. ch. "
            + str(excl_channels)[1:-1].replace("'", "")
            + ")}",
        )
    else:
        plt.title(
            r"\textbf{Thresholds from Charge Scan (parasitic injection removed)}",
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
            + "_THR_hist_inj.pdf",
        )
    )

    plt.savefig(
        os.path.join(
            ENC_THR_folder,
            "ch"
            + str(channels[0])
            + "-"
            + str(channels[len(channels) - 1])
            + "_THR_hist_inj.png",
        )
    )

    # Plot threshold derived from charge scan
    plt.clf()
    plt.plot(range(0, len(channels)), parameters[:, 0], marker="o")
    plt.xlabel("Channel")
    plt.ylabel("Threshold [keV]")
    if len(excl_channels) != 0:
        plt.title(
            r"\textbf{Thresholds from Charge Scan (parasitic injection removed, excl. ch. "
            + str(excl_channels)[1:-1].replace("'", "")
            + ")}",
        )
    else:
        plt.title(r"\textbf{Thresholds from Charge Scan (parasitic injection removed)")

    plt.grid()

    plt.savefig(
        os.path.join(
            ENC_THR_folder,
            "ch"
            + str(channels[0])
            + "-"
            + str(channels[len(channels) - 1])
            + "_THR_plot_inj.pdf",
        )
    )

    plt.savefig(
        os.path.join(
            ENC_THR_folder,
            "ch"
            + str(channels[0])
            + "-"
            + str(channels[len(channels) - 1])
            + "_THR_plot_inj.png",
        )
    )

    # Plot ENC derived from charge scan
    plt.clf()
    plt.plot(range(0, len(channels)), parameters[:, 1], marker="o")
    plt.xlabel("Channel")
    plt.ylabel("ENC [keV]")
    if len(excl_channels) != 0:
        plt.title(
            r"\textbf{ENC from Charge Scan (parasitic injection removed, excl. ch. "
            + str(excl_channels)[1:-1].replace("'", "")
            + ")}",
        )
    else:
        plt.title(r"\textbf{ENC from Charge Scan (parasitic injection removed)")

    plt.grid()

    plt.savefig(
        os.path.join(
            ENC_THR_folder,
            "ch"
            + str(channels[0])
            + "-"
            + str(channels[len(channels) - 1])
            + "_ENC_inj.pdf",
        )
    )

    plt.savefig(
        os.path.join(
            ENC_THR_folder,
            "ch"
            + str(channels[0])
            + "-"
            + str(channels[len(channels) - 1])
            + "_ENC_inj.png",
        )
    )
