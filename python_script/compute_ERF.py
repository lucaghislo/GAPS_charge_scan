from scipy import special
import matplotlib.pyplot as plt
import numpy as np
import os.path
import pandas as pd
from erf_function import *
import scipy as sp
from scipy.stats import norm

# Layer 4 Row 0 Module without fine trimming
filename_0 = "IT_L4R0M0_Gigi_m27.2C_charge_scan_THR_205_noFTH"
filename_1 = "IT_L4R0M1_Gigi_m30.4C_charge_scan_THR_205_noFTH"
filename_2 = "IT_L4R0M2_Gigi_m29.6C_charge_scan_THR_205_noFTH"
filename_3 = "IT_L4R0M3_Gigi_m29.6C_charge_scan_THR_205_noFTH"
filename_4 = "IT_L4R0M4_Gigi_m29.7C_charge_scan_THR_205_noFTH"
filename_5 = "IT_L4R0M5_Gigi_m27C_charge_scan_THR_205_noFTH"

# Layer 4 Row 0 Module with fine trimming
filename_0 = "IT_L4R0M0_Gigi_charge_scan_THR_205_FTH_MX"
filename_1 = "IT_L4R0M1_Gigi_charge_scan_THR_205_FTH_MX"  # Wrong FTH values
filename_2 = "IT_L4R0M2_Gigi_charge_scan_THR_205_FTH_MX"
filename_3 = "IT_L4R0M3_Gigi_charge_scan_THR_205_FTH_MX"
filename_4 = "IT_L4R0M4_Gigi_charge_scan_THR_205_FTH_MX"
filename_5 = "IT_L4R0M5_Gigi_charge_scan_THR_205_FTH_MX"

# Layer 4 Row 0 Module 1 with fine trimming and different threshold values
filename_200 = "IT_L4R0M1_Gigi_charge_scan_THR_200_FTH_LG"
filename_201 = "IT_L4R0M1_Gigi_charge_scan_THR_201_FTH_LG"
filename_203 = "IT_L4R0M1_Gigi_charge_scan_THR_203_FTH_LG"
filename_205 = "IT_L4R0M1_Gigi_charge_scan_THR_205_FTH_LG"
filename_207 = "IT_L4R0M1_Gigi_charge_scan_THR_207_FTH_LG"

filename = filename_203
path_in = filename + ".dat"
filepath = os.path.join("output\SSL_Berkeley\FTH\L4R0M1\data", path_in)  # \FTH\

# Open file in read mode
data = pd.read_csv(
    filepath,
    sep="\t",
)

# Get data length
x = data.iloc[:, 0]
data = data.iloc[:, 1:33]

parameters = np.zeros([1, 2])

# Analyze all channels
for ch in range(0, 32):

    # Get channel data
    ch_data = data.iloc[:, ch]

    # Interpolate erf function
    popt, pcov = sp.optimize.curve_fit(erf_function, x, ch_data / 100)

    # Plot data interpolation
    plt.plot(x, ch_data / 100)
    plt.plot(x, erf_function(x, *popt))

    # Get mean: threshold
    mu = popt[0]

    # Get standard deviation: dispersion
    sigma = popt[1] * 2.35

    if filename == filename_2:
        if ch == 7:
            mu = 0
            sigma = 0

    parameters = np.vstack([parameters, [mu, sigma]])
    print("channel " + str(ch) + " -> mu: " + str(mu) + "\tsigma: " + str(sigma))

parameters = parameters[1:, :]

# Write parameters to file
with open(
    os.path.join(
        "output\SSL_Berkeley\erf_fit_results\FTH\L4R0M1\data",
        filename + "_THR_ENC.dat",  # \FTH\
    ),
    "w",
) as filehandle:
    for i in range(0, 32):
        print("Parameter: " + str(parameters[i, 0]))
        filehandle.write("%f %f\n" % (parameters[i, 0], parameters[i, 1]))

# Plot histogram of threshold data
plt.clf()
data = parameters[:, 0]
binwidth = 10
plot_data = [int(data_i) for data_i in data]
plt.hist(
    data,
    bins=range(min(plot_data), max(plot_data) + binwidth, binwidth),
    edgecolor="black",
)
# Plot the PDF
xmin, xmax = plt.xlim()
mu, std = norm.fit(data)
x = np.linspace(xmin - 15, xmax + 15, 100)
p = norm.pdf(x, mu, std) * 320
# plt.plot(x, p, "k", linewidth=2)
plt.title(
    "Channel Threshold\n" + str(filename),
    fontweight="bold",
)
plt.xlabel("Threshold [keV]")
plt.ylabel("Count")
plt.savefig(
    "output\SSL_Berkeley\erf_fit_results\FTH\L4R0M1\\" + filename + "_thresholds.pdf"
)  # \FTH\
plt.savefig(
    "output\SSL_Berkeley\erf_fit_results\FTH\L4R0M1\\" + filename + "_thresholds.png"
)  # \FTH\


# Plot ENC derived from charge scan
plt.clf()
plt.plot(range(0, 32), parameters[:, 1], marker="o")
plt.xlabel("Channel")
plt.ylabel("ENC [keV]")
plt.title(
    "ENC from Charge Scan\n" + str(filename),
    fontweight="bold",
)
plt.savefig(
    "output\SSL_Berkeley\erf_fit_results\FTH\L4R0M1\\" + filename + "_ENC.pdf"
)  # \FTH\
plt.savefig(
    "output\SSL_Berkeley\erf_fit_results\FTH\L4R0M1\\" + filename + "_ENC.png"
)  # \FTH\
