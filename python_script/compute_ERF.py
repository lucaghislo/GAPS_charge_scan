from scipy import special
import matplotlib.pyplot as plt
import numpy as np
import os.path
import pandas as pd
from erf_function import *
import scipy as sp
from scipy.stats import norm

filename = "IT_L4R0M5_Gigi_m27C_charge_scan_THR_205_noFTH.dat"
filepath = os.path.join("output\SSL_Berkeley\data", filename)

# Open file in read mode
data = pd.read_csv(
    filepath,
    sep="\t",
)

x = data.iloc[:, 0]
data = data.iloc[:, 1:33]

parameters = np.zeros([1, 2])

for ch in range(0, 32):

    ch_data = data.iloc[:, ch]

    popt, pcov = sp.optimize.curve_fit(erf_function, x, ch_data / 100)

    plt.plot(x, ch_data / 100)
    plt.plot(x, erf_function(x, *popt))
    # plt.show()

    mu = popt[0]
    sigma = popt[1] * 2.35

    parameters = np.vstack([parameters, [mu, sigma]])

    print([mu, sigma])

    print("channel " + str(ch) + " -> mu: " + str(mu) + "\tsigma: " + str(sigma))

parameters = parameters[1:, :]

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
plt.plot(x, p, "k", linewidth=2)
plt.title(
    "Channel Threshold\n" + str(filename),
    fontweight="bold",
)
plt.xlabel("Threshold [keV]")
plt.ylabel("Count")
plt.savefig("output\SSL_Berkeley\erf_fit_results\\" + filename + "_thresholds.pdf")


# Plot ENC derived from charge scan
plt.clf()
plt.plot(range(0, 32), parameters[:, 1], marker="o")
plt.title(
    "ENC from Charge Scan\n" + str(filename),
    fontweight="bold",
)
plt.savefig("output\SSL_Berkeley\erf_fit_results\\" + str(filename) + "_ENC.pdf")
