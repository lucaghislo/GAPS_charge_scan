import os as os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from plot_config import *

# Configuration
filename = "python_script\sample_data.dat"
ch_min = 0
ch_max = 31
conv_factor = 0.841

# LaTex interpreter
plt.rcParams.update({"text.usetex": True, "font.family": "serif"})

data = pd.read_csv(filename, comment="#", sep="\t")
threshold = data.iloc[0][0]
n_events = data.iloc[0][2]

channels = range(ch_min, ch_max + 1)
dac_range = np.unique(data.iloc[:, 1])

plt.clf()
for ch in channels:
    ch_data = data[data.iloc[:, 4] == ch]
    inj_range = ch_data.iloc[:, 1]
    events = ch_data.iloc[:, 3]
    inj_range = [inj_i * conv_factor for inj_i in inj_range]
    events = [ev_i / n_events * 100 for ev_i in events]
    plt.plot(inj_range, events)

plt.title(
    r"\textbf{Charge scan (THR: " + str(threshold) + ", " + str(n_events) + " events)}"
)
plt.xlabel("Energy [keV]")
plt.ylabel("Probability [\%]")
plt.show()
