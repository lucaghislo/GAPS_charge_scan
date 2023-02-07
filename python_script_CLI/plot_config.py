import matplotlib
import matplotlib.pyplot as plt

# TEX INTERPRETER
# Interpreter set to LaTex
plt.rcParams.update({"text.usetex": True, "font.family": "serif"})
# plt.rcParams.update({"text.usetex": False})
# plt.rc("text.latex", preamble=r"\usepackage{asmath}")

# FIGURE
# Font size
matplotlib.rcParams["font.size"] = 12
# Figure size
matplotlib.rcParams["figure.figsize"] = 6.4 * 1.5, 4.8 * 1.5

# BOXPLOT
# Boxplot line width
matplotlib.rcParams["axes.linewidth"] = 0.7

# TITLE
# Title font size
matplotlib.rcParams["axes.titlesize"] = 13

# TICKS
# Tick width
matplotlib.rcParams["xtick.major.width"] = 0.7
matplotlib.rcParams["ytick.major.width"] = 0.7
# Tick direction
matplotlib.rcParams["xtick.direction"] = "inout"
matplotlib.rcParams["ytick.direction"] = "inout"
# Tick label size
matplotlib.rcParams["xtick.labelsize"] = 13
matplotlib.rcParams["ytick.labelsize"] = 13

# LEGEND
# Legend font size
matplotlib.rcParams["legend.fontsize"] = 10
matplotlib.rcParams["legend.framealpha"] = 1

# AXIS
# Label size
matplotlib.rcParams["axes.labelsize"] = 13

# GRID
# Grid transparency
matplotlib.rcParams["grid.alpha"] = 0.3
# Grid linewidth
matplotlib.rcParams["grid.linewidth"] = 0.7

# SAVEFIG
# White margins
matplotlib.rcParams["savefig.bbox"] = "tight"
matplotlib.rcParams["savefig.pad_inches"] = 0.05
