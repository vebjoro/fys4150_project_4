import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

"""Load data"""

Read energy
E = pa.mat()  Create pa.mat object (just as arma::mat in C++)
E.load("./plot/binary_data/20x20_E_mean.bin")
E = np.array(E)  Convert to numpy array

Read magnetization
M = pa.mat()  Create pa.mat object (just as arma::mat in C++)
M.load("./plot/binary_data/20x20_M_mean.bin")
M = np.array(M)  Convert to numpy array

"""
MATRIX INDICES (Temprature, Initial state)
0: T = 1, Random  |  1: T = 1, Ordered |  2: T = 2.4, Random |  3: T = 2.4, Ordered
"""
n_cycles = len(E[0])
label_fontsize = 16
ticks_fontsize = 16
legend_fontsize = 16


"""Plot"""

# ENERGY PLOT T = 1
fig = plt.figure(figsize=(6, 4.5))

plt.plot(
    np.arange(n_cycles),
    E[0],
    "-",
    color="#2c1629",
    linewidth=2.0,
    alpha=0.8,
    label=r"T = 1.0 $J / k_B$, Random",
)
plt.plot(
    np.arange(n_cycles),
    E[1],
    "-",
    color="#4d8229",
    linewidth=2.0,
    alpha=0.8,
    label=r"T = 1.0 $J / k_Ordered$, B",
)

plt.xlabel("Number of Monte Carlo cycles", fontsize=label_fontsize)
plt.ylabel(r"$\langle \epsilon \rangle$ [$J$]", fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize, rotation=45)
plt.yticks(fontsize=ticks_fontsize)
plt.legend(prop={"size": legend_fontsize})
plt.grid()
ax = plt.gca()
ax.set_facecolor("#e6e6e6")
plt.tight_layout()

plt.savefig("plot/figures/20x20_1_E.pdf")

# ENERGY PLOT T = 2.4
fig = plt.figure(figsize=(6, 4.5))
plt.plot(
    np.arange(n_cycles),
    E[2],
    "-",nn
    color="#2c1629",
    linewidth=2.0,
    alpha=0.8,
    label=r"T = 2.4 $J / k_B$, Random",
)
plt.plot(
    np.arange(n_cycles),
    E[3],
    "-",
    color="#4d8229",
    linewidth=2.0,
    alpha=0.8,
    label=r"T = 2.4 $J / k_B$, Ordered",
)

plt.xlabel("Number of Monte Carlo cycles", fontsize=label_fontsize)
plt.ylabel(r"$\langle \epsilon \rangle$ [$J$]", fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize, rotation=45)
plt.yticks(fontsize=ticks_fontsize)
plt.legend(prop={"size": legend_fontsize})
plt.grid()
ax = plt.gca()
ax.set_facecolor("#e6e6e6")
plt.tight_layout()


plt.savefig("plot/figures/20x20_24_E.pdf")


# MAGNETIZATION PLOT T = 1
fig = plt.figure(figsize=(6, 4.5))

plt.plot(
    np.arange(n_cycles),
    M[0],
    "-",
    color="#2c1629",
    linewidth=2.0,
    alpha=0.8,
    label=r"T = 1.0 $J / k_B$, Random",
)
plt.plot(
    np.arange(n_cycles),
    M[1],
    "-",
    color="#4d8229",
    linewidth=2.0,
    alpha=0.8,
    label=r"T = 1.0 $J / k_B$, Ordered",
)
plt.xlabel("Number of Monte Carlo cycles", fontsize=label_fontsize)
plt.ylabel(r"$\langle |m| \rangle$ [$-$]", fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize, rotation=45)
plt.yticks(fontsize=ticks_fontsize)
plt.legend(prop={"size": legend_fontsize})
plt.grid()
# Set background color to gray
ax = plt.gca()
ax.set_facecolor("#e6e6e6")
plt.tight_layout()


plt.savefig("plot/figures/20x20_1_M.pdf")

# MAGNETIZATION PLOT T = 2.4
fig = plt.figure(figsize=(6, 4.5))
plt.plot(
    np.arange(n_cycles),
    M[2],
    "-",
    color="#2c1629",
    linewidth=2.0,
    alpha=0.8,
    label=r"T = 2.4 $J / k_B$, Random",
)

plt.plot(
    np.arange(n_cycles),
    M[3],
    "-",
    color="#4d8229",
    linewidth=2.0,
    alpha=0.8,
    label=r"T = 2.4 $J / k_B$, Ordered",
)

plt.xlabel("Number of Monte Carlo cycles", fontsize=label_fontsize)
plt.ylabel(r"$\langle |m|\rangle$ [$-$]", fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize, rotation=45)
plt.yticks(fontsize=ticks_fontsize)
plt.legend(prop={"size": legend_fontsize})
plt.grid()
ax = plt.gca()
ax.set_facecolor("#e6e6e6")
plt.tight_layout()


plt.savefig("plot/figures/20x20_24_M.pdf")

# ENERGY PER SPIN
e_1 = pa.mat()
e_24 = pa.mat()
e_1.load("./plot/binary_data/20x20_1_e.bin")
e_24.load("./plot/binary_data/20x20_24_e.bin")
e_1 = np.array(e_1)
e_24 = np.array(e_24)

# Energy histogram T = 1
counts, bins = np.histogram(e_1, bins=4)
fig = plt.figure(figsize=(6, 4.5))
plt.stairs(counts / counts.sum(), bins, color="#4d8229")
plt.xlabel(r"$\langle \epsilon \rangle$ [$J$]", fontsize=label_fontsize)
plt.ylabel(r"Estimation of $p_{\epsilon}(\epsilon;1.0)$", fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize, rotation=45)
plt.yticks(fontsize=ticks_fontsize)
ax = plt.gca()
ax.set_facecolor("#e6e6e6")
plt.tight_layout()


plt.savefig("plot/figures/20x20_1_HIST.pdf")

# Energy histogram T = 2.4
counts, bins = np.histogram(e_24, bins=400)
fig = plt.figure(figsize=(6, 4.5))
plt.stairs(counts / counts.sum(), bins, color="#4d8229")
plt.xlabel(r"$\langle \epsilon \rangle$ [$J$]", fontsize=legend_fontsize)
plt.ylabel(r"Estimation of $p_{\epsilon}(\epsilon;2.5)$", fontsize=legend_fontsize)
plt.xticks(fontsize=ticks_fontsize, rotation=45)
plt.yticks(fontsize=ticks_fontsize)
ax = plt.gca()
ax.set_facecolor("#e6e6e6")
plt.tight_layout()

plt.savefig("plot/figures/20x20_24_HIST.pdf")
plt.show()
