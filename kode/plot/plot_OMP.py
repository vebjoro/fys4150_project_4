import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

"""Load data"""

# Temperature first scan
T1 = pa.mat()
T1.load("./plot/binary_data/OMP_T_out.bin")
T1 = np.array(T1)

# Energy first scan
e1 = pa.mat()
e1.load("./plot/binary_data/OMP_e_out.bin")
e1 = np.array(e1)

# Magnetization first scan
m1 = pa.mat()
m1.load("./plot/binary_data/OMP_m_out.bin")
m1 = np.array(m1)

# Heat capacity first scan
Cv1 = pa.mat()
Cv1.load("./plot/binary_data/OMP_Cv_out.bin")
Cv1 = np.array(Cv1)

# Susceptibility first scan
X1 = pa.mat()
X1.load("./plot/binary_data/OMP_X_out.bin")
X1 = np.array(X1)

"""____________________________________________________________________________"""

# Temperature second scan
T2 = pa.mat()
T2.load("./plot/valuable_data/OMP_T_out.bin")
T2 = np.array(T2)  # .T[0]  # Convert to numpy array and transpose

# Energy second scan
e2 = pa.mat()
e2.load("./plot/valuable_data/OMP_e_out.bin")
e2 = np.array(e2)

# Magnetization second scan
m2 = pa.mat()
m2.load("./plot/valuable_data/OMP_m_out.bin")
m2 = np.array(m2)

# Heat capacity second scan
Cv2 = pa.mat()
Cv2.load("./plot/valuable_data/OMP_Cv_out.bin")
Cv2 = np.array(Cv2)

# Susceptibility second scan
X2 = pa.mat()
X2.load("./plot/valuable_data/OMP_X_out.bin")
X2 = np.array(X2)

# Merge the two scans
T = np.append(T1, T2)
z = np.argsort(T)
T = T[z]

sorted_e = []
for i in range(len(e1)):
    sorted_e.append(np.append(e1[i], e2[i])[z])
e = np.array(sorted_e)

sorted_m = []
for i in range(len(m1)):
    sorted_m.append(np.append(m1[i], m2[i])[z])
m = np.array(sorted_m)

sorted_Cv = []
for i in range(len(Cv1)):
    sorted_Cv.append(np.append(Cv1[i], Cv2[i])[z])
Cv = np.array(sorted_Cv)

sorted_X = []
for i in range(len(X1)):
    sorted_X.append(np.append(X1[i], X2[i])[z])
X = np.array(sorted_X)


"""Plot"""

# ENERGY PLOT
fig = plt.figure(figsize=(6, 4.5))
plt.plot(T, e[0], "-", color="#8a1629", linewidth=2.0, alpha=0.8, label="40x40")

plt.plot(T, e[1], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label="60x60")
plt.plot(T, e[2], "-", color="#2c1629", linewidth=2.0, alpha=0.8, label="80x80")
plt.plot(T, e[3], "-", color="#4d8229", linewidth=2.0, alpha=0.8, label="100x100")

plt.legend()
plt.xlabel(r"T [$J/k_B$]")
plt.ylabel("Energy, expected values")
plt.grid()
ax = plt.gca()
ax.set_facecolor("#e6e6e6")
plt.savefig("plot/figures/OMP_E.pdf")

# MAGNETIZATION PLOT
fig = plt.figure(figsize=(6, 4.5))
plt.plot(T, m[0], "-", color="#8a1629", linewidth=2.0, alpha=0.8, label="40x40")
plt.plot(T, m[1], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label="60x60")
plt.plot(T, m[2], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label="80x80")
plt.plot(T, m[3], "-", color="#4d8229", linewidth=2.0, alpha=0.8, label="100x100")
plt.legend()
plt.xlabel(r"T [$J/k_B$]")
plt.ylabel("Magnetization, expected values")
plt.grid()
ax = plt.gca()
ax.set_facecolor("#e6e6e6")
plt.savefig("plot/figures/OMP_M.pdf")

# HEAT CAPACITY PLOT
fig = plt.figure(figsize=(6, 4.5))
plt.plot(T, Cv[0], "-", color="#8a1629", linewidth=2.0, alpha=0.8, label="40x40")
plt.plot(T, Cv[1], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label="60x60")
plt.plot(T, Cv[2], "-", color="#2c1629", linewidth=2.0, alpha=0.8, label="80x80")
plt.plot(T, Cv[3], "-", color="#4d8229", linewidth=2.0, alpha=0.8, label="100x100")
plt.legend()
plt.xlabel(r"$T$ [$J/k_B$]")
plt.ylabel(r"$C_V$")
plt.grid()
ax = plt.gca()
ax.set_facecolor("#e6e6e6")
plt.savefig("plot/figures/OMP_Cv.pdf")

# SUSCEPTIBILITY PLOT
fig = plt.figure(figsize=(6, 4.5))
plt.plot(T, X[0], "-", color="#8a1629", linewidth=2.0, alpha=0.8, label="40x40")
plt.plot(T, X[1], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label="60x60")
plt.plot(T, X[2], "-", color="#2c1629", linewidth=2.0, alpha=0.8, label="80x80")
plt.plot(T, X[3], "-", color="#4d8229", linewidth=2.0, alpha=0.8, label="100x100")
plt.legend()
plt.xlabel(r"T [$J/k_B$]")
plt.ylabel(r"$\chi$")
plt.grid()
ax = plt.gca()
ax.set_facecolor("#e6e6e6")
plt.savefig("plot/figures/OMP_X.pdf")
