import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

"""Load data"""
T = pa.mat()
T.load("./plot/binary_data/OMP_T_out.bin")
T = np.array(T) #.T[0]  # Convert to numpy array and transpose
print(T)

e = pa.mat()
e.load("./plot/binary_data/OMP_e_out.bin")
e = np.array(e)

m = pa.mat()
m.load("./plot/binary_data/OMP_m_out.bin")
m = np.array(m)

Cv = pa.mat()
Cv.load("./plot/binary_data/OMP_Cv_out.bin")
Cv = np.array(Cv)

X = pa.mat()
X.load("./plot/binary_data/OMP_X_out.bin")
X = np.array(X)

"""Plot"""
fig = plt.figure(figsize=(6, 4.5))
plt.plot(T, e[0], "-", color="#8a1629", linewidth=2.0, alpha=0.8, label="40x40")
plt.plot(T, e[1], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label="60x60")
plt.plot(T, e[2], "-", color="#2c1629", linewidth=2.0, alpha=0.8, label="80x80")
plt.plot(T, e[3], "-", color="#4d8229", linewidth=2.0, alpha=0.8, label="100x100")
plt.legend()
plt.xlabel("Temperature")
plt.ylabel("Energy, expected values")
plt.grid()
#plt.savefig("plot/figures/OMP_E.pdf")

fig = plt.figure(figsize=(6, 4.5))
plt.plot(T, m[0], "-", color="#8a1629", linewidth=2.0, alpha=0.8, label="40x40")
plt.plot(T, m[1], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label="60x60")
plt.plot(T, m[2], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label="80x80")
plt.plot(T, m[3], "-", color="#4d8229", linewidth=2.0, alpha=0.8, label="100x100")
plt.legend()
plt.xlabel("Temperature")
plt.ylabel("Magnetization, expected values")
plt.grid()
#plt.savefig("plot/figures/OMP_M.pdf")

fig = plt.figure(figsize=(6, 4.5))
plt.plot(T, Cv[0], "-", color="#8a1629", linewidth=2.0, alpha=0.8, label="40x40")
plt.plot(T, Cv[1], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label="60x60")
plt.plot(T, Cv[2], "-", color="#2c1629", linewidth=2.0, alpha=0.8, label="80x80")
plt.plot(T, Cv[3], "-", color="#4d8229", linewidth=2.0, alpha=0.8, label="100x100")
plt.legend()
plt.xlabel("Temperature")
plt.ylabel("Heat capacity, expected values")
plt.grid()
#plt.savefig("plot/figures/OMP_Cv.pdf")

fig = plt.figure(figsize=(6, 4.5))
plt.plot(T, X[0], "-", color="#8a1629", linewidth=2.0, alpha=0.8, label="40x40")
plt.plot(T, X[1], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label="60x60")
plt.plot(T, X[2], "-", color="#2c1629", linewidth=2.0, alpha=0.8, label="80x80")
plt.plot(T, X[3], "-", color="#4d8229", linewidth=2.0, alpha=0.8, label="100x100")
plt.legend()
plt.xlabel("Temperature")
plt.ylabel("Magnetic susceptibility, expected values")
plt.grid()
#plt.savefig("plot/figures/OMP_X.pdf")
