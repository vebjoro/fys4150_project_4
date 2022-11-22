import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa
import scipy.stats as st


"""Load data"""
T = pa.mat()
T.load("./plot/c3/data/OMP_T_out.bin")
T = np.array(T).T[0]  # Convert to numpy array and transpose

e = pa.mat()
e.load("./plot/c3/data/OMP_e_out.bin")
e = np.array(e)

m = pa.mat()
m.load("./plot/c3/data/OMP_m_out.bin")
m = np.array(m)

Cv = pa.mat()
Cv.load("./plot/c3/data/OMP_Cv_out.bin")
Cv = np.array(Cv)

X = pa.mat()
X.load("./plot/c3/data/OMP_X_out.bin")
X = np.array(X)

# Find Tc
Tc = []
for i in range(4):
    k = np.argmax(Cv[i])
    l = np.argmax(X[i])
    Tc.append((T[k] + T[l]) / 2)


L = np.array([40, 60, 80, 100])
# """Plot"""

# Fontsize
label_fontsize = 16
ticks_fontsize = 16
legend_fontsize = 16

regress = st.linregress(1 / L, Tc)
Tc_inf = regress[1]
Tc_analytic = 2 / (np.log(1 + np.sqrt(2)))

relerr = abs((Tc_analytic - Tc_inf) / Tc_analytic)
print(f"Avereged Tc's: {Tc}")
print(f"a-value: {regress[0]}")
print(f"Tc-estimate: {Tc_inf}")
print(f"Analytic Tc: {Tc_analytic}")
print(f"Relative error: {relerr}")

x_min, x_max = -0.005, 0.03
x = np.linspace(x_min, x_max)
y = regress[0] * x + Tc_inf

fig = plt.figure(figsize=(6, 4.5))
ax = fig.gca()
ax.plot(1 / L, Tc, "x", color="r")
ax.plot(x, y, "--", color="k")
ax.set_xlabel(r"$\frac{1}{L}$ [$-$]", fontsize=label_fontsize)
ax.set_ylabel(r"$T$ [$\frac{J}{k_B}$]", fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize, rotation=45)
plt.yticks(fontsize=ticks_fontsize)
ax.axhline(Tc_inf, color="k")
plt.tight_layout()
plt.savefig("./plot/figures/plot_Tc.pdf")
