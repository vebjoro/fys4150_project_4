import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa
import scipy.stats as st


"""Load data"""
T = pa.mat()
T.load("./plot/binary_data/OMP_T_out.bin")
T = np.array(T).T[0]  # Convert to numpy array and transpose


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


Tc = []
for i in range(4):
    k  = np.argmax(Cv[i])
    l = np.argmax(X[i])
    Tc.append((T[k] +T[l])/2)



L = np.array([40, 60, 80, 100])
# """Plot"""
regress = st.linregress(1/L, Tc)
Tc_inf = regress[1]
Tc_analytic = 2/(np.log(1 + np.sqrt(2)))
#print
relerr = abs((Tc_analytic - Tc_inf)/Tc_analytic)
print(f"Avereged Tc's: {Tc}")
print(f"a-value: {regress[0]}")
print(f"Tc-estimate: {Tc_inf}")
print(f"Analytic Tc: {Tc_analytic}")
print(f"Relative error: {relerr}")

x_min, x_max = -0.005, 0.03
x = np.linspace(x_min, x_max)
y = regress[0]*x + Tc_inf

 # fig = plt.figure()
# ax = plt.axes()
# ax.plot(1/L, Tc, 'x', color = 'r')
# ax.plot(x, y, '--', color  = 'k')
# ax.axhline(Tc_inf, color = 'k')
# plt.show()
