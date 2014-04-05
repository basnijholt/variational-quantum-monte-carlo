import numpy as np
import matplotlib.pyplot as plt


HO_E = np.loadtxt('VMC H2/E.dat')
HO_data = np.loadtxt('VMC H2/data.dat') #N, steps, eq_steps, N_alpha, alpha_min, alpha_max
HO_N = HO_data[0]
HO_steps = HO_data[1]
HO_eq_steps = HO_data[2]
HO_N_alpha = HO_data[3]
HO_alpha_min = HO_data[4]
HO_alpha_max = HO_data[5]
HO_alphas = np.linspace(HO_alpha_min,HO_alpha_max,HO_N_alpha)


E_avg_HO = np.mean(HO_E, axis=0)
E_error_HO = np.var(HO_E, axis=0)*HO_N




HO_min_x_idx = np.argmin(E_avg_HO)
HO_min_x = HO_alphas[HO_min_x_idx]
HO_min_y = E_avg_HO[HO_min_x_idx]


print 'HarmOs min', HO_min_x, HO_min_y

def plot():
 plt.subplot(3,1,1)
 plt.ylabel(r'Energy $\left\langle E \right\rangle$')
 plt.xlabel(r'$\alpha$')
 plt.plot(HO_alphas, E_avg_HO, 'r-p', label=r'$\left\langle E \right\rangle$ for a Harmonic oscillator')
 plt.plot(HO_min_x,HO_min_y, 'ro',markersize=10)
 #plt.errorbar(HO_alphas,E_avg_HO, E_error_HO)
 plt.legend()

 plt.subplot(3,1,2)
 plt.ylabel(r'Variance $\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$')
 plt.xlabel(r'$\alpha$')
 plt.plot(HO_alphas, E_error_HO, 'r-p', label=r'$\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$ for a Harmonic oscillator')
 plt.legend()

 plt.subplot(3,1,3)
 plt.ylabel(r'Energy $\left\langle E \right\rangle$')
 plt.xlabel(r'$\alpha$')
 plt.plot(HO_alphas, E_avg_HO, 'r-p', label=r'$\left\langle E \right\rangle$ for a Harmonic oscillator')
 plt.plot(HO_min_x,HO_min_y, 'ro',markersize=10)
 plt.errorbar(HO_alphas,E_avg_HO, E_error_HO)
 plt.legend()
 plt.show()

plot()
