import numpy as np
import matplotlib.pyplot as plt


HO_E = np.loadtxt('VMC HO/E.dat')
Hy_E = np.loadtxt('VMC Hy/E.dat')
He_E = np.loadtxt('VMC He/E.dat')
H2_E = np.loadtxt('VMC H2/E.dat')

HO_N, HO_steps, HO_eq_steps, HO_N_alpha = np.loadtxt('VMC HO/data_int.dat') #N, steps, eq_steps, N_alpha
HO_alpha_min, HO_alpha_max = np.loadtxt('VMC HO/data_float.dat')
Hy_N, Hy_steps, Hy_eq_steps, Hy_N_alpha = np.loadtxt('VMC Hy/data_int.dat') 
Hy_alpha_min, Hy_alpha_max = np.loadtxt('VMC Hy/data_float.dat')
He_N, He_steps, He_eq_steps, He_N_alpha = np.loadtxt('VMC He/data_int.dat') 
He_alpha_min, He_alpha_max = np.loadtxt('VMC He/data_float.dat') 
H2_N, H2_steps, H2_eq_steps, H2_N_alpha = np.loadtxt('VMC H2/data_int.dat') 
H2_alpha_min, H2_alpha_max = np.loadtxt('VMC H2/data_float.dat') 

HO_alphas = np.linspace(HO_alpha_min,HO_alpha_max,HO_N_alpha)
Hy_alphas = np.linspace(Hy_alpha_min,Hy_alpha_max,Hy_N_alpha)
He_alphas = np.linspace(He_alpha_min,He_alpha_max,He_N_alpha)
H2_alphas = np.linspace(H2_alpha_min,H2_alpha_max,H2_N_alpha)

######## just for HO
HO_E_analytic = lambda alpha: 0.5*alpha+1/(8*alpha)
HO_error_analytic = lambda alpha: (1-4*alpha**2)**2/(32*alpha**2)
HO_alphas_analytic = np.linspace(HO_alpha_min,HO_alpha_max, 100)
######## 


E_avg_HO = np.mean(HO_E, axis=0)
E_error_HO = np.var(HO_E, axis=0)*HO_N

E_avg_Hy = np.mean(Hy_E, axis=0)
E_error_Hy = np.var(Hy_E, axis=0)*Hy_N

E_avg_He = np.mean(He_E, axis=0)
E_error_He = np.var(He_E, axis=0)*He_N

E_avg_H2 = np.mean(H2_E, axis=0)
E_error_H2 = np.var(H2_E, axis=0)*H2_N



HO_min_x_idx = np.argmin(E_avg_HO)
HO_min_x = HO_alphas[HO_min_x_idx]
HO_min_y = E_avg_HO[HO_min_x_idx]

Hy_min_x_idx = np.argmin(E_avg_Hy)
Hy_min_x = Hy_alphas[Hy_min_x_idx]
Hy_min_y = E_avg_Hy[Hy_min_x_idx]

He_min_x_idx = np.argmin(E_avg_He)
He_min_x = He_alphas[He_min_x_idx]
He_min_y = E_avg_He[He_min_x_idx]

H2_min_x_idx = np.argmin(E_avg_H2)
H2_min_x = H2_alphas[H2_min_x_idx]
H2_min_y = E_avg_H2[H2_min_x_idx]

print 'HarmOs min', HO_min_x, HO_min_y
print 'Hydrogen min', Hy_min_x, Hy_min_y
print 'Helium min', He_min_x, He_min_y
print 'H2 min', H2_min_x, H2_min_y


def plot():
 plt.subplot(3,2,1)
 plt.ylabel(r'Energy $\left\langle E \right\rangle$')
 plt.xlabel(r'$\alpha$')
 plt.plot(HO_alphas, E_avg_HO, 'r-p', label=r'$\left\langle E \right\rangle$ for a Harmonic oscillator')
 plt.plot(HO_alphas_analytic, HO_E_analytic(HO_alphas_analytic), 'g', label=r'Analytic function for $\left\langle E \right\rangle$')
 plt.plot(HO_min_x,HO_min_y, 'ro',markersize=10)
 #plt.errorbar(HO_alphas,E_avg_HO, E_error_HO)
 plt.legend()

 plt.subplot(3,2,2)
 plt.ylabel(r'Variance $\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$')
 plt.xlabel(r'$\alpha$')
 plt.ylim(0,0.006)
 plt.plot(HO_alphas, E_error_HO, 'r-p', label=r'$\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$ for a Harmonic oscillator')
 plt.plot(HO_alphas_analytic, HO_error_analytic(HO_alphas_analytic), 'g', label=r'Analytic function for $\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$')
 plt.legend()



 plt.subplot(3,2,3)
 plt.ylabel(r'Energy $\left\langle E \right\rangle$')
 plt.xlabel(r'$\alpha$')
 plt.plot(Hy_alphas, E_avg_Hy, 'r', label=r'$\left\langle E \right\rangle$ for a hydrogen atom')
 plt.plot(Hy_min_x,Hy_min_y, 'ro',markersize=10) 
 plt.legend()

 plt.subplot(3,2,4)
 plt.ylabel(r'Variance $\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$')
 plt.xlabel(r'$\alpha$')
 plt.plot(Hy_alphas, E_error_Hy, 'r', label=r'$\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$ for a hydrogen atom')
 plt.legend()



 plt.subplot(3,2,5)
 plt.ylabel(r'Energy $\left\langle E \right\rangle$')
 plt.xlabel(r'$\alpha$')
 plt.plot(He_alphas, E_avg_He, 'r', label=r'$\left\langle E \right\rangle$ for a helium atom')
 plt.plot(He_min_x,He_min_y, 'ro',markersize=10)
 plt.legend()

 plt.subplot(3,2,6)
 plt.ylabel(r'Variance $\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$')
 plt.xlabel(r'$\alpha$')
 plt.plot(He_alphas, E_error_He, 'r', label=r'$\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$ for a helium atom')
 plt.legend()
 plt.show()


def plot_HO():
 plt.subplot(3,1,1)
 plt.title('Number of walkers: %s, MC steps: %s'%(int(HO_N), int(HO_steps-HO_eq_steps)))
 plt.ylabel(r'Energy $\left\langle E \right\rangle$')
 plt.xlabel(r'$\alpha$')
 plt.plot(HO_alphas, E_avg_HO, 'r-p', label=r'$\left\langle E \right\rangle$ for a Harmonic oscillator')
 plt.plot(HO_alphas_analytic, HO_E_analytic(HO_alphas_analytic), 'g', label=r'Analytic function for $\left\langle E \right\rangle$')
 plt.plot(HO_min_x,HO_min_y, 'ro',markersize=10)
 #plt.errorbar(HO_alphas,E_avg_HO, E_error_HO)
 plt.legend()

 plt.subplot(3,1,2)
 plt.ylabel(r'Variance $\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$')
 plt.xlabel(r'$\alpha$')
 plt.ylim(0,0.006)
 plt.plot(HO_alphas, E_error_HO, 'r-p', label=r'$\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$ for a Harmonic oscillator')
 plt.plot(HO_alphas_analytic, HO_error_analytic(HO_alphas_analytic), 'g', label=r'Analytic function for $\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$')
 plt.legend()

 plt.subplot(3,1,3)
 plt.ylabel(r'Energy $\left\langle E \right\rangle$')
 plt.xlabel(r'$\alpha$')
 plt.plot(HO_alphas, E_avg_HO, 'r-p', label=r'$\left\langle E \right\rangle$ for a Harmonic oscillator')
 plt.plot(HO_alphas_analytic, HO_E_analytic(HO_alphas_analytic), 'g', label=r'Analytic function for $\left\langle E \right\rangle$')
 plt.plot(HO_min_x,HO_min_y, 'ro',markersize=10)
 plt.errorbar(HO_alphas,E_avg_HO, E_error_HO)
 plt.legend()
 plt.show()

def plot_Hy():
 plt.subplot(3,1,1)
 plt.title('Number of walkers: %s, MC steps: %s'%(int(Hy_N), int(Hy_steps-Hy_eq_steps)))
 plt.ylabel(r'Energy $\left\langle E \right\rangle$')
 plt.xlabel(r'$\alpha$')
 plt.plot(Hy_alphas, E_avg_Hy, 'r-p', label=r'$\left\langle E \right\rangle$ for a hydrogen atom')
 plt.plot(Hy_min_x,Hy_min_y, 'ro',markersize=10)
 plt.legend()

 plt.subplot(3,1,2)
 plt.ylabel(r'Variance $\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$')
 plt.xlabel(r'$\alpha$')
 plt.plot(Hy_alphas, E_error_Hy, 'r-p', label=r'$\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$ for a hydrogen atom')
 plt.legend()

 plt.subplot(3,1,3)
 plt.ylabel(r'Energy $\left\langle E \right\rangle$')
 plt.xlabel(r'$\alpha$')
 plt.plot(Hy_alphas, E_avg_Hy, 'r-p', label=r'$\left\langle E \right\rangle$ for a hydrogen atom')
 plt.plot(Hy_min_x,Hy_min_y, 'ro',markersize=10)
 plt.errorbar(Hy_alphas,E_avg_Hy, E_error_Hy)
 plt.legend()
 plt.show()

def plot_He():
 plt.subplot(3,1,1)
 plt.title('Number of walkers: %s, MC steps: %s'%(int(He_N), int(He_steps-He_eq_steps)))
 plt.ylabel(r'Energy $\left\langle E \right\rangle$')
 plt.xlabel(r'$\alpha$')
 plt.plot(He_alphas, E_avg_He, 'r-p', label=r'$\left\langle E \right\rangle$ for a helium atom')
 plt.plot(He_min_x,He_min_y, 'ro',markersize=10)
 plt.legend()

 plt.subplot(3,1,2)
 plt.ylabel(r'Variance $\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$')
 plt.xlabel(r'$\alpha$')
 plt.plot(He_alphas, E_error_He, 'r-p', label=r'$\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$ for a helium atom')
 plt.legend()

 plt.subplot(3,1,3)
 plt.ylabel(r'Energy $\left\langle E \right\rangle$')
 plt.xlabel(r'$\alpha$')
 plt.plot(He_alphas, E_avg_He, 'r-p', label=r'$\left\langle E \right\rangle$ for a helium atom')
 plt.plot(He_min_x,He_min_y, 'ro',markersize=10)
 plt.errorbar(He_alphas,E_avg_He, E_error_He)
 plt.legend()
 plt.show()

def plot_H2():
 plt.subplot(3,1,1)
 plt.title('Number of walkers: %s, MC steps: %s'%(int(H2_N), int(H2_steps-H2_eq_steps)))
 plt.ylabel(r'Energy $\left\langle E \right\rangle$')
 plt.xlabel(r'$\alpha$')
 plt.plot(H2_alphas, E_avg_H2, 'r-p', label=r'$\left\langle E \right\rangle$ for a hydrogen molecule (H$_2$)')
 plt.plot(H2_min_x,H2_min_y, 'ro',markersize=10)
 plt.legend()

 plt.subplot(3,1,2)
 plt.ylabel(r'Variance $\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$')
 plt.xlabel(r'$\alpha$')
 plt.plot(H2_alphas, E_error_H2, 'r-p', label=r'$\left\langle E^2 \right\rangle - \left\langle E \right\rangle^2$ for a hydrogen molecule (H$_2$)')
 plt.legend()

 plt.subplot(3,1,3)
 plt.ylabel(r'Energy $\left\langle E \right\rangle$')
 plt.xlabel(r'$\alpha$')
 plt.plot(H2_alphas, E_avg_H2, 'r-p', label=r'$\left\langle E \right\rangle$ for a hydrogen molecule (H$_2$)')
 plt.plot(H2_min_x,H2_min_y, 'ro',markersize=10)
 plt.errorbar(H2_alphas,E_avg_H2, E_error_H2)
 plt.legend()
 plt.show()
plot()
