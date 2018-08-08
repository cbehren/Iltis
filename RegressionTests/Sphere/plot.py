import numpy as np
import matplotlib.pyplot as plt
import h5py

def lambda2x(lam,T):
    #convert wavelenght to dimensionless frequency
    from astropy import constants as c
    from astropy import units as u
    vthermal = 12.85*np.sqrt(T/10000.0)*u.km/u.s
    c_val = c.c
    Lambda_0 = 1215.668*u.angstrom
    Nu_0 = (c_val/Lambda_0).to(u.Hz);
    freq = c_val/((lam*u.angstrom+Lambda_0))-Nu_0;
    return (freq*c_val/Nu_0/vthermal).to(u.dimensionless_unscaled);
def analytic(x,a,tau0):
    #analytic solution as obtained by Dijkstra+2006
    factor = np.sqrt(np.pi/24.)/(a*tau0)
    return factor*(x**2./(1.+np.cosh(np.sqrt(2.*np.pi**3./27.)*np.abs(x**3.)/(a*tau0))))
def a_factor(T):
    vthermal = 12.85*np.sqrt(T/10000.0)
    return 4.7e-4*(12.85/vthermal)

with h5py.File("Output.txt.hdf5") as f:
    w = f["weight"][:]
    lam = f["lambda"][:]
T = 2e4
tau0 = 1e7
lam = lambda2x(lam,T)
H,xedges = np.histogram(lam,weights=w,bins=np.linspace(-50,50,100))
dx = (xedges[1]-xedges[0])
x = xedges[0:-1]+dx/2.
plt.plot(x,H/dx/50000.0/2/np.pi,label="Iltis")
a = a_factor(T)
print "a is",a
print "a*tau0 is",a*tau0
plt.plot(x,analytic(x,a,tau0),label="analytic")
plt.legend()
plt.xlabel("dimensionless frequency x")
plt.ylabel("normalized flux")
plt.savefig("sphere.png")
plt.show()
