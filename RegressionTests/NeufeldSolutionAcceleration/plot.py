
import numpy as np
import matplotlib.pyplot as plt
def analytic_slab(x,a,tau0,a_scale=1.0):
    #analytic solution as obtained by Neufeld+90
    #it has an optional parameter a_scale to scale a*tau0 according to our needs.
    atau0 = a_scale*a*tau0
    factor = np.sqrt(6./np.pi)/24./(atau0)
    value = factor*(x**2./(1.+np.cosh(np.sqrt(np.pi**3./54.)*np.abs(x**3.)/(atau0))))*np.pi
    #print value
    return value
def a_factor(T):
    vthermal = 12.85*np.sqrt(T/10000.0)
    return 4.7e-4*(12.85/vthermal)
a = a_factor(2e4)
tau0 = 1e7
v = np.loadtxt("drawn_frequencies.txt")
x = np.linspace(-50,50,100)
H,edges = np.histogram(v,bins=x,density=True)
xh = edges[0:-1]+(edges[1]-edges[0])/2.
plt.plot(xh,H,label="drawn",alpha=0.4)
plt.plot(x,analytic_slab(x,a,tau0,a_scale=0.71)*2*np.pi,label="analytic",alpha=0.4)
plt.legend()
plt.xlabel("frequency x")
plt.ylabel("P(x)")
plt.savefig("neufeld.png")
plt.show()

def mu_analytic(mu):
    #CDF of the distribution of mu = cos theta for direction of exit of a cell, see Tasitsiomi 2006
    return mu**2/7.*(3.+4.*mu)

import numpy as np
import matplotlib.pyplot as plt
v = np.loadtxt("drawn_direction.txt")
x = np.linspace(0,1,100)
H,edges = np.histogram(v,bins=x)
Hcum = np.cumsum(H)/np.sum(H.astype(float))
xh = edges[0:-1]+(edges[1]-edges[0])/2.
plt.plot(xh,Hcum,label="drawn")
plt.plot(x,mu_analytic(x),label="analytic",alpha=0.4)
plt.gca().set_aspect(1.0)
plt.xlabel(r"$\mu$")
plt.ylabel(r"P(<$\mu$)")
plt.savefig("directions.png")
plt.legend()
plt.show()