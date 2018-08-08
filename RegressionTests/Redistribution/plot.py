import matplotlib.pyplot as plt
import numpy as np
import sys
    
#this script plots the analytic redistribution function RII(xin,xout) from Hummer 1962 versus the numerical results.


#routines to calculate the analytic solution from Hummer 1962, eq. 3.12.1
from numpy import sqrt,arctan,exp,pi
import scipy.integrate as integrate
def fint(u,params):
    lx = params[0]
    hx = params[1]
    sigma = params[2]
    return exp(-u**2)*(arctan((lx+u)/sigma)-arctan((hx-u)/sigma))
def P(x,xprime,sigma):
    lx = min(x,xprime)
    hx = max(x,xprime)
    lint = 0.5*abs(hx-lx)
    integral = integrate.quad(fint,lint,1000.,args=[lx,hx,sigma],epsabs=1.0e-20, limit=150)
    value = pi**(-3./2.)*integral[0]
    return value
def get_sigma(T):
    local_v_thermal = 12.85*sqrt(T/1e4)
    a=4.7e-4*(12.85/(local_v_thermal))
    return a
def pdf_data(x,sigma,bins=np.linspace(-10,10,200)):
    y = np.zeros_like(bins)
    for i,xprime in enumerate(bins):
        y[i] = P(x,xprime,sigma)
    return bins,y/(y.sum()*(bins[1]-bins[0]))


#routines to calculate the analytic solution from Hummer 1962, eq. 3.12.2 for the dipole phase function.

#the solution to the inner integral of 3.12.2 (the integral over t)
def FINT_INNER(t,params):
    sigma = params[0]
    x = params[1]
    p = params[2]
    u = params[3]
    
    l1 = (t-x)*(3*p**2-3*p*t+9*p*x-3*sigma**2+t**2-2*t*x-2*u**2+x**2)
    l2 = np.log(sigma**2+t**2)*(-3*p**2*x+p*(3*sigma**2+u**2-3*x**2)+x*(3*sigma**2+u**2))
    l3 = np.arctan(t/sigma)*(p**2*(3*sigma**2+u**2-3*x**2)
                             +12*p*sigma**2*x-3*u**4+u**2*(x**2-2*sigma**2)+3*sigma**2*(x**2-sigma**2))
    r = 1./sigma/u**4*(sigma*(l1+l2)-l3)
    return r
    
def fint_outer(u,params):
    lx = params[0]
    hx = params[1]
    sigma = params[2]
    params2 = [params[2],params[3],params[4],u]
    #return exp(-u**2)*integrate.quad(fint_inner,hx-u,lx+u,args=params2)[0]
    return exp(-u**2)*(FINT_INNER(lx+u,params2)-FINT_INNER(hx-u,params2))
def Pdipole(x,xprime,sigma):
    lx = min(x,xprime)
    hx = max(x,xprime)
    lint = 0.5*abs(hx-lx)
    integral = integrate.quad(fint_outer,lint,1000.,args=[lx,hx,sigma,x,xprime])
    value = 3*pi**(-3./2.)/8*sigma*integral[0]
    return value
def get_sigma(T):
    local_v_thermal = 12.85*sqrt(T/1e4)
    a=4.7e-4*(12.85/(local_v_thermal))
    return a
def pdf_data_dipole(x,sigma,bins=np.linspace(-10,10,200)):
    y = np.zeros_like(bins)
    for i,xprime in enumerate(bins):
        y[i] = Pdipole(x,xprime,sigma)
    return bins,y/(y.sum()*(bins[1]-bins[0]))



def process_file(fname,bins):
    x=np.loadtxt(fname)
    H,_ = np.histogram(x,bins=bins,density=True)
    return H

use_dipole = False
if(len(sys.argv)>1 and "dipole" in sys.argv[1]):
    use_dipole = True
    print "Plotting data for dipole scattering. Might take a while."


bins = np.linspace(-4,10,400)
sigma = get_sigma(2.74)
colors = ["k","b","r","g","y","c"]
for i in range(0,6):
    if(use_dipole):
        fname = "x_dipole_"+repr(i)+".txt"
        label = "(dipole) x="+repr(i)
        x,y = pdf_data_dipole(float(i),sigma,bins=bins)
    else:    
        fname = "x_"+repr(i)+".txt"
        label = "x="+repr(i)
        x,y = pdf_data(float(i),sigma,bins=bins)
    H=process_file(fname,bins)
    centers = bins[0:-1]+(bins[1]-bins[0])/2.
    plt.plot(centers,H,label=label,c=colors[i])
    
    plt.plot(x,y,label="analytic x="+repr(i),c=colors[i],linestyle="--")
    plt.xlabel("x\'")
    plt.ylabel("RII(x,x\')")

plt.ylim(0,1.)
plt.legend(loc="center left")
if(use_dipole):
    plt.savefig("redistributon_dipole.png")
else:
    plt.savefig("redistributon.png")
plt.show()
