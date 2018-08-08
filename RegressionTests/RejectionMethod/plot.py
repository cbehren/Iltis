import numpy as np
from numpy import pi, exp, sqrt, log
import matplotlib.pyplot as plt

#plot the results from the rejection method versus the analytic expectation for a number of input frequencies.

def get_voigt(x,a): #following Laursen et al 2009, 08050.3153
    q=0
    z=(x*x-0.855)/(x*x+3.42)
    if(z>0):
        q=(1.0+21.0/(x*x))*(a/(pi*(x*x+1)))*(5.674*(z*z*z*z)-9.207*(z*z*z)+4.421*(z*z)+0.1117*z);
    return (q*sqrt(pi)+exp(-(x*x)));
def pdf_u_par(u,x,a):#the analytic PDF
    return a/pi*(exp(-u**2))/((x-u)**2+a**2)/get_voigt(x,a)
def get_sigma(T):
    local_v_thermal = 12.85*sqrt(T/1e4)
    a=4.7e-4*(12.85/(local_v_thermal))
    return a


colors = ["k","b","g","r","m","y"]

for i,x in enumerate([0.0,1,2,3,4,5]):
    a=get_sigma(2.74)
    us=np.linspace(-5,5,900)
    #load the data
    fname = "rejection_method_x={:1.1f}.txt".format(x)
    data = np.loadtxt(fname)
    H,xedges = np.histogram(data,bins=us,density=True)
    centers = xedges[0:-1]+(xedges[1]-xedges[0])/2.
    ps = pdf_u_par(centers,x,a)
    plt.plot(centers,H,label="numerical, x="+repr(x),c=colors[i],alpha=0.3)
    plt.plot(centers,ps,label="analytic",c=colors[i],ls="dashed")
    #plt.xlim(-2.5,-0.)
    plt.ylim(1e-4,30)
    plt.legend()
    plt.yscale("log")
plt.show()