import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as c
from astropy import units as u
import h5py

def lambda2x(lam,T):
    from astropy import constants as c
    from astropy import units as u
    vthermal = 12.85*np.sqrt(T/10000.0)*u.km/u.s
    c_val = c.c
    Lambda_0 = 1215.668*u.angstrom
    Nu_0 = (c_val/Lambda_0).to(u.Hz);
    freq = c_val/((lam*u.angstrom+Lambda_0))-Nu_0;
    return (freq*c_val/Nu_0/vthermal).to(u.dimensionless_unscaled);

names = ["Output.txt.hdf5","peeling.txt_0.hdf5","compare.hdf5"]
labels = ["Iltis","Iltis (peeling off)","reference"]
for name,label in zip(names,labels):
    with h5py.File(name) as f:
        w = f["weight"][:]
        lam = f["lambda"][:]
        if("reference" in label):
            mylam = lam
            w = None #weights are broken in the reference
        else:
            mylam = lambda2x(lam,96897.8)
        
    plt.hist(mylam,100,histtype="step",label=label,normed=True,weights=w)
plt.legend()
plt.xlabel("dimensionless frequency x")
plt.ylabel("normalized flux")
plt.savefig("shell.png")
#plt.show()
