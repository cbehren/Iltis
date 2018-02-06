import matplotlib.pyplot as plt
import numpy as np
def process_file(fname,bins):
    x=np.loadtxt(fname)
    H,_ = np.histogram(x,bins=bins,density=True)
    return H
bins = np.linspace(-10,10,200)
colors = ["k","b","r","g","y","c"]
for i in range(0,6):
    H=process_file("x_"+repr(i)+".txt",bins)
    plt.plot(bins[0:-1],H,label="x="+repr(i),c=colors[i])
    H=process_file("comparison/redistribution_test_x="+repr(i),bins)
    plt.plot(bins[0:-1],H,label="old x="+repr(i),c=colors[i],linestyle="--")

plt.legend(loc="center left")
plt.savefig("redist.pdf")
plt.show()
