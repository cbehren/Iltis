import sys
import numpy as np
sys.path.append("..")
from LLTCPlotter import LLTCPlotter
import matplotlib.pyplot as plt

fig_size = np.zeros(2)
fig_size[0] = 12
fig_size[1] = 9
plt.rcParams["figure.figsize"] = fig_size

#setup the plotter
params = {}
los = np.array([1.,0.,0.])
center = np.array([0.5,0.5,0.5])
depth = 1.0 #the total depth of the projection, i.e. we go from z = center-depth/2. to z = center+depth/2.0 along the line of sight.
width = 1.0 #total width of the image in code units
params["los"] = los
params["center"] = center
params["up"] = np.array([0.,1.,0.])
params["depth"] = depth
params["width"] = width
params["operators"] = ["massWeightedDensity","columnDensity","slice"]
#params["operators"].append("opticalDepthGrid")
#params["opticalDepthGrid.lambda_left"] = -8
#params["opticalDepthGrid.lambda_right"] = 8
#params["opticalDepthGrid.nlambda"] = 5
params["npixels"] = 64 #**2 is the number of pixels
params["oversampling"] = 1 #useless at the moment, but always set it to 1.
fname ="inputs.shell" #the paramter file we are hooking onto.
#l=LLTCPlotter(fname,params,dry=True)#do not run it, print the function call instead.
l=LLTCPlotter(fname,params,dry=False)#run it
l.store("plotter.pickle")#store to pickle
l = LLTCPlotter.from_file("plotter.pickle")#load from pickle
#optical depth grid and slice return a list/a dictionary of variables:
print l.data["slice"].keys()
for dset in l.data:
    if("slice" in dset):
        continue
    data = l.data[dset]
    vmin = None
    vmax = None
    if("Velocity" not in dset):
        data = np.log10(data)
        cmap = "inferno"
    else:
        cmap = 'coolwarm'
        #data -=data.mean()
        vmax = max(abs(data.min()),abs(data.max()))
        vmin = -vmax
        
        vmax = 400
        vmin = -400
    if("massWeightedDensity" in dset):
        vmin = -2
        vmax = 3 
    ext = l.get_extent()
    map_to_plot = np.log10(l.data["massWeightedDensity"])
    plt.imshow(data,origin="lower",cmap=cmap,extent=ext,vmin=vmin,vmax=vmax)
    plt.title(dset)
    plt.colorbar()
    plt.show()


