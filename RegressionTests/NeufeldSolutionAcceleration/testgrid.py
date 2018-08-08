#generates a uniform grid for usage with Iltis. 
#the first line of the output file should contain #<gridsize>, with gridsize**3 being the total number of cells in the grid. 
#each line of the input file should contain the following quantities:
#neutral hydrogen density in 1/cm^3
#temperature in K
#velocity x,y,z in km/s
#dust density in g/cm**3

import numpy as np
gridsize = 99
density = 77.310156544 #tau0 = 1e8
density = 15.4620313088 #tau0 = 1e6
dust_density = 0.0
velocity = 0.
temperature = 20000.0
center = [0.5,0.5,0.5]

total_size = gridsize**3
with open("grid.dat","w") as f:
    f.write("# "+repr(gridsize)+"\n")
    for x in range(0,gridsize):
        for y in range(0,gridsize):
            for z in range(0,gridsize):
                xp = x/float(gridsize)-center[0]
                yp = y/float(gridsize)-center[1]
                zp = z/float(gridsize)-center[2]
                r = np.sqrt(xp**2+yp**2+zp**2)
                if(x==gridsize/2 and y==gridsize/2 and z==gridsize/2):
                    mydensity = density
                    my_dust_density = dust_density
                else:
                    mydensity = 0.0
                    my_dust_density = 0.0
                f.write(repr(mydensity)+" "+repr(temperature)+" "+repr(velocity)+" "+repr(velocity)+" "+repr(velocity)+" "+repr(my_dust_density)+"\n")
                
