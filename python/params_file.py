import numpy as np
from collections import OrderedDict
class params_file(OrderedDict):

    def __init__(self,fname=None):
        self.float_types = ["cosmology.H0","cosmology.Omega_M","cosmology.Omega_L","boxsize","redshift","hubble_flow","tau_max","dust.albedo",
"dust.dust_metal_ratio","slab.column_density","slab.column_density_dust","slab.density","slab.density_dust","max_step",
"emission.minimum_luminosity","ramses.scale_density","ramses.scale_dust_density","shell.column_density","shell.column_density_dust",
"shell.density","shell.density_dust""emission.fixed_width","slab.temperature","shell.inner_radius","shell.outer_radius","shell.temperature","shell.outflow_velocity",]
        self.bool_types = ["no_hubble_flow","use_peeling_off","output.binary","output.do_slices","split_domain"]
        self.int_types = ["number_of_instruments","output.slices.resolution","verbosity","max_num_peeling_off_photons",
 "emission.minimum_number_of_photons_per_source","max_num_peeling_off_photons","number_of_photons","octet.rootlevel",
"octet.distribution_method.row_dir","rng.seed"]
        self.float_vector_types = ["line_of_sight","emission.bounding_box_lo","emission.bounding_box_hi"]
        self.int_vector_types = []
        
        if(fname is not None and isinstance(fname,str)):
            super(params_file,self).__init__()
            self.read(fname)
        else:
            super(params_file,self).__init__(fname)
            
            
    def update(self):
        self.write(self.original_name)
    def read(self,fname):
        self.original_name = fname
        self.clear()
        los_count = 0
        with open(fname,"r") as f:
            lines = f.readlines()
            for line in lines:
                line=line.rstrip()
                line=line.lstrip()
                if(len(line)<1 or line[0] is '#'):
                    continue
                line=line.split('=')
                if(len(line)>1):
                    key=line[0].strip()
                    value = line[1].split("#")[0]
                    value = self.type_conversion(key,value.strip())
                    if("line_of_sight" in key):
                        key += "%"+repr(los_count)
                        los_count+=1
                    self[key]=value
    def write(self,fname,comment=None):
        with open(fname,"w") as f:
            f.write("#file written by params_file class\n")
            if(comment is not None):
                f.write("#"+comment+"\n")
            for key,value in  self.iteritems():
                if("%" in key):
                    key = key.split("%")[0]
                if(isinstance(value,str)):
                    v = value
                elif(isinstance(value,np.ndarray)):
                    v = ""
                    for d in value:
                        v+=repr(d)+" "
                else:
                    v = repr(value)
                if(isinstance(value,bool)):
                    v = v.lower()
                f.write(key+"="+v+"\n")
    def type_conversion(self,key,value):
        if key in self.float_types:
            return float(value)
        if key in self.int_types:
            try:
                v = int(value)
            except ValueError:
                print "Warning:",key,"should be integer, but is",value
                v= int(round(float(value)))
            
            return v
        if key in self.float_vector_types:
            return np.array([float(v) for v in value.split()])
        if key in self.int_vector_types:
            return np.array([int(v) for v in value.split()])
        if key in self.bool_types:
            if(value in ["true","True"]):
                   ret = True
            elif(value in ["false","False"]):
                   ret = False
            else:
                raise NameError("Could not convert parameter to boolean!")
            return ret
        
        return value
    
    
