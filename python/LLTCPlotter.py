import sys
from params_file import params_file
import warnings
try:
  from pymses import RamsesOutput
except ImportError:
  warnings.warn('pymses not imported')
import numpy as np
import subprocess
import os
import tempfile
import astropy.units as u
#needs to refer to the path where plotter.exe is located
ABSPATH = os.path.dirname(os.path.realpath(__file__))
ABSPATH = ABSPATH[0:ABSPATH.rfind('python')]
if(not os.path.exists(ABSPATH+'/plotter.exe')):
  warnings.warn('plotter.exe not found at '+ABSPATH)

class LLTCPlotter(object):
    def __init__(self,finputs,params,dry=False):
        self.check_sufficient(params)
        self.finputs = finputs
        self.finputs_params = params_file(finputs)
        #if the input is using a ramses dataset, we do not know the boxlength from the inputs file, 
        #but we will need it to get the extent of the plottted region.
        if("dataset_type" in self.finputs_params and "Ramses" in self.finputs_params["dataset_type"]):
            s = self.finputs_params['ramses.root']
            p, f = os.path.split(s)
            path,_ = os.path.split(p)
            snapnum = int(f.split("_")[-1].split(".")[0])
            ro = RamsesOutput(path,snapnum)
            self.finputs_params["boxsize"] = ro.info["unit_length"].val*100.0 #cm  
            #TODO: if we want to plot things in terms of angular extent, we need to figure out the redshift as well.
        params["los"] = params["los"]
        self.params = params
        self.dir= tempfile.mkdtemp()
        self.params["output_prefix"] = os.path.join(self.dir,"python_")
        self.data = {}
        if(dry):
            args = self.construct_args()
            string = ""
            for arg in args:
                string +=arg+" "
            print string
        else:
            self.run()
    def construct_args(self):
        args = []
        args.append(os.path.join(ABSPATH,"plotter.exe"))
        args.append(self.finputs)
        for key in self.params:
            v = self.params[key]
            if(isinstance(v,int) or isinstance(v,float)):
                value = repr(v)
            elif(isinstance(v,str)):
                 value = v
            else:
                value = ""
                for w in v:
                    value+=repr(w)+" "
            if(key is "operators"):
                for operator in self.params["operators"]:
                    args.append("plotter.operator="+operator)
            if(key is "override"):
                for parameter in self.params["override"]:
                    args.append(parameter+"="+self.params["override"][parameter])
            else:
                args.append("plotter."+key+"="+value)
        return args        
            
    def check_sufficient(self,params):
        required_params = ["los","center","up","width","npixels","oversampling","depth","operators"]
        for p in required_params:
            if(p not in params):
                raise NameError("You need to specify "+p+" in the parameters!")
    def run(self):
        args = self.construct_args()
        #print args
        p = subprocess.Popen(args,shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdout,stderr = p.communicate()
        #print "stdout",stdout
        if(stderr is not ""):
            print "stdout",stdout
            print "stderr",stderr
        #collect the data.
        for operator in self.params["operators"]:
            if("opticalDepthGrid" in operator):
                grid = []
                for i in range(0,self.params[operator+".nlambda"]):
                    path = self.params["output_prefix"]+operator+"_"+repr(i)
                    #CAUTION need to transpose here
                    d = np.loadtxt(path).T
                    grid.append(d)
                d = grid
            elif("slice" in operator):
                d = {}
                for var in ["density","temperature","velocity_x","velocity_y","velocity_z","dust_density","dx"]:
                    path = self.params["output_prefix"]+operator+"_"+var
                    #CAUTION need to transpose here
                    d[var] = np.loadtxt(path).T   
            else:
                path = self.params["output_prefix"]+operator
                #CAUTION need to transpose here
                d = np.loadtxt(path).T
            self.data[operator] = d
    def get_extent(self,centered=True):
        boxsize = self.finputs_params["boxsize"]*u.cm
        imgsize = (self.params["width"]*boxsize.to(u.kpc)).value
        if(centered):
            ext = [-imgsize/2.0,+imgsize/2.0,-imgsize/2.0,+imgsize/2.0]
        else:
            ext = [0,imgsize,0.0,+imgsize]
        return ext
    def store(self,fname):
        import pickle
        with open(fname,"w") as f:
            pickle.dump(self,f)
    def get(key):
        return self.data[key]
    @staticmethod
    def from_file(fname):
        import pickle
        return pickle.load(open(fname))
        
            
        
        
