from AdvDiffPyTools.io import save_bin_file
from AdvDiffPyTools.core import KeyData, KeyNorthPole, KeySouthPole, KeyTime
from AdvDiffPyTools.calc import calc_grid_coords
import  numpy as np 

def save_all(out_folder, phi=None, ff=None, sigma=None, mulon=None, mulat=None, u=None, v=None, layout="F", 
             custom_names={"phi":"phi", "ff":"ff", "sigma":"sigma", "mulon":"mulon", "mulat":"mulat", "u":"u", "v":"v"} ):
    from os.path import join
    if phi   :
        oname="phi"
        if oname in custom_names.keys(): oname = custom_names[oname] 
        save_bin_file(join(out_folder, oname), phi, layout=layout)        
    if ff    : 
        oname="ff"
        if oname in custom_names.keys(): oname = custom_names[oname] 
        save_bin_file(join(out_folder, oname), ff, layout=layout)
    if sigma : 
        oname="sigma"
        if oname in custom_names.keys(): oname = custom_names[oname] 
        save_bin_file(join(out_folder, oname), sigma, layout=layout)    
    if mulon : 
        oname="mulon"
        if oname in custom_names.keys(): oname = custom_names[oname] 
        save_bin_file(join(out_folder, oname), mulon, layout=layout)    
    if mulat : 
        oname="mulat"
        if oname in custom_names.keys(): oname = custom_names[oname] 
        save_bin_file(join(out_folder, oname), mulat, layout=layout)    
    if u     : 
        oname="u"
        if oname in custom_names.keys(): oname = custom_names[oname] 
        save_bin_file(join(out_folder, oname), u, layout=layout)
    if v     : 
        oname="v"
        if oname in custom_names.keys(): oname = custom_names[oname] 
        save_bin_file(join(out_folder, oname), v, layout=layout)
    
def empty_problem(I, J):
    lats, lons = calc_grid_coords(nlats=J,  nlons=I, return_grid=True)    
    phi = {KeyNorthPole: 0.0, KeyData: np.zeros_like(lats), KeySouthPole: 0.0, KeyTime: 0.0} 
    fff = {KeyNorthPole: 0.0, KeyData: np.zeros_like(lats), KeySouthPole: 0.0, KeyTime: 0.0} 
    sigma = {KeyNorthPole: 0.0, KeyData: np.zeros_like(lats), KeySouthPole: 0.0, KeyTime: 0.0} 

    shape = (J,I)
    vshape = (J+1,I)
    u = np.zeros(shape); u = {KeyData: u}
    v = np.zeros(vshape); v = {KeyData: v}
    
    mu_lon = np.zeros(shape); mu_lon = {KeyData: mu_lon}
    mu_lat = np.zeros(vshape); mu_lat = {KeyData: mu_lat}
        
    return lats, lons, phi, fff, sigma, mu_lon, mu_lat, u, v 

class GenNamelist:
    def __init__(self, GS : bool = False) -> None:
        self.keys=["nlats", "nlons", "n_species", "nsave", "dt", "tf", "method", "standard_grid", "radius", 
                   "phi0_path", "ff0_path", "sigma0_path", "mulon0_path", "mulat0_path", "output_folder", "npar"]
        self.data={}
        self.data["nlats"]="359"
        self.data["nlons"]="720" 
        self.data["n_species"]="1" 
        self.data["nsave"]="1" 
        self.data["dt"]="0.1" 
        self.data["tf"]="1" 
        self.data["method"]="1" 
        self.data["standard_grid"]=".true." 
        self.data["radius"]="1.0", 
        self.data["phi0_path"]="\"\"" 
        self.data["ff0_path"]="\"\"" 
        self.data["sigma0_path"]="\"\"" 
        self.data["mulon0_path"]="\"\"" 
        self.data["mulat0_path"]="\"\"" 
        self.data["output_folder"]="\"\""
        self.data["npar"]="1"
        self.GS = GS
        if self.GS:
            self.data["fparam"]="0.30"
            self.data["kparam"]="0.60"     
            self.keys.append("fparam")
            self.keys.append("kparam")   


    def set(self, key, value):
        if key not in self.keys:
            raise("Unknow key = "+key)
        
        if key in ("phi0_path", "ff0_path", "sigma0_path", "mulon0_path", "mulat0_path", "output_folder"):
            value = "\"{}\"".format(value)
        self.data[key]=value        

    def gen_namelist(self):
        result="&input\n"
        for k in self.keys:
            value = self.data[k]
            result+="{}={}\n".format(k, value)
        result+='/\n'
        return result

def gen_parser():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--I', type=int, default=720, help="nlons") 
    parser.add_argument('--J', type=int, default=359, help='nlats')
    parser.add_argument('--dt', type=float, default=0.005, help='Time step')
    parser.add_argument('--radius', type=float, default=1.0, help='Radius')
    parser.add_argument('--nsave', type=int, default=50, help='nsave')
    parser.add_argument('--mt', type=int, default=3, help='method')
    parser.add_argument('--tf', type=float, default=10.0, help='tf')
    return parser

