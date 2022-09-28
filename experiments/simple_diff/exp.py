from audioop import mul
import sys     
sys.path.append("../../python/")

from AdvDiffPyTools.core import KeyData, KeyNorthPole, KeySouthPole
import  numpy as np 
from math import pi, radians   
from  experiments import empty_problem, save_all, GenNamelist


if __name__ == "__main__":
    import os 
    np.random.seed(231)

    I, J = 720,359
    mt = 1  
    dt = 5e-2  
    grid = True 
    npar=6
    nsave=100
    radius=1

    w = 15.0  
    lats, lons, phi, ff, _, mulon, mulat = empty_problem(I=I,J=J)    
    mu = 0.0005
    mulon[KeyData][...]=mu 
    mulat[KeyData][...]=mu 

    phi[KeyNorthPole] = 100.0
    save_all(out_folder=".", phi=phi, mulat=mulat, mulon=mulon, ff=ff, 
             custom_names={"phi":"phi", "mulon":"mulon",  "mulat":"mulat", "ff":"ff"})

    name_lst = GenNamelist(GS=False)
    name_lst.set("nlats", J)
    name_lst.set("nlons", I)
    name_lst.set("n_species", 1)
    name_lst.set("nsave", nsave)
    name_lst.set("dt", dt)
    name_lst.set("tf", 100000.0)
    name_lst.set("method", mt)
    name_lst.set("standard_grid", grid)
    name_lst.set("radius", radius)
    name_lst.set("phi0_path", "../experiments/simple_diff/phi")
    name_lst.set("mulon0_path", "../experiments/simple_diff/mulon") 
    name_lst.set("mulat0_path", "../experiments/simple_diff/mulat")   
    name_lst.set("ff0_path", "../experiments/simple_diff/ff")
    name_lst.set("output_folder", "../experiments/simple_diff/")
    name_lst.set("npar", npar)
    with open("namelist", "w") as ofile:
        ofile.write(name_lst.gen_namelist())

    os.system("cp namelist ../../build/")
