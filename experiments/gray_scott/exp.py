import sys     
sys.path.append("../../python/")

from AdvDiffPyTools.core import KeyData, KeyNorthPole, KeySouthPole
from AdvDiffPyTools.calc import calc_grid_coords
import  numpy as np 
from math import pi, radians   
from  experiments import empty_problem, save_all, GenNamelist

def add_randon_noise(U, V, mask, ER):
    from numpy.random import uniform 
    s = U[KeyData][~mask].size; U[KeyData][~mask] = 1.0  + uniform(-ER, ER, s)
    s = U[KeyData][mask].size;  U[KeyData][mask]  = 0.50 + uniform(-ER, ER, s)
    
    s = V[KeyData][~mask].size; V[KeyData][~mask] = 0.0  + uniform(0.0, ER, s)       
    s = V[KeyData][mask].size;  V[KeyData][mask]  = 0.25 + uniform(-ER, ER, s)    
    return U, V

def pattern(I, J, U, V):
    from numpy.random import uniform 
    lats, lons, = calc_grid_coords(nlats=J, nlons=I, return_grid=True, geo_latlon=False)    
    DA = radians(5.0)
    ER = 0.05
    mask = (lats>0.5*pi-DA)&(lats<0.5*pi+DA)&(lons>pi-DA)&(lons<pi+DA)   
    add_randon_noise(U,V,mask,ER)

    U[KeyNorthPole] = 1.0 + uniform(-ER, ER, 1)[0]
    U[KeySouthPole] = 1.0 + uniform(-ER, ER, 1)[0]         
    V[KeyNorthPole] = 0.0 + uniform(0.0, ER, 1)[0]
    V[KeySouthPole] = 0.0 + uniform(0.0, ER, 1)[0]    


if __name__ == "__main__":
    import os 
    np.random.seed(231)

    I, J = 1440,719
   
    fparam = 0.046
    kparam = 0.0594 
    mt = 2 # 
    dt = 5.0 # 
    grid = 1 # 
    grid = False 
    npar=6

    lats, lons, phiU, _, _, mulonU, mulatU, = empty_problem(I=I,J=J)    
    _,       _, phiV, _, _, mulonV, mulatV, = empty_problem(I=I,J=J)    

    pattern(I=I, J=J, U=phiU, V=phiV)
      
    mu=0.16
    mulonU[KeyData][...] = mu
    mulatU[KeyData][...] = mu

    mu=0.08
    mulonV[KeyData][...] = mu
    mulatV[KeyData][...] = mu

    save_all(out_folder=".", phi=phiU, mulat=mulatU, mulon=mulonU, 
             custom_names={"phi":"phiU", "mulon":"mulonU",  "mulat":"mulatU"})
    save_all(out_folder=".", phi=phiV, mulat=mulatV, mulon=mulonV, 
             custom_names={"phi":"phiV", "mulon":"mulonV",  "mulat":"mulatV"})



    name_lst = GenNamelist(GS=True)
    name_lst.set("nlats", J)
    name_lst.set("nlons", I)
    name_lst.set("n_species", 2)
    name_lst.set("nsave", 1000)
    name_lst.set("dt", dt)
    name_lst.set("tf", 100000.0)
    name_lst.set("method", mt)
    name_lst.set("standard_grid", grid)
    name_lst.set("radius", 230.0)
    name_lst.set("phi0_path", "../experiments/gray_scott/phi")
    name_lst.set("mulon0_path", "../experiments/gray_scott/mulon") 
    name_lst.set("mulat0_path", "../experiments/gray_scott/mulat")   
    name_lst.set("output_folder", "../experiments/gray_scott/")
    name_lst.set("npar", npar)
    name_lst.set("fparam", fparam)
    name_lst.set("kparam", kparam)
    with open("namelist", "w") as ofile:
        ofile.write(name_lst.gen_namelist())

    os.system("cp namelist ../../build/")
