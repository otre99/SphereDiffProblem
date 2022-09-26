import numpy as np 
from math import pi
from math import sin 
from .calc import calc_grid_coords, calc_pole_area, calc_cell_area, calc_mayavi_mesh


KeyTime="Time"
KeyData="Data"
KeyNorthPole="NorthPole"
KeySouthPole="SouthPole"
KeyVMin="vmin"
KeyVMax="vmax"
KeyCutLat="cutlat"
KeyCutLon="cutlon"
KeyLonGrid="longrid"
KeyLatGrid="latgrid"

class AdvDiffSolution:
    def __init__(self, sol=None, n_sol=None, s_sol=None, a=1.0, time=0.0, nlats=None, nlons=None):       
        
        if callable(sol):
            self.I = nlons
            self.J = nlats           
            lats, lons = self.grid(return_grid=True, geo_latlon=False)
            self.sol, self.n_sol, self.s_sol = sol(lats, lons)
        elif isinstance(sol, np.ndarray):            
            self.sol   = sol.copy() 
            self.n_sol = n_sol
            self.s_sol = s_sol 
            self.J, self.I = self.sol.shape                         
        else:
            self.I = nlons
            self.J = nlats           
            self.sol = np.zeros(shape=(nlats, nlons))
            self.n_sol = 0.0
            self.s_sol = 0.0

        self.time = time        
        self.dlon = 2*pi/self.I
        self.dlat =   pi/(self.J+1)
        self.a = a

    def cell_value(self, jlat, ilon):
        return self.sol[jlat, ilon]
        
    def north_pole_value(self):
        return self.n_sol
        
    def south_pole_value(self):
        return self.s_sol

    def grid(self, return_grid=True, geo_latlon=False):
        return calc_grid_coords(nlats=self.J, nlons=self.I, return_grid=return_grid, geo_latlon=geo_latlon)
    
    def l2_norm(self):
        tmp = calc_pole_area(self.dlat, self.a)*self.n_sol**2
        lats, _ = calc_grid_coords(nlats=self.J, nlons=1, return_grid=False)
        for j, lat in enumerate(lats):
            area = calc_cell_area(lat=lat, dlat=self.dlat, dlon=self.dlon, a=self.a)
            tmp += area*(self.sol[j,:]**2).sum()
        tmp += calc_pole_area(self.dlat, self.a)*self.s_sol**2
        return tmp 

    def mass(self):
        tmp = calc_pole_area(self.dlat, self.a)*self.n_sol
        lats, _ = calc_grid_coords(nlats=self.J, nlons=1)

        for j, lat in enumerate(lats):
            area = calc_cell_area(lat=lat, dlat=self.dlat, dlon=self.dlon, a=self.a)
            tmp += area*self.sol[j,:].sum()
        tmp += calc_pole_area(self.dlat, self.a)*self.s_sol
        return tmp         

    def mayavi_mesh(self, a=None):
        if a is None:
            a = self.a 
        return calc_mayavi_mesh(n_sol=self.n_sol, sol = self.sol, s_sol = self.s_sol,  a=a)

    def save_to_vtk(self, fname, a=None):
        from .io import vtk_save_file
        if a is None:
            a = self.a 
        dd = {}
        dd["NorthPole"] = self.n_sol
        dd["Data"] = self.sol
        dd["SouthPole"] = self.s_sol
        vtk_save_file(dd, fname, R=a)

    def save(self, fname, time=None):
        with open(fname, "wt") as ofile:
            if time is None:
                time = self.time
            ofile.write("Time: {}\n".format(time))
            ofile.write("{} {}\n".format(self.n_sol, self.s_sol))
            np.savetxt(ofile, self.sol)

class AdvDiffWind:
    def __init__(self, uv=None, u=None, v=None, a=1.0, time=0.0, nlats=None, nlons=None):       
        
        if callable(uv):
            self.I = nlons
            self.J = nlats
            coords = self.uv_grid(return_grid=True)            
            self.u, self.v = uv(*coords)
        else:            
            self.u = u.copy() 
            self.v = v.copy()
            self.J, self.I = self.u                                     

        self.time = time        
        self.dlon = 2*pi/self.I
        self.dlat =   pi/(self.J+1)
        self.a = a

    def to_center(self):               
        u = np.empty(self.u.shape)    
        u[:,:-1] = 0.5*( self.u[:,:-1] + self.u[:,1:] )
        u[:,-1]  = 0.5*( self.u[:,-1]  + self.u[:,0] )

        v = 0.5*(self.v[:-1,:] + self.v[1:,:])
        return u, v

    def uv_grid(self, return_grid=False):
        from .calc import calc_grid_coords_uv
        return calc_grid_coords_uv(nlats=self.J, nlons=self.I,return_grid=return_grid)        

    def save_to_vtk(self, fname, a=None):
        from .io import vtk_save_file_uv
        if a is None:
            a = self.a 
        dd = {}
        dd["u"], dd["v"] = self.to_center()
        vtk_save_file_uv(dd, fname, R=a)

    def divergence(self):
        from .calc import get_divergence
        d = get_divergence(self.u, self.v, a=self.a)
        return d["Data"], d["NorthPole"], d["SouthPole"]


    def save(self, fname, time=None):
        from .io import save_uv_matrices
        if time is None:
            time = self.time
        save_uv_matrices(fname,u=self.u, v=self.v, time=time)

