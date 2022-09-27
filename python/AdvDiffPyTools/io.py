from math import pi
import numpy as np 

from  .calc import  latlon_to_xyz, calc_grid_coords, wind2d_to_wind3d
from . import core as CORE 

def get_saved_file_list(dir_path):
    import os
    from os.path import join
    fl = os.listdir(dir_path)
    
    def is_ok(w):
        """
        0123456789012345678901234
        sp00_0000000000000497.out
        """        
        if len(w) == 25 and w[:2] == "sp" and w[-4:] =='.out': return True
        return False    
    l = []
    for fname in fl:
        if  is_ok( fname ):
            l.append(join(dir_path, fname))
    l.sort()
    return l 

def vtk_write_points(dd, fout, a=100):   
    nlats, nlons = dd[CORE.KeyData].shape
    dlon = 2*pi/nlons
    dlat = pi/(nlats+1)
    lats = np.linspace(0.5*dlat, pi-0.5*dlat, nlats+1)
    lons = np.linspace(0, 2*pi-dlon, nlons)
    lons, lats = np.meshgrid(lons, lats)
    x, y, z = latlon_to_xyz(lats=lats, lons=lons, a=a)
    fout.write("POINTS {} float\n".format(x.size) )
    for i in range(x.size):
        fout.write("{} {} {}\n".format(x.flat[i], y.flat[i], z.flat[i]))

def vtk_write_cell_polygons(dd, fout, a=100):
        _, nlons = dd[CORE.KeyData].shape
        s = dd[CORE.KeyData]
        fout.write("POLYGONS {} {}\n".format(s.size+2, 5*s.size + 2*(nlons+1)))
        fout.write("{} ".format(nlons))
        for i in range(nlons-1, -1, -1):
            fout.write(" {}".format(i))
        fout.write("\n")
        for i in range(s.size):
            col = i%nlons
            i1 = i
            i2 = i1+1 if col<nlons-1 else i1-col
            i3 = i1+nlons
            i4 = i3+1 if col<nlons-1 else i3-col
            fout.write("4 {} {} {} {}\n".format(i1,i2,i4,i3))        
        fout.write("{} ".format(nlons))
        for i in range(nlons):
            fout.write(" {}".format(s.size+i))        
        fout.write("\n")

def vtk_write_cell_scalars(dd, fout, name):        
    fout.write("SCALARS {} float\n".format(name))
    fout.write("LOOKUP_TABLE default\n")
    s = dd[CORE.KeyData]
    
    fout.write("{}\n".format(dd[CORE.KeyNorthPole]))
    for i in range(s.size):
        fout.write("{}\n".format(s.flat[i]))
    fout.write("{}\n".format(dd[CORE.KeySouthPole]))            

def vtk_write_cell_normals(dd, fout):
        fout.write("NORMALS cell_normals float\n")
        
        la, lo = calc_grid_coords(nlats=dd[CORE.KeyData].shape[0], nlons=dd[CORE.KeyData].shape[1], return_grid=True)
        lo, la = np.meshgrid(lo,la)
        nx, ny, nz = latlon_to_xyz(la, lo, 1.0)
        s = dd[CORE.KeyData]
        fout.write("0.0 0.0 1.0\n")
        for i in range(s.size):
            fout.write("{} {} {}\n".format(nx.flat[i],ny.flat[i],nz.flat[i]))
        fout.write("0.0 0.0 -1.0\n")    

def vtk_write_cell_vectors(dd, fout, R=100):
        fout.write("VECTORS cell_vector float\n")
        lats, lons = calc_grid_coords(nlats = dd['u'].shape[0], nlons = dd['u'].shape[1], return_grid=True)

        _, _, _, vx, vy, vz = wind2d_to_wind3d(lats=lats, lons=lons, u=dd['u'], v=dd['v'], a=R)
        fout.write("0.0 0.0 0.0\n")
        for i in range(vx.size):
            c1, c2, c3 = vx.flat[i],vy.flat[i],vz.flat[i]
            fout.write("{} {} {}\n".format(c1, c2, c3))
        fout.write("0.0 0.0 0.0\n")        

def vtk_save_file_uv(dd, fname, R=100, with_cell_normals=False):
   
    with open(fname, 'wt') as fout:
        fout.write("# vtk DataFile Version 2.0\n")
        fout.write("Sphere Surface\n")
        fout.write("ASCII\n")            
        fout.write("DATASET POLYDATA\n")
        vtk_write_points({CORE.KeyData:dd['u']}, fout, R)
        vtk_write_cell_polygons({CORE.KeyData:dd['u']}, fout, R)

        fout.write("CELL_DATA {}\n".format(dd['u'].size+2))
        vtk_write_cell_vectors(dd, fout, R=R)
        if with_cell_normals:
            vtk_write_cell_normals(dd, fout)

def vtk_save_file(dd, fname, R=100, with_cell_normals=False):   
    s = dd[CORE.KeyData]
    with open(fname, 'wt') as fout:
        fout.write("# vtk DataFile Version 2.0\n")
        fout.write("Sphere Surface\n")
        fout.write("ASCII\n")            
        fout.write("DATASET POLYDATA\n")
        vtk_write_points(dd, fout, R)
        vtk_write_cell_polygons(dd, fout, R)

        fout.write("CELL_DATA {}\n".format(s.size+2))
        vtk_write_cell_scalars(dd, fout, "phi")
        if with_cell_normals:
            vtk_write_cell_normals(dd, fout)

def read_data_from_file(fname, a=1.0):    
    from .core import AdvDiffSolution
    """
    Read the content of a single file saved with program
    """
    f = open(fname)
    result = {} 
    for l in f:
        if ':' in l:
            key, value = l.split(':')
            result[key.strip()] = float(value)
        else:
            nPoleSol, sPoleSol  = l.split()
            nPoleSol, sPoleSol  = float(nPoleSol), float(sPoleSol)                
            
            

            result[CORE.KeyNorthPole] = nPoleSol
            result[CORE.KeySouthPole] = sPoleSol    
            break
    result[CORE.KeyData] = np.loadtxt(f)    

    return AdvDiffSolution(n_sol=result[CORE.KeyNorthPole], 
                           sol=result[CORE.KeyData], 
                           s_sol=result[CORE.KeySouthPole],
                           a=a,
                           time=result[CORE.KeyTime])
    
def save_uv_matrices(fname, u, v, time=0.0):
    if u.shape[0] != v.shape[0]-1 or u.shape[1] != v.shape[1]:
        print("Error with matrix shapes. u.shape = {} and v.shape = {}".format(u.shape, v.shape))
        return 
    with open(fname, "wt") as ofile:
        ofile.write("Time: {}\n".format(time))
        np.savetxt(ofile, u)
        np.savetxt(ofile, v)

def save_matrix(fname, data, time=0.0):
    with open(fname, "wt") as ofile:
        ofile.write("Time: {}\n".format(time))
        np.savetxt(ofile, data)


def save_initial_cond(folder, phi, f=0.0,  uv=0.0, mu=0.0, time=None):
    from .core import AdvDiffSolution, AdvDiffWind
    from os.path import join 
    if time is None:
        time = 0.0 

    phi.save(join(folder, "PHI_0000.TXT"), time=time)


    if isinstance(f, float):
        J, I = phi.sol.shape 
        ff = np.empty(shape=(J,I))
        ff[:,:] = f
        f = AdvDiffSolution(n_sol=f, sol=ff, s_sol=f, time=time)
    f.save(join(folder, "SOURCES_0000.TXT"), time=time)

    if isinstance(uv, float):
        J, I = phi.sol.shape 
        u, v = np.empty(shape=(J,I)), np.empty(shape=(J+1,I))        
        u[:,:] = uv 
        v[:,:] = uv 
        uv = AdvDiffWind(u=u, v=v)
    uv.save(join(folder, "WIND_0000.TXT"))

    if isinstance(mu, float):
        J, I = phi.sol.shape 
        m1, m2 = np.empty(shape=(J,I)), np.empty(shape=(J+1,I))
        m1[:,:] = mu 
        m2[:,:] = mu 
        mu = (m1, m2)
    save_uv_matrices(join(folder, "MU_0000.TXT"), u=mu[0], v=mu[1], time=time)



def read_bin_output(fname, layout="F"):
    from .core import KeyTime, KeyData, KeySouthPole, KeyNorthPole

    with open(fname, "rb") as ifile:
        data = np.fromfile(ifile, dtype=np.float64)
        tt, rows, cols, N = data[0], int(data[1]), int(data[2]), int(data[3])
    r={}
    if N-2 != rows*cols:
        r[KeyTime] = tt
        if layout == "F":
            r[KeyData] =   data[4:].reshape(cols, rows).T          
        else:
            r[KeyData] =   data[4:].reshape(rows, cols)                      
        return r
    
    r[KeyTime] = tt
    r[KeyNorthPole] = data[4]
    if layout=="F":
        r[KeyData] = data[5:-1].reshape(cols,rows).T 
    else:
        r[KeyData] = data[5:-1].reshape(rows,cols) 
    r[KeySouthPole] = data[-1] 

    return r

def save_bin_file(fname, r, layout="F"):
    from .core import KeyTime, KeyData, KeySouthPole, KeyNorthPole
    J, I = r[KeyData].shape    
    T = 0.0
    if KeyTime in r.keys(): T = r[KeyTime]   
  
    tmp = r[KeyData]
    if layout == "F": tmp = r[KeyData].T

    if KeySouthPole in r.keys():
        raw_data = np.concatenate( [  [ r[KeyNorthPole] ], tmp.ravel(), [ r[KeySouthPole] ]  ] )
    else:
        raw_data = tmp.ravel()
 
    header=[T, J, I, raw_data.size]
    raw_data = np.concatenate([header, raw_data]).astype(np.float64)
    raw_data.tofile(fname)


def read_bin_output_limited(fname):
    from .core import KeyTime, KeyData

    with open(fname, "rb") as ifile:
        data = np.fromfile(ifile, dtype=np.float64)
        tt, rows, cols = data[0], int(data[1]), int(data[2])

    r={}
    r[KeyTime] = tt
    r[KeyData] =   data[3:].reshape(rows, cols)  
    return r

def get_info(output_folder, a=1.0):
    r={}

    t=[]
    mass=[]
    l2norm=[]
    vmin=[]
    vmax=[]
    s_pole=[]
    n_pole=[]
    flist = get_saved_file_list(output_folder)
    for fname in flist:
        dd = read_data_from_file(fname,a=a)
        mass.append( dd.mass() )
        l2norm.append(dd.l2_norm() )
        s_pole.append( dd.s_sol )
        n_pole.append( dd.n_sol )
        vmin.append( dd.sol.min() )
        vmax.append( dd.sol.max() )
        t.append( dd.time )
    r["t"] = np.array(t)
    r["mass"] = np.array(mass)
    r["l2norm"] = np.array(l2norm)
    r["s_pole"] = np.array(s_pole)
    r["n_pole"] = np.array(n_pole)
    r["vmin"] = np.array(vmin)
    r["vmax"] = np.array(vmax)
    return r 
