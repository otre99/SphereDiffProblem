import numpy as np 
from math import sin, cos, exp, pi      

def calc_grid_coords(nlats:int, nlons:int, return_grid=False, geo_latlon=False):
    dlon = 2*pi/nlons
    dlat = pi/(nlats+1)   
    lats = np.linspace(dlat, pi-dlat, nlats)
    lons = np.linspace(0.5*dlon, 2*pi-0.5*dlon, nlons)    
    if return_grid:
        lons, lats = np.meshgrid(lons, lats)

    if geo_latlon:
        lats = np.degrees(lats)
        lons = np.degrees(lons)
        lats = 90.0-lats        
        lons[ lons < 0.0 ] += 360.0
    return lats, lons 

def calc_cell_area(lat, dlat, dlon, a=1.0):
    return np.sin(lat)*dlat*dlon*a**2

def calc_pole_area(dlat, a):
    return (dlat**2)*(a**2)*pi/4

def calc_mayavi_mesh(n_sol, sol, s_sol, a):
    lats, lons = calc_grid_coords(nlats=sol.shape[0], nlons=sol.shape[1], return_grid=False)
    lats = np.concatenate(([0.0], lats ,[pi]))
    lons = np.concatenate((lons, [lons[0]]))    
    lons, lats = np.meshgrid(lons,lats) 

    x, y, z = latlon_to_xyz(lats, lons, a)

    S = np.empty( shape=(lats.shape[0], lats.shape[1]) )
    S[0,:]  = n_sol
    S[-1,:] = s_sol
    S[1:-1,:-1] = sol
    S[1:-1,-1]  = sol[:,0] 
    return x, y, z, S

def latlon_to_xyz(lats, lons, a):
    x = a*np.sin(lats)*np.cos(lons)
    y = a*np.sin(lats)*np.sin(lons)
    z = a*np.cos(lats)
    return x,  y,  z    

def cell_average(clat, clon, a, dlat, dlon, funct):
    from scipy.integrate import dblquad
    f = lambda lat, lon: a**2*funct(lat, lon)
    
    lat1, lat2 = clat-dlat, clat+dlat
    clon1, clon2 = clon-dlon, clon+dlon
    
    r = dblquad(func=f, a=lat1, b=lat2, gfun=lambda x: clon1, hfun=lambda x: clon2)
    return r

def np_cell_average(a, dlat, dlon, funct):
    from scipy.integrate import dblquad
    f = lambda lat, lon: a**2*funct(lat, lon)
    
    lat1, lat2 = 0, dlat
    clon1, clon2 = 0, 2*pi
    
    r = dblquad(func=f, a=lat1, b=lat2, gfun=lambda x: clon1, hfun=lambda x: clon2)
    return r


def sp_cell_average(a, dlat, dlon, funct):
    from scipy.integrate import dblquad
    f = lambda lat, lon: a**2*funct(lat, lon)
    
    lat1, lat2 = pi-dlat, pi
    clon1, clon2 = 0, 2*pi
    
    r = dblquad(func=f, a=lat1, b=lat2, gfun=lambda x: clon1, hfun=lambda x: clon2)
    return r

def staggereduv_to_centeruv(U, V):
    u = np.empty(U.shape)    
    u[:,:-1] = 0.5*( U[:,:-1] + U[:,1:] )
    u[:,-1]  = 0.5*( U[:,-1]  + U[:,0]  )
    v = 0.5*(V[:-1,:] + V[1:,:])
    return u, v

def dt_from_cfl_cond(I, J, a, max_vel, lat_th=None):
    dlat = pi/(J+1)
    dlon = 2*pi/I 
    if lat_th is None:
        dx = min( a*sin(0.5*dlat)*dlon, a*dlat) 
    else:
        dx = min( a*sin(lat_th)*dlon, a*dlat)
    return dx/max_vel
    
def spherical_distance(lats, lons, a, ref_pt):
    rlat, rlon = ref_pt
        
    tmp = sin(rlat)*np.sin(lats)*np.cos(rlon-lons)+cos(rlat)*np.cos(lats) 
    tmp = np.minimum( np.maximum(tmp, -1), 1.0)
    
    sol   = a*np.arccos( tmp ) #sin(rlat)*np.sin(lats)*np.cos(rlon-lons)+cos(rlat)*np.cos(lats))   
    n_sol = a*rlat 
    s_sol = a*(pi-rlat)
    return n_sol, sol, s_sol     

def gauss_spot(lats, lons, a, center_pt, C, F, kexp=2.0):
    sol_n, s, sol_s = spherical_distance(lats=lats, lons=lons, a=a, ref_pt=center_pt)
    s = C*np.exp( -np.abs(F*s**kexp) )
    sol_n = C*exp( -abs(F*sol_n**kexp) )
    sol_s = C*exp( -abs(F*sol_s**kexp) )
    return sol_n, s, sol_s



def get_divergence(u, v, a=1.0):
    from .core import KeyNorthPole, KeySouthPole, KeyData 
    from math import sin, cos   
    J, I = u.shape

    data={}

    dlat = pi / (J+1)
    dlon = 2*pi / I

    lats, _ = calc_grid_coords(nlats=J,nlons=I,return_grid=False)

    data[KeyNorthPole] = (a*dlon*sin(0.5*dlat)*v[0,:]).sum()    
    data[KeyData] = np.empty(shape=(J,I))
    for j in range(J):
        upper_border  = a*sin( lats[j]-0.5*dlat )*dlon
        botton_border = a*sin( lats[j]+0.5*dlat )*dlon
        lr_border = a*dlat
        for i in range(I-1):
            data[KeyData][j,i] = lr_border*( u[j, i+1] - u[j, i  ] ) + botton_border*v[j+1,i  ] - upper_border*v[j,i  ]
        data[KeyData][j,I-1]   = lr_border*( u[j, 0  ] - u[j, I-1] ) + botton_border*v[j+1,I-1] - upper_border*v[j,I-1]
    data[KeySouthPole] = -(a*dlon*sin(0.5*dlat)*v[-1,:]).sum()    
    return data

def wind2d_to_wind3d(lats, lons, u, v, a):
    x, y, z = latlon_to_xyz(lats, lons, a)
    
    vx = -u*np.sin(lons) + v*np.cos(lats)*np.cos(lons)
    vy =  u*np.cos(lons) + v*np.cos(lats)*np.sin(lons)
    vz = -v*np.sin(lats)
       
    return x.flatten(), y.flatten(), z.flatten(), vx.flatten(), vy.flatten(), vz.flatten()


def sol_diff(sol1, sol2, a=1.0):
    from .core import KeyNorthPole, KeySouthPole, KeyData 
    diff_npole = sol1[KeyNorthPole]-sol2[KeyNorthPole]    
    diff_spole = sol1[KeySouthPole]-sol2[KeySouthPole]     
    S = sol1[KeyData] - sol2[KeyData]
    return {KeyNorthPole: diff_npole, KeyData: S, KeySouthPole: diff_spole} 

def get_dlat_dlon(I:int, J:int):
    dlat = pi / (J+1)
    dlon = 2*pi / I
    return dlat, dlon    
        
def sol_diff_l2norm(sol1, sol2, a=1.0):
    from .core import KeyNorthPole, KeySouthPole, KeyData 
    J, I = sol1[KeyData].shape 
    dlat, dlon = get_dlat_dlon(I=I,J=J)
    lats, lons = calc_grid_coords(nlats=J, nlons=I, return_grid=True, geo_latlon=False)
    
    area_pole = calc_pole_area(dlat=dlat,a=a)
    area_non_pole = calc_cell_area(lats, dlat, dlon, a=a)

    D = sol_diff(sol1, sol2, a=a)
    return (area_non_pole*D[KeyData]**2).sum() + (D[KeyNorthPole]**2 + D[KeySouthPole]**2)*area_pole
    
    
def calc_mass(sol, a=1.0):    
    from .core import KeyNorthPole, KeySouthPole, KeyData 
    J, I = sol[KeyData].shape 
    dlat, dlon = get_dlat_dlon(I=I,J=J)
    lats, _ = calc_grid_coords(nlats=J, nlons=I, return_grid=True, geo_latlon=False)
    
    area_pole = calc_pole_area(dlat=dlat,a=a)
    area_non_pole = calc_cell_area(lats, dlat, dlon, a=a)

    return (area_non_pole*sol[KeyData]).sum() + (sol[KeyNorthPole] + sol[KeySouthPole])*area_pole    

def calc_l2(sol, a=1.0):    
    from .core import KeyNorthPole, KeySouthPole, KeyData 
    J, I = sol[KeyData].shape 
    dlat, dlon = get_dlat_dlon(I=I,J=J)
    lats, _ = calc_grid_coords(nlats=J, nlons=I, return_grid=True, geo_latlon=False)
    
    area_pole = calc_pole_area(dlat=dlat,a=a)
    area_non_pole = calc_cell_area(lats, dlat, dlon, a=a)

    return (area_non_pole*sol[KeyData]**2).sum() + (sol[KeyNorthPole]**2 + sol[KeySouthPole]**2)*area_pole   

def get_min_max(sol):
    from .core import KeyNorthPole, KeySouthPole, KeyData 
    rmin = min( sol[KeyNorthPole], sol[KeyData].min(), sol[KeySouthPole]  )
    rmax = max( sol[KeyNorthPole], sol[KeyData].max(), sol[KeySouthPole]  )    
    return rmin, rmax