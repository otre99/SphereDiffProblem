#from .core import sphere_distance, gen_lonlat_u, gen_lonlat_v, gen_lonlat
from .calc import calc_grid_coords, calc_grid_coords_uv

import numpy as np 
from math import pi
import matplotlib.pyplot as plt  
from math import sin, cos
import sys 
from scipy.interpolate import RectSphereBivariateSpline

class WFUNCT:
    def __init__(self, NP, D, SP):
        self.NP   = NP    
        self.DATA = D
        self.SP   = SP 
        self.J, self.I  = D.shape 

    def __call__(self, j, i):
        if j==-1:
            return self.NP
        elif j == self.J:
            return self.SP
        elif i==-1:
            return self.DATA[j,self.I-1]
        elif i==self.I:
            return self.DATA[j,0]
        else:
            return self.DATA[j,i]

    def get_flatten(self):
        return np.hstack( [ [self.NP], self.DATA.flatten(), [self.SP] ] )

def get_Q(u, v, R):
    J, I = u.shape
    #_, lats = gen_lonlat(nlats=J, nlons=I)
    lats, _ = calc_grid_coords(nlats=J, nlons=I,return_grid=False,geo_latlon=False)
    
    dlat =   pi / (J+1)
    dlon = 2*pi / I
        
    Q_NP = -R*dlon*v[0,:].sum()
    Q    = np.empty(u.shape)
    for j in range(J):
        for i in range(I-1):
            Q[j, i] = -R*( dlon*( v[j+1,i  ]*sin(lats[j]+0.5*dlat) -v[j,i   ]*sin(lats[j]-0.5*dlat)  ) + dlat*(u[j, i+1]-u[j, i  ]) )
        Q[j,I-1]    = -R*( dlon*( v[j+1,I-1]*sin(lats[j]+0.5*dlat) -v[j, I-1]*sin(lats[j]-0.5*dlat)  ) + dlat*(u[j,   0]-u[j, I-1]) )
    Q_SP = -R*dlon*v[-1,:].sum()

    return WFUNCT(Q_NP, Q, Q_SP)

def get_L_B_T_C(I, J):
    lats, _ = calc_grid_coords(nlats=J, nlons=I,return_grid=False, geo_latlon=False)
    #_, lats = gen_lonlat(nlats=J, nlons=I)
    dlat =   pi / (J+1)
    dlon = 2*pi / I

    L = dlat/(dlon*np.sin(lats))
    B = (dlon/dlat)*np.sin( lats+0.5*dlat )
    T = (dlon/dlat)*np.sin( lats-0.5*dlat )
    C = (dlon/dlat)*(np.sin(lats+0.5*dlat) + np.sin(lats-0.5*dlat) ) + (2*dlat)/(dlon*np.sin(lats))
    return L, B, T, C


def wind_to_advdiff_grid(lons, lats, u, v, I, J, u_poles=[None, None], v_poles=[None, None], method="rbf"):

    u_lats, u_lons, v_lats, v_lons = calc_grid_coords_uv(nlats=J, nlons=I, return_grid=True)

    if method == "rbf":

        from scipy.interpolate import Rbf
        from scipy.spatial.distance import cdist
        
        flat_lons = lons.flatten()
        flat_lats = lats.flatten()
        flat_u    = u.flatten() 
        flat_v    = v.flatten()

    
        def sphere_distance(p1, p2, R=1):
            from math import sin, cos, acos    
            lon1, lat1 = p1
            lon2, lat2 = p2 
            return R*acos( sin(lat1)*sin(lat2)*cos(lon1-lon2)+cos(lat1)*cos(lat2))

        def norm(p1, p2):
            return cdist(p1.reshape(-1,2), p2.reshape(-1,2), metric=sphere_distance)

        u_rbf = Rbf(flat_lons, flat_lats, flat_u, function='linear', norm=norm)
        v_rbf = Rbf(flat_lons, flat_lats, flat_v, function='linear', norm=norm)

        u = u_rbf(u_lons,u_lats)
        v = v_rbf(v_lons,v_lats)

        return u, v 
    elif method == "spline":
        from scipy.interpolate import RectSphereBivariateSpline

        u_sp = RectSphereBivariateSpline(lats, lons, u, pole_values=u_poles)
        v_sp = RectSphereBivariateSpline(lats, lons, v, pole_values=v_poles)

        u = u_sp.ev(u_lats, u_lons)
        v = v_sp.ev(v_lats, v_lons)
        return u, v 

def get_uv_from_scalar_potential(north, data, south, R=1.0):
    J, I = data.shape
    lats, _ = calc_grid_coords(nlats=J, nlons=I, return_grid=False,geo_latlon=False)
    dlat =   pi / (J+1)
    dlon = 2*pi / I

    u = np.empty(shape=(J  , I))
    v = np.empty(shape=(J+1, I))

    # U 
    u[:,0] = ( data[:,0]-data[:,I-1] )/( dlon*np.sin(lats)*R )                
    for i in range(1,I):
        u[:,i] = ( data[:,i]-data[:,i-1])/( dlon*np.sin(lats)*R )        

    v[0,:]    =  (data[0,:]-north)/(R*dlat)
    v[1:-1,:] =  (data[1:,:]-data[0:-1,:])/(R*dlat)
    v[-1,:]    =  (south - data[-1,:])/(R*dlat)

    return u, v 

def get_sparse_matrix(I, J, L, B, T, C):
    from scipy.sparse import csr_matrix
    N = I*J+2
    indices = []
    indptr  = []
    data    = []

    # fisrt row 
    row = 0
    indptr.append( len(indices) )
    indptr.append( len(indices) + I + 1  )     
    indices += np.arange(I+1).tolist()
    data    += np.ones(shape=(I+1)).tolist(); data[0] = -I    
    row =+1

    # second row
    for i in range(I):
        indptr.append( len(indices) + 5 )
        if i == 0:
            indices += [0   ,    row, row+1, row+I-1, row+I]
            data    += [T[0],  -C[0],  L[0],    L[0],  B[0]]
        elif i == I-1:
            indices += [0   ,    1, row-1,  row,  row+I]
            data    += [T[0], L[0],  L[0], -C[0],   B[0]]
        else:
            indices += [0   , row-1,  row,  row+1, row+I]
            data    += [T[0],  L[0], -C[0],   L[0],  B[0]]
        row+=1 
    
    # midlw rows 
    for j in range(1, J-1):
        for i in range(I):
            indptr.append( len(indices) + 5 )
            if i == 0:
                indices += [row-I,    row,  row+1, row+I-1, row+I]
                data    += [T[j] ,  -C[j],   L[j],    L[j],  B[j]]
            elif i == I-1:
                indices += [row-I, row-I+1, row-1,  row,  row+I]
                data    += [ T[j],    L[j],  L[j], -C[j],   B[j]]
            else:
                indices += [ row-I, row-1, row,   row+1, row+I ]
                data    += [ T[j],  L[j], -C[j],   L[j],  B[j] ]
            row+=1 

    # row before end 
    for i in range(I):
        indptr.append( len(indices) + 5  )
        if i == 0:
            indices += [row-I ,     row,   row+1,   row+I-1,      N-1]
            data    += [T[J-1],  -C[J-1],  L[J-1],    L[J-1],   B[J-1]]
        elif i == I-1:
            indices += [row-I  ,   row-I+1,   row-1,    row,       N-1]
            data    += [ T[J-1],    L[J-1],  L[J-1], -C[J-1],   B[J-1]]
        else:
            indices += [row-I,   row-1,     row,   row+1,    N-1 ]
            data    += [T[J-1],  L[J-1], -C[J-1],  L[J-1],  B[J-1]]
        row+=1 


    #last row 
    indptr.append( len(indices) + I + 1  )     
    indices += np.arange(N-I-1,N).tolist()

    data    += (-np.ones(shape=(I+1))).tolist(); data[-1]=I    
    sparse_matrix = csr_matrix( (data, indices, indptr), shape=(N,N) )
    return sparse_matrix

def get_non_divergent_wind_to_delete(u, v, R, ERROR_TH=1e-3):
    from math import sin, cos 

    J, I = u.shape

    PHI = WFUNCT(0.0, np.zeros(u.shape), 0.0) 
    Q   = get_Q(u, v, R)
    L, B, T, C = get_L_B_T_C(I=I, J=J)

    ER=0
    while True:
        TMP = -( Q.NP - PHI.DATA[0,:].sum() ) / I          
        ER = abs( PHI.NP - TMP )
        PHI.NP = TMP

        for j in range(J):
            for i in range(I):
                TMP = -( Q(j,i) - B[j]*PHI(j+1,i)-T[j]*PHI(j-1,i) - L[j]*(PHI(j,i+1)+ PHI(j,i-1) ) )/C[j]
                ER = max(ER, abs( PHI(j,i)-TMP)  )
                PHI.DATA[j,i] = TMP

        TMP = ( Q.SP + PHI.DATA[-1,:].sum() ) / I
        ER = max(ER, abs( PHI.SP-TMP)  )
        PHI.SP = TMP
        
        if ER < ERROR_TH:
            break

    return PHI.get_flatten()

    #out_u, out_v = get_uv_from_scalar_potential( PHI.NP, PHI.DATA, PHI.SP, R)
    #out_u += u 
    #out_v += v 


    #return out_u, out_v

def get_non_divergent_wind(u, v, R, tol=1e-8, atol=1e-8, maxiter=99999, display=False, method="lgmres"):
    import scipy.sparse.linalg as spla
    import scipy.sparse as sp
    import scipy.linalg as la

    J, I = u.shape
    Q   = get_Q(u=u, v=v, R=R)
    L, B, T, C = get_L_B_T_C(I=I,J=J)

    class Display:
        def __init__(self, A, b):
            self.A = A
            self.b = b
            self.last_sol=None
        def __call__(self, res):
            if self.last_sol is None:
                self.last_sol = res.copy()
            else:
                print( "   ||b-Ax||^2", la.norm(self.A.dot(res) - self.b) )
                print( "||x_2-x_1||^2", la.norm(self.last_sol-res))
                print("-"*16)
                self.last_sol = res.copy()            

    sparse_matrix = get_sparse_matrix(I=I, J=J, L=L, B=B, T=T, C=C)
    callback = Display(sparse_matrix, Q.get_flatten()) if display else None

    M2 = spla.spilu( sp.csc_matrix(sparse_matrix) )
    M_x = lambda x: M2.solve(x)
    M = spla.LinearOperator(M2.shape,M_x)

    
    if method=="lgmres":
        PHI = spla.lgmres(sparse_matrix, Q.get_flatten(), maxiter=maxiter, tol=tol, atol=atol, callback=callback)#, M=M)
    if method=="bicgstab":
        PHI = spla.bicgstab(A=sparse_matrix, b=Q.get_flatten(),  tol=tol, atol=atol, callback=callback, maxiter=maxiter)
    #if method=="gmres":
    #    PHI = spla.gmres(A=sparse_matrix, b=Q.get_flatten(),  tol=tol, atol=atol, callback=callback, maxiter=maxiter)
    if method=="minres":
        PHI = spla.minres(A=sparse_matrix, b=Q.get_flatten(),  tol=tol, callback=callback, maxiter=maxiter)
    if method=="bicg":
        PHI = spla.bicg(sparse_matrix, Q.get_flatten(), maxiter=maxiter, tol=tol, atol=atol, callback=callback) #, M=M)
    if method=="qmr":
        PHI = spla.qmr(sparse_matrix, Q.get_flatten(), maxiter=maxiter, tol=tol, atol=atol, callback=callback) #, M=M)
    if method=="gcrotmk":
        PHI = spla.gcrotmk(sparse_matrix, Q.get_flatten(), maxiter=maxiter, tol=tol, atol=atol, callback=callback) #, M=M)
    if method=="gc":
        PHI = spla.cg(sparse_matrix, Q.get_flatten(), maxiter=maxiter, tol=tol, atol=atol, callback=callback) #, M=M)
    if method=="gcs":
        PHI = spla.cgs(sparse_matrix, Q.get_flatten(), maxiter=maxiter, tol=tol, atol=atol, callback=callback) #, M=M)        
    if display:
        if PHI[1] == 0:
            print("Successful")
        elif PHI[1]>0:
            print("Convergence to tolerance not achieved, number of iterations: ", PHI[1])
        else:
            print("Illegal input or breakdown")

    PHI = PHI[0]
    out_u, out_v = get_uv_from_scalar_potential(PHI[0], PHI[1:-1].reshape(u.shape[0], u.shape[1]), PHI[-1], R=R)
    out_u += u 
    out_v += v 
    return out_u, out_v, PHI[1:]

    
    
    



