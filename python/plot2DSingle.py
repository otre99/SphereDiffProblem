#!/home/programs/Anaconda/bin/python3
import numpy as np 
import matplotlib.pyplot as plt 
from AdvDiffPyTools import core as CORE
from AdvDiffPyTools import io as IO 
from AdvDiffPyTools import calc as CALC
import argparse

def draw_poles(parser):
    FLAGS, _ = parser.parse_known_args()
    dd = IO.read_bin_output(FLAGS.input_file)
    data = dd[CORE.KeyData]
    J,I = data.shape
    lats, lons = CALC.calc_grid_coords(nlats=J, nlons=I, return_grid=True, geo_latlon=False)
    
    lats, lons = np.degrees(lats), np.degrees(lons)
    N = 150
    sinlats = np.sin(lats[:N,0]); sinlats = np.concatenate([ [0.0], sinlats ])

    r     = np.linspace(0,1.0,sinlats.size)
    theta = lons[0,:]   

    values = np.empty( (r.size, theta.size) )
    values[0,:]   = dd[CORE.KeyNorthPole]
    values[1:,:] = data[:N,:]
    
    #-- Plot... ------------------------------------------------
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    cb = ax.contourf(theta, r, values, cmap=plt.cm.hot)    
    plt.colorbar(cb)
    plt.grid()
    plt.show()


def get_info(fname, FLAGS):
    dd = IO.read_bin_output(fname)
    data = dd[CORE.KeyData]
    J,I = data.shape
    lats, lons = CALC.calc_grid_coords(nlats=J, nlons=I, return_grid=True, geo_latlon=False)   
    lats, lons = np.degrees(lats), np.degrees(lons)

    info={}
    info[CORE.KeyNorthPole]=dd[CORE.KeyNorthPole]
    info[CORE.KeySouthPole]=dd[CORE.KeySouthPole]
    info[CORE.KeyVMin]=data.min()        
    info[CORE.KeyVMax]=data.max()        

    info["mass"] = CALC.calc_mass(dd,a=FLAGS.radius)
    info["l2"] = CALC.calc_l2(dd,a=FLAGS.radius)
    print("mass", info["mass"])
    print("l2", info["l2"])
    print("Min", info[CORE.KeyVMin])
    print("Max", info[CORE.KeyVMax])
    
    # LAT CUT 
    latAng = FLAGS.cut_lat
    jlat = np.argmin( np.abs( lats[:,0]-latAng) )
    cut_lat = lats[jlat,0]
    ylats = data[jlat,:]
    info[CORE.KeyCutLat] = (cut_lat, lons[0,:], ylats)       



    # LON CUT 
    lonAng = FLAGS.cut_lon
    ilon = np.argmin( np.abs(lons[0,:]-lonAng) )
    cut_lon = lons[0,ilon]        
    if ilon >= I//2: ilon-=(I//2)
            

    n_p = [info[CORE.KeyNorthPole]]        
    r1 = data[:,ilon]        
    s_p = [info[CORE.KeySouthPole]]        
    r2 = data[-1::-1,ilon+I//2]
    r = np.concatenate([n_p, r1, s_p, r2])         
    x = np.arange(r.size)
    info[CORE.KeyCutLon] = (0, J+1, cut_lon, r)       

    info[CORE.KeyData] = data
    info[CORE.KeyLatGrid] = lats    
    info[CORE.KeyLonGrid] = lons
    info[CORE.KeyTime] = dd[CORE.KeyTime]
    return info    
    
def draw_default(FLAGS):  
    FLAGS, _ = parser.parse_known_args()
    info = get_info(FLAGS.input_file[0], FLAGS)

    plt.subplot(2,1,1)       
    X,Y,Z = info[CORE.KeyLonGrid], info[CORE.KeyLatGrid], info[CORE.KeyData]
    plt.pcolormesh(X,Y,Z, cmap=plt.cm.hot)
    #plt.contourf(X,Y,Z, cmap=plt.cm.jet)
    plt.colorbar()
    plt.xlabel("lons")
    plt.ylabel("lats")
    plt.gca().invert_yaxis()
    

    # LAT CUT 
    plt.subplot(2,2,3)        
    ang, x, y = info[CORE.KeyCutLat]
    plt.plot(x,y)
    plt.xlabel("lons")
    plt.ylabel(r"$\phi$")
    plt.title("lat = {:.4f}".format(ang))

    # LON CUT 
    plt.subplot(2,2,4)        
    j0, j1, ang, y = info[CORE.KeyCutLon]
    x = np.arange(len(y))
    plt.plot(x, y)
    plt.xticks([j0, j1], ["np", "sp"])
    plt.title("lon = {:.4f}".format(ang))


    # INFO 
    npole, spole, vmin, vmax = info[CORE.KeyNorthPole], info[CORE.KeySouthPole], info[CORE.KeyVMin], info[CORE.KeyVMax]
    ttime = info[CORE.KeyTime]
    mass, l2 = info["mass"], info["l2"]
    plt.suptitle("Time = {:0.3f}\nVmin = {:04f} Vmax = {:04f} Np = {:04f} Sp = {:04f}\nMass = {:04f} L2 = {:04f}".format(ttime, vmin, vmax, npole, spole, mass, l2))        
    plt.tight_layout()
    plt.show()

def draw_diff(parser):  
    info1 = get_info(FLAGS.input_file[0], FLAGS)
    info2 = get_info(FLAGS.input_file[1], FLAGS)
    

    l2 = CALC.calc_l2(sol=info1, a=1.0)
    er = CALC.sol_diff_l2norm(sol1=info1, sol2=info2, a=1.0)
    print("Error", np.sqrt(er/l2)*100.0)


    #plt.subplot(2,1,1)       
    #X,Y,Z1 = info1[CORE.KeyLonGrid], info1[CORE.KeyLatGrid], info1[CORE.KeyData]
    #_,_,Z2 = info2[CORE.KeyLonGrid], info2[CORE.KeyLatGrid], info2[CORE.KeyData]

    #Z = (Z1-Z2)
    ##Z = 100.0*(Z1-Z2)/( 0.5*(np.abs(Z1)+np.abs(Z2)) )
    #plt.pcolormesh(X,Y,Z, cmap=plt.cm.jet)
    ##plt.contourf(X,Y,Z, cmap=plt.cm.jet)
    #plt.colorbar()
    #plt.xlabel("lons")
    #plt.ylabel("lats")
    #plt.gca().invert_yaxis()

    # LAT CUT 
    plt.subplot(1,2,1) #plt.subplot(2,2,3)        
    ang, x1, y1 = info1[CORE.KeyCutLat]
    ang, x2, y2 = info2[CORE.KeyCutLat]
    
    plt.plot(x1, y1, label="Initial Condition")
    plt.plot(x2, y2, label="Sol. t = {:0.4f}".format(info2[CORE.KeyTime]))    
    #plt.plot(x1,y1, label="file 01")
    #plt.plot(x2,y2, label="file 02")
    
    plt.legend()    
    plt.xlabel("lons")
    plt.ylabel(r"$\phi$")
    plt.title("lat = {:.4f}".format(ang))


    # LON CUT 
    plt.subplot(1,2,2) #plt.subplot(2,2,4)        
    j0, j1, ang, y1 = info1[CORE.KeyCutLon]
    _ ,  _,   _, y2 = info2[CORE.KeyCutLon]

    x = np.arange(len(y1))
    plt.plot(x, y1, label="Initial Condition")
    plt.plot(x, y2, label="Sol. t = {:0.4f}".format(info2[CORE.KeyTime]))
    #plt.plot(x, y1, label="file 01")
    #plt.plot(x, y2, label="file 02")

    plt.legend()    
    plt.xticks([j0, j1], ["np", "sp"])
    plt.title("lon = {:.4f}".format(ang))


    # INFO 
    #npole, spole, vmin, vmax = info[CORE.KeyNorthPole], info[CORE.KeySouthPole], info[CORE.KeyVMin], info[CORE.KeyVMax]
    #ttime = info[CORE.KeyTime]
    #plt.suptitle("Time = {:0.3f}\nVmin = {:04f} Vmax = {:04f} Np = {:04f} Sp = {:04f}".format(ttime, vmin, vmax, npole, spole))        

    fig = plt.gcf()
    fig.set_size_inches(10, 5)
    plt.tight_layout()
    #fig.savefig('{:10.4f}.png'.format(info2[CORE.KeyTime]), dpi=100)
    plt.show()

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', type=str, nargs='+', help="Input file") 
    parser.add_argument('--cut_lat', type=float, default=90.0, help='Plot cut with fixed lat')
    parser.add_argument('--cut_lon', type=float, default=180.0, help='Plot cut with fixed lon')
    parser.add_argument('--radius', type=float, default=1.0, help='Radius')
    FLAGS, _ = parser.parse_known_args()
    
    if len(FLAGS.input_file) == 1:
        draw_default(FLAGS)
    if len(FLAGS.input_file) == 2:
        draw_diff(FLAGS)
        
