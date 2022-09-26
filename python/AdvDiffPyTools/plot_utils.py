import numpy as np 
from math import pi 
from . import core as CORE
from . import calc as CALC 


try:
    import matplotlib.pyplot as plt 
    import matplotlib.colors as mpl_colors
    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18)
    plt.rc('font', size=18)
    plt.rc('figure', autolayout=True, figsize=(9,6))
    plt.rc('axes', titlesize=16, labelsize=18, linewidth=1.0)
    plt.rc('lines', linewidth=2, markersize=6)
    plt.rc('legend',fontsize=16)
    plt.rc('mathtext', fontset='stix')
    plt.rc('font', family='STIXGeneral')
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False    
    
    
except ImportError:
    mplt_available = False
else:
    mplt_available = True
    
try:
    import cartopy.crs as ccrs
    import cartopy
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
except ImportError:
    cartopy_available = False
else:
    cartopy_available = True
    
   
def mayavi_draw_grid(nlats, nlons, a=100.0, resLat=128, resLon=128):
    from mayavi import mlab 
   
    lons = np.arange(0, 2*pi, 2*pi/nlons)
    dlat = pi/(nlats+1)
    lats = np.linspace(dlat/2, pi-dlat/2, nlats+1)
    
    sslats = np.linspace(lats[0], lats[-1], resLat)
    sslons = np.linspace(0, 2*pi, resLon)
    
    for ll in lons:
        x = a*np.sin(sslats)*np.cos(ll)
        y = a*np.sin(sslats)*np.sin(ll)
        z = a*np.cos(sslats)
        mlab.plot3d(x,y,z, color=(0,0,0), line_width=a*0.01, 
                    tube_radius=0.1)

    for ll in lats:
        x = a*np.sin(ll)*np.cos(sslons)
        y = a*np.sin(ll)*np.sin(sslons)
        z = a*np.cos(ll)
        zz = x.copy()
        zz[:] = z
        mlab.plot3d(x,y,zz, color=(0,0,0), line_width=a*0.01, tube_radius=0.1)

def add_ticks(ax, xticks, yticks):
    crs=ccrs.PlateCarree()    
    ax.set_xticks(xticks, crs=crs)
    ax.set_yticks(yticks, crs=crs)
    lon_formatter = LongitudeFormatter(number_format='.1f', degree_symbol='', dateline_direction_label=True)
    lat_formatter = LatitudeFormatter(number_format='.1f', degree_symbol='')

    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)    

def plot_streamplot(u, v, title="", cb_label="", cmap_str="autumn_r", density=2):
    J, I = u.shape
    if u.shape != v.shape:        
        u, v = CALC.staggereduv_to_centeruv(U=u, V=v)            
    J, I = u.shape
    lats, lons = CALC.calc_grid_coords(nlats=J, nlons=I, return_grid=False, geo_latlon=True)
    
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    ax.set_global()
    ax.coastlines(linewidth=0.5)
    magnitude = np.sqrt(u**2+v**2)
    cb = ax.streamplot(lons, lats, u, v, transform=ccrs.PlateCarree(), 
    linewidth=1.2, density=density, arrowsize=0.7, color=magnitude,
    cmap=plt.cm.get_cmap(cmap_str))
    plt.colorbar(cb.lines, label=cb_label)
    ax.set_title(title)
    add_ticks(ax, [-180, -120, -60, 0, 60, 120, 180], [-78.5, -60, -25.5, 25.5, 60, 80])
    return ax 


def plot_quiver(u, v, title="", cb_label="", rows=30, cols=60):    
    if u.shape != v.shape: u, v = CALC.staggereduv_to_centeruv(U=u, V=v)            
    J, I = u.shape
    lats, lons = CALC.calc_grid_coords(nlats=J, nlons=I, return_grid=False, geo_latlon=False)

    u, v = CALC.wind_interpolate_spline(lons=lons, lats=lats, u=u, v=v, I=cols, J=rows)

    J, I = u.shape
    lats, lons = CALC.calc_grid_coords(nlats=J, nlons=I, return_grid=False, geo_latlon=True)
    
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    ax.set_global()
    ax.coastlines(linewidth=0.6)
    ax.quiver(lons, lats, u, v, transform=ccrs.PlateCarree(),headwidth=2, width=0.002)
    ax.set_title(title)
    add_ticks(ax, [-180, -120, -60, 0, 60, 120, 180], [-78.5, -60, -25.5, 25.5, 60, 80])
    return ax 

def plot_scalar(dd, title="", cb_label="", 
                cmap_str="bwr", 
                mode="contourf", 
                del_neg=False, 
                levels=None,
                bg_color=None):

    import matplotlib.pyplot as plt 
    import cartopy.crs as ccrs
    import cartopy
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


    if type(dd) is dict:
        data = dd[CORE.KeyData]
        lats, lons = CALC.calc_grid_coords(nlons=data.shape[1], nlats=data.shape[0], return_grid=True, geo_latlon=True)
    else:
        data = dd
        lats, lons = CALC.calc_grid_coords(nlons=data.shape[1], nlats=data.shape[0], return_grid=True, geo_latlon=True)


    if del_neg: data[data<0.0] = 0.0

    lons[lons<0]+=360 
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    ax.set_global()
    ax.coastlines(linewidth=0.3,resolution='110m',)
    if mode=="contourf":
        cmap = plt.cm.get_cmap(cmap_str)

        if levels is not None:                        
            colors = []                            
            max_val = max(levels)
            if bg_color:
                colors = [bg_color]                
                for i in range(1,len(levels)-1): colors.append(  cmap( 0.5*(levels[i]+levels[i])/max_val ) )
            else:
                for i in range(0,len(levels)-1): colors.append(  cmap( 0.5*(levels[i]+levels[i]) )/max_val )

            cb = ax.contourf(lons, lats, data, transform=ccrs.PlateCarree(), colors=colors, levels=levels)
        else:
            cb = ax.contourf(lons, lats, data, transform=ccrs.PlateCarree(), cmap = cmap, levels=levels)#, norm=mpl_colors.LogNorm(data.min(), data.max()))

    elif mode == "contour":
        cb = ax.contour(lons, lats, data, transform=ccrs.PlateCarree(), 
                         cmap=plt.cm.get_cmap(cmap_str))
    elif mode == "pcolormesh":
        cb = ax.pcolormesh(lons, lats, data, transform=ccrs.PlateCarree(), cmap=plt.cm.get_cmap(cmap_str))#, norm=mpl_colors.LogNorm(data.min(), data.max()))               
    else:
        print("Error with this parameter: mode={}".format(mode))
        sys.exit(-1)



    #cb = plt.colorbar(cb, label=cb_label, pad=0.05, orientation="horizontal", shrink=1.0)    
    cb = plt.colorbar(cb, label=cb_label, pad=0.01, orientation="vertical")    
    #add_ticks(ax, [-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180], [-78.5, -60, -25.5, 25.5, 60, 80])
    add_ticks(ax, [-180, -90, 0, 90, 180], [-80.0, -60, -20, 20, 60, 80])
    ax.set_title(title, fontsize=28)
    return ax 



    
    
