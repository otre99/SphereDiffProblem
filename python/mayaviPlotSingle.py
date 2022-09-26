#!/home/programs/Anaconda/envs/mayavi/bin/python3

from mayavi import mlab
import numpy as np
import sys 
import argparse
import AdvDiffPyTools.core as CORE
import AdvDiffPyTools.io as IO
import AdvDiffPyTools.calc as CALC 
from AdvDiffPyTools.plot_utils import mayavi_draw_grid 
from AdvDiffPyTools.misc import sec_to_strtime

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', type=str, help="Input file") 
    parser.add_argument('--radius', type=float, default=1.0, help='Radius of the sphere')
    parser.add_argument('--fixed_scale', action="store_true", help='Use a fixed color scaled')
    parser.add_argument('--show_grid', action="store_true", help='Show grid')
    parser.add_argument('--nlats', type=int, help="N lats") 
    parser.add_argument('--nlons',  type=int,help="N lons") 
    return parser    

def main():
    parser = get_parser()
    FLAGS, unparsed = parser.parse_known_args()
    if len(unparsed):
        print("Warning: unknow arguments {}".format(unparsed))
        sys.exit(2)

    a = FLAGS.radius    
    dd = IO.read_bin_output(FLAGS.input_file)
    t = dd[CORE.KeyTime]
   
    x, y, z, s = CALC.calc_mayavi_mesh(n_sol=dd[CORE.KeyNorthPole], 
                                       sol=dd[CORE.KeyData], 
                                       s_sol=dd[CORE.KeySouthPole], a=a )
   

    #mlab.options.offscreen = True
    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(800, 800))
    if FLAGS.show_grid:
        nlons = FLAGS.nlons if FLAGS.nlons else  dd.I
        nlats = FLAGS.nlats if FLAGS.nlats else  dd.J
        mayavi_draw_grid(nlats=nlats, nlons=nlons,a=a*1.01)

    
    # Make sphere, choose colors
    # Display
    print(x.shape, y.shape, z.shape, s.shape)
    m = mlab.mesh(x, y, z, scalars=s, resolution=32, vmin=s.min(), vmax=s.max(), representation="surface")
    mlab.title("{}".format(sec_to_strtime(t)), height=0.93, size=0.5)
    mlab.colorbar(m, label_fmt="%.2f", orientation="vertical")
    mlab.view(azimuth=0, elevation=45, distance=4.3*a)
    #mlab.savefig("/home/rbt/p.png", (1024,768))
    mlab.show()


if __name__ == '__main__':
    main()
