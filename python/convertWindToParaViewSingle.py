#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 15:31:43 2018

@author: rbt
"""

import numpy as np
import sys 
from math import pi 
import argparse
import AdvDiffPyTools.io as IO 
import AdvDiffPyTools.core as CORE 
import AdvDiffPyTools.calc as CALC


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('u',  help="U input file") 
    parser.add_argument('v',  help="V input file") 

    parser.add_argument('--scale_factor', default=1.0, type=float, help='Scale vector size')
    parser.add_argument('--radius', type=float, default=100.0, help='Sphere radius')
    parser.add_argument('--normalize', action="store_true", help='Normalize velocity')
    parser.add_argument('--output_file', help='Output file (*.vtk file)')
    return parser    


def gen_uv(I,J, alpha):
    ulats, ulons, vlats, vlons = calc_grid_coords_uv(nlats=J, nlons=I, return_grid=True)
    
    u =  U0*( np.sin(ulats)*cos(alpha) - np.cos(ulats)*np.cos(ulons)*sin(alpha) )
    v = -U0*np.sin(vlons)*sin(alpha)
    return u, v 

if __name__ == '__main__':

    parser = get_parser()
    FLAGS, unparsed = parser.parse_known_args()
    if len(unparsed):
        print("Warning: unknow arguments {}".format(unparsed))
        sys.exit(2)

    R = FLAGS.radius 

    u = IO.read_bin_output(FLAGS.u)
    v = IO.read_bin_output(FLAGS.v)
    scale_factor = FLAGS.scale_factor

    u, v = u[CORE.KeyData], v[CORE.KeyData]
    nlats, nlons = u.shape
    
    if (nlats+1,nlons) != v.shape:
        print("Error with u and v dimensions")
        sys.exit(1)
    
    u, v = CALC.staggereduv_to_centeruv(U=u, V=v)

    dd={ "u": u*scale_factor, "v": v*scale_factor }    

    output_file = "out.vtk" if FLAGS.output_file is None else FLAGS.output_file
    IO.vtk_save_file_uv(dd=dd, fname=output_file, R=R, with_cell_normals=False)
    
