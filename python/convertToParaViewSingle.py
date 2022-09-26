#!/home/programs/Anaconda/bin/python3
import numpy as np
import sys 
from math import pi 
import argparse
import AdvDiffPyTools.core as CORE
import AdvDiffPyTools.io as IO
from os.path import split, splitext, join


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', type=str, help="Input file") 
    parser.add_argument('output_file', type=str, help="Output file") 
    parser.add_argument('--folder', action="store_true", help='Convert all files in folder')
    parser.add_argument('--radius', type=float, default=100.0, help='Radius of the sphere')
    return parser    


def convert_in_folder(input_folder, output_folder, R):
    lst_files = IO.get_saved_file_list(input_folder)
    for fname in lst_files:
        _, base_name = split(fname) 
        base_name, _ = splitext(base_name)                       
        oname = join(output_folder, base_name+".vtk")
    
        dd = IO.read_bin_output(fname)     
        IO.vtk_save_file(dd,oname, R=R)
        print(dd[CORE.KeyTime], fname)

def main():
    parser = get_parser()
    FLAGS, unparsed = parser.parse_known_args()
    if len(unparsed):
        print("Warning: unknow arguments {}".format(unparsed))
        sys.exit(2)
    R = FLAGS.radius
    if FLAGS.folder:
        convert_in_folder(FLAGS.input_file, FLAGS.output_file, R)
    else:
        dd = IO.read_bin_output(FLAGS.input_file)     
        IO.vtk_save_file(dd,FLAGS.output_file, R=R)

    
    
if __name__ == '__main__':
    main()
    
    #r, c = 20, 30
    #dd = {CORE.KeyNorthPole: 100, 
    #      CORE.KeyData: np.zeros((r,c)),
    #      CORE.KeySouthPole: 0}

    #R=101.
    #IO.vtk_save_file(dd, "GRIDR_{}x{}_{}.vtk".format(r,c,R),R=R)
    
    
