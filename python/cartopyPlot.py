#!/home/programs/Anaconda/envs/maria/bin/python3
import AdvDiffPyTools.core as CORE 
import AdvDiffPyTools.plot_utils as PLOT_UTILS
import AdvDiffPyTools.io as IO
import matplotlib.pyplot as plt 
import sys 
import argparse
from AdvDiffPyTools.misc import sec_to_strtime
import numpy as np

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file',  help="model output file") 
    parser.add_argument('--cmap_str', default="Blues")
    parser.add_argument('--del_neg', action="store_true", help='Replace neg. values by zero')
    return parser    

if __name__ == "__main__":
    parser = get_parser()
    FLAGS, unparsed = parser.parse_known_args()
    if len(unparsed):
        print("Warning: unknow arguments {}".format(unparsed))
        sys.exit(2)

    data = IO.read_bin_output(FLAGS.input_file, layout="F")
    print("INFO")
    for key in data.keys():
        if key == "Data": continue
        print(key, "-->", data[key])
    vmax = max(data[CORE.KeySouthPole], data[CORE.KeyNorthPole], data[CORE.KeyData].max())
    print("Min. value: ", min(data[CORE.KeySouthPole], data[CORE.KeyNorthPole], data[CORE.KeyData].min())  )
    print("Max. value: ", vmax  )




    #levels = [0] + list(np.linspace(0.0, vmax, 10)[1:])
    levels = np.linspace(0.0, vmax, 10)
    print(levels)
    ax=PLOT_UTILS.plot_scalar(data, cmap_str=FLAGS.cmap_str,
                           title=sec_to_strtime(data["Time"]),
                           del_neg=FLAGS.del_neg,
                           levels=levels, #np.linspace(0.0, data[CORE.KeyData].max(), 10),
                           #levels=[0.0, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.1],
                           #levels=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.05],                           
                           #levels=np.linspace(0, 0.05, 8),                           
                           bg_color="white",                            
                           #mode="pcolormesh"
                           mode="contourf"
                           )
    #plt.gcf().set_size_inches(9,4)
    plt.gcf().set_size_inches(16,7)
    plt.subplots_adjust(top=0.912, bottom=0.086, left=0.017, right=0.983, hspace=0.2, wspace=0.2)
    plt.tight_layout()
    #ax.set_xlim(-75,75)
    #ax.set_ylim(-45,45)
    plt.savefig(FLAGS.input_file+".eps")

    plt.show()
