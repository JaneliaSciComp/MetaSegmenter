#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

# Purpose: compare 2D segmentation stacks using the VI criterion
# as defined in:
# Marina Meila, "Comparing clusterings.an information based distance"
# Journal of Multivariate Analysis 98 (2007) 873 - 895
#

import os, sys, re, h5py
import json
import numpy
import math
import MS_LIB_Util

ms_home = os.environ['MS_HOME'] # where source code is located
ms_data = os.environ['MS_DATA'] # dir for initial input data and final output
ms_temp = os.environ['MS_TEMP'] # dir for intermediate data files

# -------------------------------------------------------------------------------

if __name__ == "__main__":

    print "len(sys.argv)=", len(sys.argv)
    if len(sys.argv) == 6:
        label1  =     sys.argv[1]
        label2  =     sys.argv[2]
        zmin    = int(sys.argv[3])
        zmax    = int(sys.argv[4])
        outpath =     sys.argv[5]
    else:
        sys.exit("\nusage: MS_CS_ComputeVI2D.py label1 label2 zmin zmax \n")

    z  = zmin
    VI = numpy.zeros([int(zmax) - int(zmin)+1, 4], dtype=numpy.float32)

    # Extract counts
    while z < zmax:                
        print "\n...current layer=", z

        input_path1  = os.path.join(ms_temp, "VI_" + label1 + "_" + label2 + "_N1_counts_z"  + str(z) + ".json")
        with open(input_path1, 'r') as fp:
            N1 = json.load(fp)

        input_path2  = os.path.join(ms_temp, "VI_" + label1 + "_" + label2 + "_N2_counts_z"  + str(z) + ".json")
        with open(input_path2, 'r') as fp:
            N2 = json.load(fp)

        input_path12 = os.path.join(ms_temp, "VI_" + label1 + "_" + label2 + "_N12_counts_z"  + str(z) + ".json")
        with open(input_path12, 'r') as fp:
            N12 = json.load(fp)

        print "z=", z, " num_N1_values=",  len(N1.values()),  " N1.values()=",  sorted(N1.values(),  reverse=True)
        print "z=", z, " num_N2_values=",  len(N2.values()),  " N2.values()=",  sorted(N2.values(),  reverse=True)
        print "z=", z, " num_N12_values=", len(N12.values()), " N12.values()=", sorted(N12.values(), reverse=True)

        curr_VI,VI_U,VI_O = MS_LIB_Util.counts2VI2D(N1, N2, N12)

        VI[z - zmin, 0] = z
        VI[z - zmin, 1] = curr_VI
        VI[z - zmin, 2] = VI_U
        VI[z - zmin, 3] = VI_O
        print "VI, VI_U, VI_O for layer ", z, "= ", VI[z - zmin, 1], VI[z - zmin, 2],  VI[z - zmin, 3]
        z = z + 1
       
    # Output VI to file
    print "Output VI to file:", outpath
    numpy.savetxt(outpath, VI, fmt='%10.5f', delimiter='\t')

