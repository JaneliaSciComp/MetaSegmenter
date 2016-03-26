#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

# Purpose: compare 3D segmentation stacks using the VI criterion
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

    if len(sys.argv) == 6:
        label1  =     sys.argv[1]
        label2  =     sys.argv[2]
        zmin    = int(sys.argv[3])
        zmax    = int(sys.argv[4])
        outpath =     sys.argv[5]
    else:
        sys.exit("\nusage: MS_CS_ComputeVI2D.py label1 label2 zmin zmax \n")

    z   = zmin
    N1  = {}
    N2  = {}
    N12 = {}

    # Extract counts
    while z < zmax:                
        print "\n...current layer=", z

        input_path1  = os.path.join(ms_temp, "VI_" + label1 + "_" + label2 + "_N1_counts_z"  + str(z+1) + ".json")
        with open(input_path1, 'r') as fp:
            N1_z = json.load(fp)

        input_path2  = os.path.join(ms_temp, "VI_" + label1 + "_" + label2 + "_N2_counts_z"  + str(z+1) + ".json")
        with open(input_path2, 'r') as fp:
            N2_z = json.load(fp)

        input_path12 = os.path.join(ms_temp, "VI_" + label1 + "_" + label2 + "_N12_counts_z"  + str(z+1) + ".json")
        with open(input_path12, 'r') as fp:
            N12_z = json.load(fp)

#       print "N1_z.values=", N1_z.values()
#       print "N2_z.values=", N2_z.values()

        print "Computing N1  ..."
        for k1 in N1_z.keys():
            if k1 in N1.keys():
                N1[k1] += N1_z[k1]
            else:
                N1[k1]  = N1_z[k1]

        print "Computing N2  ..."
        for k2 in N2_z.keys():
            if k2 in N2.keys():
                N2[k2] += N2_z[k2]
            else:
                N2[k2]  = N2_z[k2]

        print "Computing N12 ..." 
        for k12 in N12_z.keys():
            try:
                N12[k12] += N12_z[k12]
            except:
                N12[k12]  = N12_z[k12]
        z = z + 1

    VI3D, VI3D_U, VI3D_O = MS_LIB_Util.counts2VI3D(N1, N2, N12)

    # Output VI2D to file
    f = open(outpath, 'wt')
    f.write("VI3D = " + str(VI3D) + " VI3D_U=", VI3D_U, \
            " VI3D_O=", VI3D_O, " for nodes " + str(zmin) + "-" + str(zmax))
    print "VI3D = ", VI3D, " VI3D_U=", VI3D_U, \
          " VI3D_O=", VI3D_O, " for nodes ", zmin, "-", zmax
