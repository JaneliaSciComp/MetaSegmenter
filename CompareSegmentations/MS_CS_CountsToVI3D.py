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

        # Compute probabilities
    P1   = {}
    P2   = {}
    P12  = {}
    N1s  = 0
    N2s  = 0
    N12s = 0

    print "Computing P1  ..."
    for k1 in N1.keys():
        N1s     += N1[k1]
    for k1 in N1.keys():
        P1[k1]   = float(N1[k1])/float(N1s)
#   print "k1=", k1, " N1=", N1[k1], " P1[k1]=", P1[k1]

    print "Computing P2  ..."
    for k2 in N2.keys():
        N2s     += N2[k2]
    for k2 in N2.keys():
        P2[k2]   = float(N2[k2])/float(N2s)
#   print "k2=", k2, " N2=", N2[k2], " P2[k2]=", P1[k2]

    print "Computing P12  ..."
    for k12 in N12.keys():
        N12s    += N12[k12]
    for k12 in N12.keys():
        P12[k12] = float(N12[k12])/float(N12s)
#   print "k12=", k12, " N12=", N12[k12], " P12[k12]=", P12[k12]
    print "N1s=", N1s, " N2s=", N2s, " N12s=", N12s
    print "len(N1.keys())=", len(N1.keys())
    print "len(N2.keys())=", len(N2.keys())
    print "len(N12.keys())=", len(N12.keys())

    # Compute VI3D
    H1 = 0
    H2 = 0
    I  = 0

   
    print "\nComputing H1  ..." 
    for k in N1.keys():
        if P1[k] > 0:
            H1 += -P1[k]*math.log(P1[k])
        else:
            print "Warning: P1[", k, "]=", P1[k]
    print "H1=", H1

    print "\nComputing H2  ..."
    for k in N2.keys():
        if P2[k] > 0:
            H2 += -P2[k]*math.log(P2[k])
        else:
            print "Warning: P2[", k, "]=", P2[k]
    print "H2=", H2

    print "\nComputing I  ..."
    for k1 in N1.keys():
        for k2 in N2.keys():
            k12 = str(k1) + "_" + str(k2)
            try:                   
                if P1[k1] > 0 and P2[k2] > 0 and P12[k12] > 0:
                    I += P12[k12]*math.log(P12[k12]/P1[k1]/P2[k2])
            except:
                continue
    print "I=", I

    VI3D = H1 + H2 - 2*I

    # Output VI2D to file
    f = open(outpath, 'wt')
    f.write("VI3D = " + str(VI3D) + " for nodes " + str(zmin) + "-" + str(zmax))
    print "VI3D = ", VI3D, " for nodes ", zmin, "-", zmax
