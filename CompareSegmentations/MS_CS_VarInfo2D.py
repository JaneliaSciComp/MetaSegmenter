#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

import os, sys, re, h5py
import tifffile as tiff
import numpy
from skimage.morphology import label
from scipy.ndimage import measurements,morphology
import matplotlib.pyplot as plt
from PIL import Image
import math

def read_input(stack_file_name):
    if re.search(".h5", stack_file_name):
        f = h5py.File(stack_file_name, 'r')
        key = f.keys()[0]
        data = numpy.uint64(numpy.transpose(f[key]))
    elif re.search(".tif", stack_file_name):
        data = numpy.uint64(numpy.transpose(tiff.imread(stack_file_name)))
    else:
        sys.exit("Unrecognized stack type")
    return data

# -------------------------------------------------------------------------------

if __name__ == "__main__":

    if len(sys.argv) == 3:
        stack_file_name1 = sys.argv[1]
        stack_file_name2 = sys.argv[2]
    else:
        sys.exit("\nusage: MS_CS_CompareSegmentations3D.py seg_stack1 seg_stack2 \n")

    print "stack_file_name1=", stack_file_name1
    print "stack_file_name2=", stack_file_name2
    data1 = read_input(stack_file_name1)
    print "data1_shape =", data1.shape, " num_layers=", data1.shape[2], " dtype=", data1.dtype

    data2 = read_input(stack_file_name2)
    print "data2_shape =", data2.shape, " num_layers=", data2.shape[2], " dtype=", data2.dtype
    if not data1.shape == data2.shape:
            sys.exit("\nShapes of the two stacks are different")

    layer_id = 1
    num_layers = data1.shape[2]
    VI = []

    # Extract counts
    while layer_id <= num_layers:
        print "\n...current layer=", layer_id
        im1 = data1[:,:,layer_id-1]
        im2 = data2[:,:,layer_id-1]
        max_label1 = int(numpy.max(im1))
        max_label2 = int(numpy.max(im2))
        N1  = {}
        N2  = {}
        N12 = {}
        labels1 = []
        labels2 = []
        print "   ... handling labels1"
        for k in range(1, max_label1+1):
            count = (im1 == k).sum()
            if count > 0:
                labels1.append(k)
        print "   ... handling labels2"
        for k in range(1, max_label2+1):
            count = (im2 == k).sum()
            if count > 0:
                labels2.append(k)
        print "   ... handling N1"
        for k in range(0, len(labels1)):
            k1 = labels1[k]
            count = (im1 == k1).sum()
            if not k1 in N1.keys():
                N1[k1]  = count
            else:
                N1[k1] += count
        print "   ... handling N2"
        for k in range(0, len(labels2)):
            k2 = labels2[k]
            count = (im2 == k2).sum()
            if not k2 in N2.keys():
                N2[k2]  = count
            else:
                N2[k2] += count
        print "   ... handling N12"
        for i in range(0, len(labels1)):
            k1 = labels1[i]
            for j in range(0, len(labels2)):
                k2 = labels2[j]
                k12 = str(k1) + "_" + str(k2)
                count = ((im1 == k1) & (im2 == k2)).sum()
                if not k12 in N12.keys():
                    N12[k12]  = count
                else:
                    N12[k12] += count
        print "Layer", layer_id, " len(labels1)=", len(labels1), " len(keys1)=", len(N1.keys()), " len(labels2)=", len(labels2), " len(keys2)=", len(N2.keys())
#       print "N1.values=", N1.values()
#       print "N2.values=", N2.values()

        # Compute probabilities
        P1   = {}
        P2   = {}
        P12  = {}
        N1s  = 0
        N2s  = 0
        N12s = 0
        for k1 in N1.keys():
            N1s     += N1[k1]
        for k1 in N1.keys():
            P1[k1]   = float(N1[k1])/float(N1s)
#       print "k1=", k1, " N1=", N1[k1], " P1[k1]=", P1[k1]
        for k2 in N2.keys():
            N2s     += N2[k2]
        for k2 in N2.keys():
            P2[k2]   = float(N2[k2])/float(N2s)
#       print "k2=", k2, " N2=", N2[k2], " P2[k2]=", P1[k2]
        for k12 in N12.keys():
            N12s    += N12[k12]
        for k12 in N12.keys():
            P12[k12] = float(N12[k12])/float(N12s)
#       print "k12=", k12, " N12=", N12[k12], " P12[k12]=", P12[k12]
        print "N1s=", N1s, " N2s=", N2s, " N12s=", N12s
        print "len(N1.keys())=", len(N1.keys())
        print "len(N2.keys())=", len(N2.keys())
        print "len(N12.keys())=", len(N12.keys())


        # Compute VI
        H1 = 0
        H2 = 0
        I  = 0
        for k in N1.keys():
            if P1[k] > 0:
                H1 += P1[k]*math.log(P1[k])
            else:
                print "Warning: P1[", k, "]=", P1[k]
        for k in N2.keys():
            if P2[k] > 0:
                H2 += P2[k]*math.log(P2[k])
            else:
                print "Warning: P2[", k, "]=", P2[k]
        for k1 in N1.keys():
            for k2 in N2.keys():
                k12 = str(k1) + "_" + str(k2)
                if k12 in P12.keys():
                    if P1[k1] > 0 and P2[k2] > 0 and P12[k12] > 0:
                        I += P12[k12]*math.log(P12[k12]/P1[k1]/P2[k2])
        print "VI for layer ", layer_id, "= ", H1 + H2 - 2*I
        VI.append(H1 + H2 - 2*I)
        layer_id = layer_id + 1

print "\nVI=", VI



