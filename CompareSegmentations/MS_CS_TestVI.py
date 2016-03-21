#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

# Purpose: compare 2D segmentation stacks using the VI criterion
# as defined in:
# Marina Meila, "Comparing clusterings.an information based distance"
# Journal of Multivariate Analysis 98 (2007) 873 - 895
#

import os, sys, re
import numpy
import math
import h5py
from PIL import Image
from skimage import data, io, filters
import MS_LIB_IO
import MS_LIB_Util

ms_home = os.environ['MS_HOME'] # where source code is located
ms_data = os.environ['MS_DATA'] # dir for initial input data and final output
ms_temp = os.environ['MS_TEMP'] # dir for intermediate data files

# -------------------------------------------------------------------------------

if __name__ == "__main__":

    if len(sys.argv) == 1:
        # Fake data: two images with N identical regions
        # Regions in image 1 are assigned labels 1,...,N
        # Regions in image 2 have the same labels, except the last region, which has label N-1
        # In other words, the last two regions in image 2 aree merged inbto one region
        show_data = 0
        N = 3
        imdata1 = numpy.zeros([N, N])
        imdata2 = numpy.zeros([N, N])
        for n in range(0, N):
            imdata1[n,:] = n + 1
            imdata2[n,:] = n + 1
        imdata2[N-1,:]   = N - 1
        print "imdata1=", imdata1
        print "imdata2=", imdata2
        # Example from Jan Funke's paper on Tolerant Distance:
#       for n in range(0, N):
#           if n < 3:
#               imdata1[n,:] = 1
#           else:
#               imdata1[n,:] = 2
#           if n < 2:
#               imdata2[n,:] = 3          
#           else:
#               imdata2[n,:] = 4
    elif len(sys.argv) in [4,5,6]:
        # Real data from two stacks
        seg_stack1_path = sys.argv[1]
        seg_stack2_path = sys.argv[2]
        layer_id        = int(sys.argv[3])
        transpose2 = 0
        if len(sys.argv) >= 5:
            transpose2 = int(sys.argv[4])
        show_data = 0
        if len(sys.argv) >= 6:
            show_data = int(sys.argv[5])

        imdata1 =     numpy.squeeze(numpy.transpose(MS_LIB_IO.read_input1(seg_stack1_path, layer_id)))
        if transpose2 > 0:
            print "Transposing data2 ..."
            imdata2 = numpy.squeeze(numpy.transpose(MS_LIB_IO.read_input1(seg_stack2_path, layer_id)))
        else:
            print "Do not transpose data2 ..."
            imdata2 = numpy.squeeze(                MS_LIB_IO.read_input1(seg_stack2_path, layer_id))
    else:
        sys.exit("\nusage: MS_CS_TestVI.py [ seg_stack_path1 seg_stack_path2 layer_if [ transpose [ show_data ]]] \n")

    # Visualize the data
    print "imdata1.shape=", imdata1.shape
    print "imdata2.shape=", imdata2.shape

    if show_data:
        imdata10 = imdata1.copy()
        edges1 = filters.sobel(numpy.float64(imdata10))
        imdata10[edges1 > 0] = 0
        imdata10[imdata10 > 0] = 255
        img1 = Image.fromarray(numpy.uint8(imdata10))
        img1.thumbnail(imdata10.shape, Image.ANTIALIAS)
        img1.show()
  
        imdata20 = imdata2.copy() 
        edges2 = filters.sobel(numpy.float64(imdata20))
        imdata20[edges2 > 0] = 0
        imdata20[imdata20 > 0] = 255 
        img2 = Image.fromarray(numpy.uint8(imdata20))
        img2.thumbnail(imdata20.shape, Image.ANTIALIAS)
        img2.show()

    # Compute probabilities
    N1, N2, N12 = MS_LIB_Util.get_VI_counts(imdata1, imdata2)
    print "\nnum_N1_values=",  len(N1.values()),  " N1.values()=",  sorted(N1.values())
    print "\nnum_N2_values=",  len(N2.values()),  " N2.values()=",  sorted(N2.values())
    print "\nnum_N12_values=", len(N12.values()), " N12.values()=", sorted(N12.values())

    VI = MS_LIB_Util.counts2VI2D(N1, N2, N12)

    print "VI= ", VI

