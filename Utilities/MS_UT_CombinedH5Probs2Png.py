#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

from PIL import Image
import h5py
import numpy
import sys, os

if __name__ == "__main__":

    if len(sys.argv) == 4:   
        h5_dir  = sys.argv[1]
        h5_dir2 = sys.argv[2]
        png_dir = sys.argv[3]
        if not os.path.isdir(png_dir):
            os.mkdir(png_dir)
        infiles  = []
        infiles2 = []
        outfiles = []
        for file in os.listdir(h5_dir):
            infiles.append(os.path.join( h5_dir,  file))
            infiles2.append(os.path.join(h5_dir2, file))
            outfiles.append(os.path.join(png_dir, file.split(".")[0] + ".png"))
    else:
        sys.exit("usage: MS_UT_CombinedH5Probs2Png.py <h5_dir> <h5_dir2> <png_dir> \n")

    for i in range(0, len(infiles)):
        f  = h5py.File(infiles[i],  'r')   
        f2 = h5py.File(infiles2[i], 'r')
        data1 = numpy.array(f["probabilities"])
        data2 = numpy.array(f2["probabilities"])
        data  = data1
        data[data2 > data1] = data2[data2 > data1]
        img  = Image.fromarray(numpy.uint8(data/numpy.max(data)*255.))
        img.save(outfiles[i])                 
