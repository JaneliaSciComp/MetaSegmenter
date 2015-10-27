#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

from PIL import Image
from skimage.morphology import disk, binary_erosion
import numpy
import sys, os

if __name__ == "__main__":

    if len(sys.argv) == 3:
        indir  = sys.argv[1] # target = white, background = black    
        outdir = sys.argv[2] # target = black, background = white   
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
    else:
        sys.exit("usage: MS_UT_InvertLabels.py  <dir_labels> <dir_inv_labels>   \n")

    for infile in os.listdir(indir):
        inpath  = os.path.join(indir, infile)
        suffix  = infile.split(".")[len(infile.split("."))-1]
        outfile = infile[:(len(infile) - len(suffix) -1)] + ".png"
        outpath = os.path.join(outdir, outfile)                               
        data_orig = numpy.asarray(Image.open(inpath, 'r')) 
        data_inv  = numpy.zeros(data_orig.shape, dtype=numpy.uint8)
        data_inv  = 255 - data_orig  
        img = Image.fromarray(numpy.transpose(numpy.uint8(data_inv)))
        img.save(outpath)

