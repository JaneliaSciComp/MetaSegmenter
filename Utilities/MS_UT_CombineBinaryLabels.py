#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

from PIL import Image
import numpy
import sys, os

if __name__ == "__main__":

    if len(sys.argv) == 4:
        indir1 = sys.argv[1] # cellular regions (membranes = black)
        indir2 = sys.argv[2] # labels of mitochondria membranes
        outdir = sys.argv[3] # final regions (lells + mitochondria)
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
    else:
        sys.exit("usage: MS_UT_CombineBinaryLabels.py  <dir_labels1> <dir_labels2> <dir_out>  \n")

    for file in os.listdir(indir1):
        inpath1 = os.path.join(indir1, file)
        inpath2 = os.path.join(indir2, file.split(".")[0] + ".png")
        outpath = os.path.join(outdir, file.split(".")[0] + ".png")
        data1  =              numpy.asarray(Image.open(inpath1, 'r'))  
        data2  = numpy.invert(numpy.asarray(Image.open(inpath2, 'r'))) # mito membranes = now black
        data   = numpy.zeros(data1.shape, dtype="int")
        data[:] = data1[:]
        data[data2 == 0] = 0.
        img = Image.fromarray(numpy.uint8(data/numpy.max(data)*255.))
        img.save(outpath)

