#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

from PIL import Image
from skimage.morphology import disk, binary_erosion
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
        sys.exit("usage: MS_UT_CombineBinaryLabels.py  <dir_memb_labels> <dir_mito_labels> <dir_comb_labels>  \n")

    for file in os.listdir(indir1):
        inpath1 = os.path.join(indir1, file)
        inpath2 = os.path.join(indir2, file)
        outpath = os.path.join(outdir, file.split(".")[0] + ".png")
        
        data_memb = numpy.asarray(Image.open(inpath1, 'r')) # membranes = black, the rest = white
        data_mito = numpy.asarray(Image.open(inpath2, 'r')) # mitochondria = white, the rest = black
        strel     = disk(3)
        data_mito_1 = numpy.zeros(data_mito.shape,  dtype=numpy.uint8)
        data_mito_1[data_mito > 0] = 1    # having 1's instead of 255's is necessary for applying the binary_erosion operation
        data_mito_e = binary_erosion(data_mito_1, strel)
        data_mito_memb = data_mito_1.copy()
        data_mito_memb[data_mito_e > 0] = 0
        print "(data_mito_1 ==   1).sum()=", (data_mito_1 ==   1).sum()
        print "(data_mito_e ==   1).sum()=", (data_mito_e ==   1).sum()
        data      = numpy.zeros(data_memb.shape, dtype=numpy.uint8)
        data[data_memb == 0]      = 1 # cellular membranes
        data[data_mito_e > 0]     = 2 # mitochondria-interior
        data[data_mito_memb > 0 ] = 3 # mito-membranes
        img = Image.fromarray(numpy.uint8(data))
        img.save(outpath)

