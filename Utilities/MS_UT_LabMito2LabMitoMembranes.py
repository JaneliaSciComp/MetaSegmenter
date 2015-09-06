#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

from PIL import Image
import numpy
from skimage.morphology import disk, binary_erosion
import sys, os

if __name__ == "__main__":

    if len(sys.argv) == 3:
        if os.path.isdir(sys.argv[1]):
            indir    = sys.argv[1]
            outdir   = sys.argv[2]
            if not os.path.isdir(outdir):
                os.mkdir(outdir)
            infiles  = []
            outfiles = []
            for file in os.listdir(indir):
                infiles.append(os.path.join(indir,file))
                outfiles.append(os.path.join(outdir, file.split(".")[0] + ".png"))
        else:
            sys.exit("Input type is not recognized")
    else:
        sys.exit("usage: MS_UT_LabMito2LabMitoMembranes.py  <dir of mito labels> <dir of mito membrane labels> \n")

    for i in range(0, len(infiles)):
        I  = numpy.asarray(Image.open(infiles[i], 'r'))  
        strel = disk(2)
        Ie = binary_erosion(I, strel)
        Im = numpy.zeros(I.shape, dtype="int")
        Im[:] = I[:]
        Im[Ie > 0] = 0
        img = Image.fromarray(numpy.uint8(Im/numpy.max(Im)*255.))
        img.save(outfiles[i])                 
