#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

from PIL import Image
import numpy
from skimage.measure import label
from skimage.color import label2rgb
from skimage.segmentation import clear_border
import sys, os

if __name__ == "__main__":

    if len(sys.argv) == 3:
        indir  = sys.argv[1] # binary regions (lells + mitochondria)
        outdir = sys.argv[2] # labeled regions                            
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
    else:
        sys.exit("usage: MS_UT_ProduceRGBLabels.py  <dir_binary_regions> <dir_out>  \n")

    for file in os.listdir(indir):
        inpath  = os.path.join(indir, file)
        outpath = os.path.join(outdir, file)
        indata  = numpy.asarray(Image.open(inpath, 'r'))  
        # Clear border
        bwdata = indata.copy()
        b = 2
        bwdata[0:b,:] = 0
        bwdata[-b:,:] = 0
        bwdata[:,0:b] = 0
        bwdata[:,-b:] = 0

        outdata = label(bwdata)

        img = Image.fromarray(numpy.uint8(outdata))
        img.save(outpath)

