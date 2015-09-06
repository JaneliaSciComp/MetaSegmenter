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

# This script takes a dir of labels files and outputs RGB files

if __name__ == "__main__":

    if len(sys.argv) == 3:
        indir  = sys.argv[1] # dir of segmented labels files
        outdir = sys.argv[2] # dir of files with labeled RGB regions                            
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
    else:
        sys.exit("usage: MS_UT_SegLabels2Binary.py  <dir_segmentation_labels> <dir_BW>  \n")

    for file in os.listdir(indir):
        inpath  = os.path.join(indir, file)
        outpath = os.path.join(outdir, file)
        label_image = numpy.asarray(Image.open(inpath, 'r'))  
        bwdata = label_image.copy()
        bwdata[label_image > 0] = 1

        img = Image.fromarray(numpy.uint8(bwdata/numpy.max(bwdata)*255.))
        img.save(outpath)

