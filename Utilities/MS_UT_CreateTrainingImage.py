#! /usr/local/python-2.7.6/bin/python

import os, re, sys
import numpy
from PIL import Image

if __name__ == "__main__":
    if len(sys.argv) in [4,5]:
        raw_image_file     = sys.argv[1]
        pos_BW_image_file  = sys.argv[2]
        out_RGB_image_file = sys.argv[3]
        neg_BW_image_file  = ""
        if len(sys.argv) == 5:
            neg_BW_image_file = sys.argv[4]
    else:
        sys.exit("\nUsage: MS_UT_CreateTrainingImage.py <raw_file> <pos_BW_file> <out_RGB_file> [ <neg_BW_file> ]\n")

pos_BW_data  = numpy.array(Image.open(pos_BW_image_file))
pos_BW_shape = pos_BW_data.shape
RGB_shape    = (pos_BW_shape[0],pos_BW_shape[1],3)
RGB_data     = numpy.zeros(RGB_shape, 'uint8')
if re.search("membrane", pos_BW_image_file) and len(neg_BW_image_file) == 0:      
    # train membranes vs the rest: membranes=black, the rest =white
    RGB_data[pos_BW_data==  0,1] = 255   # positive examples: black -> green
    RGB_data[pos_BW_data==255,0] = 255   # negatice examples: white -> red
elif len(neg_BW_image_file) == 0:    
    # train mito or mito-membranes vs the rest
    RGB_data[pos_BW_data==255,1] = 255   # positive examples: white -> green
    RGB_data[pos_BW_data==  0,0] = 255   # negatice examples: black -> red
else:
    # train membranes vs mito
    neg_BW_data   = numpy.array(Image.open(neg_BW_image_file))
    RGB_data[pos_BW_data==  0,1] = 255   # positive examples: white -> green
    RGB_data[neg_BW_data==255,0] = 255   # negatice examples: white -> red    
RGB_image = Image.fromarray(RGB_data)
parent_dir = os.path.dirname(os.path.normpath(out_RGB_image_file))
print "out_RGB_image_file=", out_RGB_image_file, " parent_dir=", parent_dir
if not os.path.isdir(parent_dir):
    os.mkdir(parent_dir)
RGB_image.save(out_RGB_image_file)

