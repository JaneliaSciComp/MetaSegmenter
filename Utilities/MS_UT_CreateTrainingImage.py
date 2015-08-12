#! /usr/local/python-2.7.6/bin/python

import os, re, sys
import numpy
from PIL import Image

if __name__ == "__main__":
    if len(sys.argv) == 4:
        raw_image_file  = sys.argv[1]
        BW_image_file   = sys.argv[2]
        RGB_image_file  = sys.argv[3]
    else:
        sys.exit("\nUsage: MS_UT_CreateTrainingImage.py <raw_image_file> <BW_file> <RGB_file>\n")

BW_data = numpy.array(Image.open(BW_image_file))
BW_shape = BW_data.shape
RGB_shape = (BW_shape[0],BW_shape[1],3)
RGB_data = numpy.zeros(RGB_shape, 'uint8')
if re.search("membrane", BW_image_file):      
    # membranes=black, cytoplasm=white
    RGB_data[BW_data==  0,1] = 255   # positive examples: black -> green
    RGB_data[BW_data==255,0] = 255   # negatice examples: white -> red
else:
    RGB_data[BW_data==255,1] = 255   # positive examples: white -> green
    RGB_data[BW_data==  0,0] = 255   # negatice examples: black -> red
RGB_image = Image.fromarray(RGB_data)
parent_dir = os.path.dirname(os.path.normpath(RGB_image_file))
print "RGB_image_file=", RGB_image_file, " parent_dir=", parent_dir
if not os.path.isdir(parent_dir):
    os.mkdir(parent_dir)
RGB_image.save(RGB_image_file)

