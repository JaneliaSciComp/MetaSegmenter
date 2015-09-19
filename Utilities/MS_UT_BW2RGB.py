#! /usr/local/python-2.7.6/bin/python

import re, sys
import h5py
from scipy import misc
import numpy
from PIL import Image

if __name__ == "__main__":
    if len(sys.argv) == 3:
        BW_image_file  = sys.argv[1]
        RGB_image_file = sys.argv[2]
    else:
        sys.exit("\nUsage: MS_UT_BW2RGB.py <BW_file_name> <RGB_file_name>\n")

BW_data = numpy.array(Image.open(BW_image_file))
BW_shape = BW_data.shape
RGB_shape = (BW_shape[0],BW_shape[1],3)
RGB_data = numpy.zeros(RGB_shape, 'uint8')

#RGB_data[BW_data==0,:]   = 255  # inverse black->white, white-> black

#RGB_data[BW_data==0,0]   = 255  # make black -> red

#RGB_data[BW_data==0,1]   = 255  # make black -> green

#RGB_data[BW_data==255,1] = 255  # make white -> red

RGB_data[BW_data==255,1] = 255  # make white -> green

RGB_image = Image.fromarray(RGB_data)
RGB_image.save(RGB_image_file)

