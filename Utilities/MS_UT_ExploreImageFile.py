#! /usr/local/python-2.7.6/bin/python

import re, sys
import h5py
from scipy import misc
import tifffile as tiff

if __name__ == "__main__":
    if len(sys.argv) == 2:
        image_file = sys.argv[1]
    else:
        sys.exit("\nusage: show_image.py image_file\n")

if re.search(".png", image_file) or re.search(".jpeg", image_file):
    image_data = misc.imread(image_file)
elif re.search(".tif", image_file):
    image_data = tiff.imread(image_file)
elif re.search(".h5", image_file):
    f  = h5py.File(image_file, 'r')
    key = f.keys()[0]
    image_data = f[key]
else:
    sys.exit("Unsupported file format")
print "Image shape=", image_data.shape

