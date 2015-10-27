#! /usr/local/python-2.7.6/bin/python

import re, sys
import h5py
from scipy import misc
import tifffile as tiff
import numpy

if __name__ == "__main__":
    if len(sys.argv) == 2:
        image_file = sys.argv[1]
    else:
        sys.exit("\nUsage: MS_UT_ExploreImageFile.py <image_file>\n")

if re.search(".png", image_file) or re.search(".jpeg", image_file):
    image_data = misc.imread(image_file)
elif re.search(".tif", image_file):
    image_data = tiff.imread(image_file)
elif re.search(".h5", image_file):
    f  = h5py.File(image_file, 'r')
    key = f.keys()[0]
    image_data = numpy.array(f[key])
else:
    sys.exit("Unsupported file format")
print "Image shape=", image_data.shape
print "dtype=", image_data.dtype
max_value = int(image_data.max())
min_value = int(image_data.min())
print "max=", max_value, " min=", min_value
if round(max_value) == max_value and round(min_value) == min_value and max_value - min_value > 1:
    num_values = 0
    for i in range(min_value, max_value + 1):
        if (image_data == i).sum() > 0:
            num_values += 1
    print "num_values=", num_values

