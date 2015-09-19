#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

import os, sys, numpy
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pyplot import figure, ion, show, draw
ion()

if __name__ == "__main__":

    if len(sys.argv) == 2:            
        image_file = sys.argv[1]
    else:
        sys.exit("\nusage: MS_UT_ShowImage image_file\n")

image_data = numpy.asarray(Image.open(image_file))
if len(image_data.shape) == 2:
    plt.imshow(image_data, cmap = cm.Greys_r)
    plt.show()
else:
    plt.imshow(image_data)
    show()
try:
    id = input("Press enter to remove image")
    i = int(id)
    if i == 0:
        sys.exit("Exiting ...")
except SyntaxError:
    pass
    



