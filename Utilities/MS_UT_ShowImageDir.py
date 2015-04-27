#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

import os, sys, numpy
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pyplot import figure, ion, show
ion()

if __name__ == "__main__":

    if len(sys.argv) in [2, 3]:
        dir_path = sys.argv[1]
        file_id = 1
        if len(sys.argv) > 2:
            file_id = int(sys.argv[2]) -1
            if  file_id <= 0:
                sys.exit("Exiting ...")
    else:
        sys.exit("\nusage: MS_UT_ShowImageDir.py dir_path [ first_file# (>= 1)]\n")

files = sorted(os.listdir(dir_path))     
num_files = len(files) 
print "num_files=", num_files               
i = file_id-1
while i < num_files:
    print "File#=", i+1, " name=", files[i]
    myfile = os.path.join(dir_path, files[i])
    image_data = numpy.asarray(Image.open(myfile))
    if len(image_data.shape) == 2:
        plt.imshow(image_data, cmap = cm.Greys_r)
        show()
    else:
        plt.imshow(image_data)
        plt.show()
    try:
        id = input("Enter layer_id (=0 to exit) or press enter to continue to next layer")
        print "id=", id
        i = int(id)
        if i == 0:
            sys.exit("Exiting ...")
    except SyntaxError:
        i = i + 1
        pass
    



