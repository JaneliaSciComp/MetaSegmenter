#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

import os, sys, numpy
from PIL import Image

if __name__ == "__main__":

    if len(sys.argv) == 2:
        file_path = sys.argv[1]
    else:
        sys.exit("\nusage: MS_UT_ShowImage.py file_path \n")

size = 1250, 1250
img = Image.open(file_path)
img.thumbnail(size, Image.ANTIALIAS)
img.show()

