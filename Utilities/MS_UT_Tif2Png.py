#! /usr/local/python-2.7.6/bin/python

import re, sys, os
import tifffile as tiff
from PIL import Image

if __name__ == "__main__":
    if len(sys.argv) == 3:
        tif_image_dir  = sys.argv[1]
        png_image_dir  = sys.argv[2]
        if not os.path.isdir(png_image_dir):
            os.mkdir(png_image_dir)
    else:
        sys.exit("\nUsage: MS_UT_Tif2Png.py <tif_image_dir> <png_image_dir>\n")

for file in os.listdir(tif_image_dir):
    old_path = os.path.join(tif_image_dir, file)
    new_path = os.path.join(png_image_dir, file.split(".")[0] + ".png")
    im = Image.open(old_path)
    im.save(new_path)

