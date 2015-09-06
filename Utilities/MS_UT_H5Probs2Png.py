#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

from PIL import Image
import h5py
import numpy
import sys, os

if __name__ == "__main__":

    if len(sys.argv) in [2, 3]:
        h5_name =     sys.argv[1]
        if os.path.isfile(h5_name):
            infiles  = [ h5_name ]
            outfiles = [ h5_name.split(".")[0] + ".png"]
            if len(sys.argv) == 3:
                outfiles = [ sys.argv[2] ]
        elif os.path.isdir(h5_name):
            outdir   = h5_name
            if len(sys.argv) == 3:
                outdir = sys.argv[2]
                if os.path.isfile(outdir):
                    sys.exit("The output directory cannot be existing file")
                if not os.path.isdir(outdir):
                    os.mkdir(outdir)
            infiles  = []
            outfiles = []
            for file in os.listdir(h5_name):
                infiles.append(os.path.join(h5_name, file))
                outfiles.append(os.path.join(outdir, file.split(".")[0] + ".png"))
        else:
            sys.exit("Input type is not recognized")
    else:
        sys.exit("usage: MS_UT_H5ToPng.py <h5file or dir_of_h5files> [ output_png_file or dir ] \n")

    for i in range(0, len(infiles)):
        f  = h5py.File(infiles[i], 'r')   
        data = f["probabilities"]
        img  = Image.fromarray(numpy.uint8(data/numpy.max(data)*255.))
        img.save(outfiles[i])                 
