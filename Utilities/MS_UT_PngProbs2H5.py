#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

from PIL import Image
import h5py
import numpy
import sys, os

if __name__ == "__main__":

    transpose = 1
    
    inverse_probs = 1

    if len(sys.argv) in [2, 3]:
        png_name =     sys.argv[1]
        if os.path.isfile(png_name):
            infiles  = [ png_name ]
            outfiles = [ png_name.split(".")[0] + ".h5"]
            if len(sys.argv) == 3:
                outfiles = [ sys.argv[2] ]
        elif os.path.isdir(png_name):
            outdir   = png_name
            if len(sys.argv) == 3:
                outdir = sys.argv[2]
                if os.path.isfile(outdir):
                    sys.exit("The output directory cannot be existing file")
                if not os.path.isdir(outdir):
                    os.mkdir(outdir)
            infiles  = []
            outfiles = []
            for infile in os.listdir(png_name):
                infiles.append(os.path.join(png_name, infile))
                suffix = infile.split(".")[len(infile.split("."))-1]
                outfile = outfile = infile[:(len(infile) - len(suffix) -1)] + ".h5"
                outfiles.append(os.path.join(outdir, outfile))
        else:
            sys.exit("Input type is not recognized")
    else:
        sys.exit("usage: MS_UT_PngProbs2H5.py <png_file or dir_of_png_files> [ output_h5_file or dir ] \n")

    for i in range(0, len(infiles)):
        data = numpy.double(numpy.array(Image.open(infiles[i])))
        if transpose:
            data = numpy.transpose(data)
        if inverse_probs:
            data = 1. - data
        data = data/numpy.max(data)
        f = h5py.File(outfiles[i], 'w')
        f.create_dataset('probabilities', data.shape, data = data, dtype = numpy.double)
