#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

from PIL import Image
import h5py
import numpy
import sys, os

if __name__ == "__main__":
    print "len(sys.argv) = ", len(sys.argv)
    if len(sys.argv) in [3, 4]:
        h5_name =     sys.argv[1]
        if os.path.isfile(h5_name):
            infiles  = [ h5_name ]
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
            for infile in os.listdir(h5_name):
                infiles.append(os.path.join(h5_name, infile))
                suffix = infile.split(".")[len(infile.split("."))-1]
                outfile = infile[:(len(infile) - len(suffix) -1)] + ".png"
                outfiles.append(os.path.join(outdir, outfile))
        else:
            sys.exit("Input type is not recognized")
        prob_id  = 0
        if len(sys.argv) == 4:
            prob_id = int(sys.argv[3])
    else:
        sys.exit("usage: MS_UT_H5ProbsToPng.py <h5file or dir_of_h5files> [ output_png_file or dir ] \n")

    for i in range(0, len(infiles)):
        print "\nConverting ", infiles[i], " to ", outfiles[i]
        f  = h5py.File(infiles[i], 'r')   
        try:
            data0 = numpy.array(f["volume/predictions"])
        except:
            try:
                data0 = f["probabilities"]
            except:
                data0 = f["stack"]
        print "data0.shape=", data0.shape

        if len(data0.shape) == 2:
            data = data0
        elif len(data0.shape) == 3:
            data = numpy.squeeze(data0[:,:,  prob_id])
        elif len(data0.shape) == 4:
            data = numpy.squeeze(data0[:,:,:,prob_id])
        else:
            sys.exit("Cannot process this data type")
        print "data.shape=", data.shape, " max(data)=", numpy.max(data), " data.dtype=", data.dtype
        max_data = numpy.max(data)
        if max_data > 0:
            img  = Image.fromarray(numpy.uint8(numpy.round(data/max_data * 255.)))
            img.save(outfiles[i])                 
        else:
            print "Skipping file ", infiles[i], " since max_data=", max_data
