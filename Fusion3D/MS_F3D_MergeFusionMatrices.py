#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#
import os, sys, re, optparse
import h5py
import numpy
from   zibopt import scip

ms_home = os.environ['MS_HOME']
ms_data = os.environ['MS_DATA']
ms_temp = os.environ['MS_TEMP']

# -----------------------------------------------------------------------

def merge_command_line_parser(parser):
    parser.add_option("-D", "--debug",dest="debug", help="debugging; don't delete shell scripts", action="store_true", default=False)
    parser.add_option("-n", "--node",   dest="node", help="id of the cluster node to be used", metavar="node",  default=0)
    parser.add_option("-o", "--output_path",dest="output_path",help="output path",metavar="output_path",default="")
    parser.add_option("-v", "--verbose",action="store_true",dest="verbose",help="increase the verbosity level of output",default=False)
    parser.add_option("-X", "--nx",  dest="nx",  help="# of subsections in x direction", metavar="nx", default=1)
    parser.add_option("-Y", "--ny",  dest="ny",  help="# of subsections in y direction", metavar="ny", default=1)
    parser.add_option("-z", "--zmin",dest="zmin",help="# of subsections in z direction", metavar="zmin", default=0)
    parser.add_option("-Z", "--zmax",dest="zmax",help="# of subsections in z direction", metavar="zmax", default=sys.maxint)
    return parser

# -----------------------------------------------------------------------

def read_input_data(input_data, input_type, zmin, zmax, options):
    if input_type == "file":
        input_path = os.path.join(ms_data, input_data)
        f  = h5py.File(input_path, 'r')
        key = f.keys()[0]
        image_data1 = numpy.transpose(f[key][zmin]) # read one chunk
        image_data2 = numpy.transpose(f[key][zmax])
        data_shape = image_data1.shape
    elif input_type == "directory":
        files = []
        i = 0
        dir_path = os.path.join(ms_data,input_data)
        files_list = os.listdir(dir_path)
        files_list.sort()
        files = []
        for file in files_list:
            file_path = os.path.join(dir_path, file)
            if i in range(int(zmin), int(zmax)):
                files.append(os.path.join(input_data, file_path))
            i = i+1
        im1 =  misc.imread(files[0])
        data_shape = misc.imread(files[0]).shape
        image_data1 = numpy.zeros(data_shape, dtype="float")
        image_data2 = numpy.zeros(data_shape, dtype="float")
        i = 0
        for file in files:
            im = misc.imread(file)
            image_data1[:,:,i] = im[:,:,1]
            i = i + 1
    elif input_type == "DVID":
        connection = httplib.HTTPConnection("emdata2.int.janelia.org:80",\
                     timeout=5.0 )
        dvid_volume = voxels.VoxelsAccessor(connection, options.uuid, \
                                            options.dataset)
        image_data1 = numpy.asarray(dvid_volume[:,:,:, zmin])
        image_data2 = numpy.asarray(dvid_volume[:,:,:, zmax])
        # Put data in C-contiguous order (last-index varies the fastest).
        image_data1 = image_data1.copy(order='C')
        image_data2 = image_data2.copy(order='C')
    return (image_data1, image_data2)              

# -----------------------------------------------------------------------

def produce_unified_matrices(input_data, input_type, options):  

    input_label    = "ms_DVID"
    if input_type == "directory":
        input_label = input_data
    else:
        input_label = input_data.split('.')[0]

    node = int(options.node)
    if options.verbose:
        print "node=", node, " type=", type(node)

    # Extract image data
    image_data1, image_data2 = \
        read_input_data(input_data, input_type, \
        int(options.zmin), int(options.zmax), options)

    z        = int(options.node)
    nx       = int(options.nx)
    ny       = int(options.ny)
    matrix1  = numpy.matrix(image_data1)
    matrix2  = numpy.matrix(image_data2)
    nlabels1 = int(matrix1.max())
    nlabels2 = int(matrix2.max())
    print "nlabels1=", nlabels1, " nlabels2=", nlabels2
 
    # Overkap matrix
    unified_omatrix = numpy.zeros((nlabels1+1,nlabels2+1), dtype="int")
    print "unified_omatrix.shape=", unified_omatrix.shape
    unified_omatrix[0,1:] = range(1, nlabels2+1)
    unified_omatrix[1:,0] = range(1, nlabels1+1)

    # Fusion matrix
    unified_fmatrix = numpy.zeros((nlabels1+1,nlabels2+1), dtype="int")
    unified_fmatrix[0,1:] = range(1, nlabels2+1)
    unified_fmatrix[1:,0] = range(1, nlabels1+1)

    output_fpath = os.path.join(ms_temp, input_label + "_fusion" +\
                                "_z" + str(z) + ".txt")
    output_opath = os.path.join(ms_temp, input_label + "_overlap" +\
                                "_z" + str(z) + ".txt")
    for y in range(1, ny):
        for x in range(1, nx):
            input_opath = os.path.join(ms_temp, input_label + "_overlap" +\
                                       "_y" + str(y) +\
                                       "_x" + str(x) + "_z" + str(z) + ".txt")
            input_fpath = os.path.join(ms_temp, input_label + "_fusion" +\
                                       "_y" + str(y) +\
                                       "_x" + str(x) + "_z" + str(z) + ".txt")    
            omatrix = numpy.loadtxt(input_opath, delimiter='\t')
            fmatrix = numpy.loadtxt(input_fpath, delimiter='\t') 
            print "fmatrix.shape=", fmatrix.shape
            print "omatrix.shape=", omatrix.shape
            my_shape  = fmatrix.shape
            for k in range(1, my_shape[0]):
                for l in range(1, my_shape[1]):
                    # Coordinates in the unified matrix:
                    i = int(fmatrix[k][0])
                    j = int(fmatrix[0][l])
#                   print "x=", x, " y=", y, " i=", i, " j=", j, " k=", k, " l=", l
                    # Update the unified matrix only if the weight
                    # (=overlap) of element in the input matrix 
                    # is higher than in current unified matrix
                    if  unified_omatrix[i][j] < omatrix[k][l]:
                        unified_fmatrix[i][j] = fmatrix[k][l]
                        unified_omatrix[i][j] = omatrix[k][l]
    if options.verbose:
        print "Saving the fusion matrix in file: ", output_fpath
    numpy.savetxt(output_fpath, unified_fmatrix, fmt='%10u', delimiter='\t')
    if options.debug:
        if options.verbose:
            print "Saving the fusion matrix in file: ", output_fpath
        numpy.savetxt(output_fpath, unified_fmatrix, fmt='%10u', delimiter='\t')

# -----------------------------------------------------------------------

if __name__ == "__main__":

    usage = "\nUsage: %prog input_data input_type [options (-h to list)]\n"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = merge_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if len(args) == 2:
        input_data, input_type = args
        produce_unified_matrices(input_data, input_type, options)
    else:
        print usage
        sys.exit(2)


