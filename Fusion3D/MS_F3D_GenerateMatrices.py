#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#
import os
import sys, re, optparse
import numpy
from scipy import misc
import httplib
import h5py
from pydvid import voxels, general
from   zibopt import scip

ms_home = os.environ['MS_HOME']
ms_data = os.environ['MS_DATA']
ms_temp = os.environ['MS_TEMP']

# -----------------------------------------------------------------------

def get_data_command_line_parser(parser):
    parser.add_option("-d", "--dataset", dest="dataset", help="dataset in DVID store",default='test1')
    parser.add_option("-D", "--debug",dest="debug", help="debugging; don't delete shell scripts", action="store_true", default=False)
    parser.add_option("-i", "--uuid",dest="uuid",help="uuid for DVID store",default='fe7')
    parser.add_option("-L", "--length", dest="ydim", help="y-size (length) of an image", metavar="ydim", default="")
    parser.add_option("-H", "--height", dest="zdim", help="z-size (height, or # layers) of an image stack", metavar="zdim", default="")
    parser.add_option("-n", "--node",   dest="node", help="id of the cluster node to be used", metavar="node",  default=0)
    parser.add_option("-o", "--output_path",dest="output_path",help="output path",metavar="output_path",default="")
    parser.add_option("-r", "--uri",dest="uri",help="uri for DVID store",default='/api/repo/fe7/instance')
    parser.add_option("-v", "--verbose",action="store_true",dest="verbose",help="increase the verbosity level of output",default=False)
    parser.add_option("-W", "--width", dest="xdim", help="x-size (width) of an image", metavar="xdim", default="")
    parser.add_option("-X", "--nx",  dest="nx",  help="# of subsections in x direction (default=1)", metavar="nx", default=1)
    parser.add_option("-Y", "--ny",  dest="ny",  help="# of subsections in y direction (default=1)", metavar="ny", default=1)
    parser.add_option("-x", "--dx",  dest="dx",  help="# of scans for image overlap in x direction", metavar="dx", default=50)
    parser.add_option("-y", "--dy",  dest="dy",  help="# of scans for image overlap in y direction", metavar="dy", default=50)
    parser.add_option("-z", "--zmin",dest="zmin",help="# of subsections in z direction", metavar="zmin", default=0)
    parser.add_option("-Z", "--zmax",dest="zmax",help="# of subsections in z direction", metavar="zmax", default=sys.maxint)
    return parser

# ----------------------------------------------------------------------

# Given node id, determine z-layer id, ymin, ymax, xmin and xmax
def map_node_to_xyz(input_dim, options):
    z_range = range(max(0,           int(options.zmin)),\
                    min(input_dim[2],int(options.zmax)))
    y_range = range(0, int(options.ny))
    x_range = range(0, int(options.nx))

    y_size = int(round(float(input_dim[0])/float(options.ny)))
    x_size = int(round(float(input_dim[1])/float(options.nx)))
    if options.verbose:
        print "y_size=", y_size, " x_size=", x_size
    if y_size <= 2*int(options.dy):
        print "y overlap ", options.dy, " exceeds half of y_size ", y_size
        sys.exit()
    if x_size <= 2*int(options.dx):
        print "x overlap ", options.dx, " exceeds half of x_size ", x_size
        sys.exit()
    dict_node_xyz = {}
    node = 0
    for y in y_range:
        if y == 0:
            ymin =  0
            ymax =       y_size 
            if options.ny > 1:
                ymax =   y_size + int(options.dy)
        elif y > 0 and y < max(y_range):
            ymin =  y   *y_size - int(options.dy)
            ymax = (y+1)*y_size + int(options.dy)
        else:
            ymin =  y   *y_size - int(options.dy)
            ymax =  input_dim[0]
        for x in x_range:
            if x == 0:
                xmin =  0
                xmax =       x_size
                if options.nx > 1:
                    xmax =   x_size + int(options.dx)
            elif x > 0 and x < max(x_range):
                xmin =  x   *x_size - int(options.dx)
                xmax = (x+1)*x_size + int(options.dx)
            else:
                xmin =  x   *x_size - int(options.dx)
                xmax =  input_dim[1]
            for z in z_range:
                node += 1
                dict_node_xyz[node] = \
                    [y, ymin, ymax, x, xmin, xmax, z, z+1]
    return (node, dict_node_xyz)

# -----------------------------------------------------------------------

def generate_overlap_matrix(image_fragment1, image_fragment2, options):
    fragment_matrix1 = numpy.matrix(image_fragment1)      
    fragment_matrix2 = numpy.matrix(image_fragment2)
    max_label1 = fragment_matrix1.max()
    max_label2 = fragment_matrix2.max()
    if options.verbose:
        print "max_label1=", max_label1, " max_label2=", max_label2    
    labels1 = []
    labels2 = []
    for i in range(1, int(max_label1)+1):
        count = (fragment_matrix1 == i).sum()
        if count > 0:
            labels1.append(i)
    if options.verbose:
        print "len(labels1)=", len(labels1)
    for i in range(1, int(max_label2)+1):
        count = (fragment_matrix2 == i).sum()
        if count > 0:
            labels2.append(i)
    overlap_matrix = numpy.zeros((len(labels1)+1, len(labels2)+1), dtype=int)
    if options.verbose:
        print "len(labels2)=", len(labels2)
    overlap_matrix[0 , 0] = 0
    overlap_matrix[1:, 0] = labels1
    overlap_matrix[0 ,1:] = labels2 
    for i in range(0, len(labels1)):
        for j in range(0, len(labels2)):
            overlap_matrix[i+1, j+1] = \
                ((fragment_matrix1 == int(labels1[i]))  & \
                 (fragment_matrix2 == int(labels2[j]))).sum()
            if options.verbose and options.debug and overlap_matrix[i+1, j+1] > 0:
                print "i=", i, "/", len(labels1), " j=", j, "/", len(labels2), \
                      " overlap=", overlap_matrix[i+1, j+1]
    return overlap_matrix

# -----------------------------------------------------------------------

def read_input_data(input_data, input_type, ymin, ymax, xmin, xmax, \
                    zmin, zmax, options):
    if input_type == "file":
        input_path = os.path.join(ms_data, input_data)
        f  = h5py.File(input_path, 'r')
        key = f.keys()[0]
        image_data1 = numpy.transpose(f[key][zmin]) # read one chunk
        image_data2 = numpy.transpose(f[key][zmax])
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
            image_data1[:,:,i] = im[ymin:ymax, xmin:xmax,1]
            i = i + 1
    elif input_type == "DVID":
        connection = httplib.HTTPConnection("emdata2.int.janelia.org:80",\
                     timeout=5.0 )
        dvid_volume = voxels.VoxelsAccessor(connection, options.uuid, \
                                            options.dataset)
        image_data1 = numpy.asarray(dvid_volume[:, ymin:ymax, xmin:xmax, zmin])
        image_data2 = numpy.asarray(dvid_volume[:, ymin:ymax, xmin:xmax, zmax])
        # Put data in C-contiguous order (last-index varies the fastest).
        image_data1 = image_data1.copy(order='C')
        image_data2 = image_data2.copy(order='C')
    image_fragment1 = image_data1[ymin:ymax, xmin:xmax]
    image_fragment2 = image_data2[ymin:ymax, xmin:xmax]
    return (image_fragment1, image_fragment2)

# -----------------------------------------------------------------------

def output_overlap_matrix(overlap_matrix, input_label, x, y, z, options):
    # Extract output path
    if len(options.output_path) == 0:
        output_path = os.path.join(ms_temp, input_label + "_overlap" +\
                                   "_y" + str(y+1) +\
                                   "_x" + str(x+1) + "_z" + str(z+1) + ".txt")
    else:
        output_path = options.output_path

    # Ourput data
    if options.verbose:
        print "Writing overlap_matrix to file", output_path

    if options.verbose:
        print "Storing the overlap matrix in file ", output_path
    numpy.savetxt(output_path, overlap_matrix, fmt='%10u', delimiter='\t')
    if options.debug:
        a = numpy.loadtxt(output_path, delimiter='\t')
        print "overlap_matrix[0:4,0:4]=", a[0:4,0:4]

# -----------------------------------------------------------------------

def generate_fusion_matrix(overlap_matrix, options):
    solver = scip.solver()# quiet=False)
    matrix_shape = overlap_matrix.shape
    fusion_matrix = numpy.zeros(matrix_shape)
    # Setting the boundary elements, which are label ids
    for i in range(0, matrix_shape[0]):
        fusion_matrix[i][0] = overlap_matrix[i][0]
    for j in range(0, matrix_shape[1]):
        fusion_matrix[0][j] = overlap_matrix[0][j]

    # Define SCIP variables
    dict_map = {}
    variables = []
    k = 0
    for i in range(1, matrix_shape[0]):
        for j in range(1, matrix_shape[1]):
            if overlap_matrix[i][j] > 0:
                variables.append(solver.variable(scip.BINARY))
                dict_map[k] = [i, j]
                k = k + 1

    # Define SCIP constraints
    # 1) sum of elements of fusion_matrix in each row is <= 2
    for i in range(1, matrix_shape[0]):
        nz_linear_inds = []
        for c in range(1, matrix_shape[1]):
            if overlap_matrix[i][c] == 0:
                continue
            for key, value in dict_map.iteritems():
                if [i, c] == value:
                    nz_linear_inds.append(key)
        if len(nz_linear_inds) > 0:
            if options.verbose and options.debug:
                print "i=", i, " nz_linear_inds=", nz_linear_inds
            solver.constraint( sum(variables[k] for k in nz_linear_inds) <= 2 )# allow binary branching

    # 2) sum of elements of fusion_matrix in each column is <= 2
    t_overlap_matrix = numpy.matrix.transpose(overlap_matrix)
    for i in range(1, matrix_shape[1]):
        nz_linear_inds = []
        for r in range(1, matrix_shape[0]):
            if t_overlap_matrix[i][r] == 0:
                continue
            for key, value in dict_map.iteritems():
                if [i, r] == value:
                    nz_linear_inds.append(key)
        if len(nz_linear_inds) > 0:
            if options.verbose and options.debug:
                print "j=", i, " nz_linear_inds=", nz_linear_inds
            solver.constraint( sum(variables[k] for k in nz_linear_inds) <= 2 ) # allow binary branching

    # Define the objective
#   print "dict_map.keys()=", dict_map.keys()
#   print "dict_map.values()=", dict_map.values()
    if options.verbose and options.debug:
        print "len(dict_map.keys())=", len(dict_map.keys())
    solution = solver.maximize(objective =  sum(overlap_matrix[dict_map[k][0]][dict_map[k][1]]*variables[k] \
                                            for k in dict_map.keys()))

    # Restore the fusion_matrix
    for i in range(1, matrix_shape[0]):
        for j in range(1, matrix_shape[1]):
            if overlap_matrix[i][j] > 0:
                for key, value in dict_map.iteritems():
                    if [i, j] == value:
                        fusion_matrix[i][j] = solution[variables[key]]
                if fusion_matrix[i][j] > 0:
                    if options.verbose and options.debug:
                        print "i=", i, " j=", j, " fusion_matrix=", fusion_matrix[i][j]
    return fusion_matrix

# -----------------------------------------------------------------------

def output_fusion_matrix(fusion_matrix, input_label, x, y, z, options):

    if len(options.output_path) == 0:
        output_path = os.path.join(ms_temp, input_label + "_fusion" +\
                                   "_y" + str(y+1) +\
                                   "_x" + str(x+1) + "_z" + str(z+1) + ".txt")
    else:
        output_path = options.output_path
    numpy.savetxt(output_path, fusion_matrix, fmt='%10u', delimiter='\t')
    if options.verbose:
        print "Saving the fusion matrix in file: ", output_path
    numpy.savetxt(output_path, fusion_matrix, fmt='%10u', delimiter='\t')

# -----------------------------------------------------------------------

def process_inputs(input_data, input_type, dict_node_xyz, options):
    if input_type == 'DVID':
        input_label = "ms_DVID" 
    elif input_type == "directory":
        input_label = input_data 
    else:
        input_label = input_data.split('.')[0]

    node = int(options.node)
    if options.verbose:
        print "node=", node, " type=", type(node)
    if node > 0:
        y, ymin, ymax, x, xmin, xmax, zmin, zmax = dict_node_xyz[node]
    else:
        ymin = 0
        ymax = int(options.ydim)
        xmin = 0
        xmax = int(options.xdim)
        zmin = int(options.zmin)
        zmax = int(options.zmax)

    if options.verbose:
        print "    node=", node, "y, ymin, ymax, x, xmin, xmax, zmin, zmax=",\
              [y, ymin, ymax, x, xmin, xmax, zmin, zmax]

    # Extract image data
    image_fragment1, image_fragment2 = \
        read_input_data(input_data, input_type, \
        ymin, ymax, xmin, xmax, zmin, zmax, options)

    if options.verbose:
        print "image_fragment1.shape=", image_fragment1.shape, \
             " image_fragment2.shape=", image_fragment2.shape
        print "\n...Generating overlap matrix..."
    overlap_matrix = generate_overlap_matrix(image_fragment1, image_fragment2, \
                                             options) 
    output_overlap_matrix(overlap_matrix, input_label, x, y, zmin, options)

    fusion_matrix = generate_fusion_matrix(overlap_matrix, options)

    output_fusion_matrix(fusion_matrix, input_label, x, y, zmin, options)

# -----------------------------------------------------------------------

if __name__ == "__main__":

    usage = "\nUsage: %prog input_data input_type [ options ]\n"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = get_data_command_line_parser(parser)
    (options, args) = parser.parse_args()
    
    if len(args) == 2:
        input_data, input_type = args
        input_dim = [int(options.ydim), int(options.xdim), int(options.zdim)]
        num_nodes, dict_node_xyz = \
            map_node_to_xyz(input_dim, options)
        if options.verbose:
            print "In MS_F3D_GenerateOverlapMatrix.py: num_nodes=", num_nodes
        if len(options.output_path) == 0 and len(options.node) == 0:
            print "\nPlease, specify the output name (with option -o)"
            sys.exit(2) 
        process_inputs(input_data, input_type, dict_node_xyz, options)
    else:
        print usage
        sys.exit(2)

