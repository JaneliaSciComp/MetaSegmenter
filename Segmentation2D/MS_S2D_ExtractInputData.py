#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#
import os
import sys, re, optparse
import tifffile as tiff
import numpy
from scipy import misc
import httplib
from pydvid import voxels, general

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
    parser.add_option("-n", "--node",   dest="node", help="id of the cluster node to be used", metavar="node",  default="")
    parser.add_option("-o", "--output_path",dest="output_path",help="output path",metavar="output_path",default="")
    parser.add_option("-r", "--uri",dest="uri",help="uri for DVID store",default='/api/repo/fe7/instance')
    parser.add_option("-U", "--unprocessed",action="store_true",dest="unprocessed", help="reprocess only the data for which an output file does not exist",default=False)
    parser.add_option("-v", "--verbose",action="store_true",dest="verbose",help="increase the verbosity level of output",default=False)
    parser.add_option("-W", "--width", dest="xdim", help="x-size (width) of an image", metavar="xdim", default="")
    parser.add_option("-X", "--nx",  dest="nx",  help="# of subsections in x direction", metavar="nx", default=1)
    parser.add_option("-Y", "--ny",  dest="ny",  help="# of subsections in y direction", metavar="ny", default=1)
    parser.add_option("-x", "--dx",  dest="dx",  help="# of scans for image overlap in x direction", metavar="dx", default=50)
    parser.add_option("-y", "--dy",  dest="dy",  help="# of scans for image overlap in y direction", metavar="dy", default=50)
    parser.add_option("-z", "--zmin",dest="zmin",help="# of subsections in z direction", metavar="zmin", default=0)
    parser.add_option("-Z", "--zmax",dest="zmax",help="# of subsections in z direction", metavar="zmax", default=sys.maxint)
    return parser

# ----------------------------------------------------------------------

# Given node id, determine z-layer id, ymin, ymax, xmin and xmax
def map_node_to_xyz(input_dim, input_label, options):
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
            ymax =       y_size + int(options.dy)
        elif y > 0 and y < max(y_range):
            ymin =  y   *y_size - int(options.dy)
            ymax = (y+1)*y_size + int(options.dy)
        else:
            ymin =  y   *y_size - int(options.dy)
            ymax =  input_dim[0]
        for x in x_range:
            if x == 0:
                xmin =  0
                xmax =       x_size + int(options.dx)
            elif x > 0 and x < max(x_range):
                xmin =  x   *x_size - int(options.dx)
                xmax = (x+1)*x_size + int(options.dx)
            else:
                xmin =  x   *x_size - int(options.dx)
                xmax =  input_dim[1]
            for z in z_range:
                if not options.unprocessed:
                    node += 1
                    dict_node_xyz[node] = \
                        [y, ymin, ymax, x, xmin, xmax, z]
                else:
                    output_file = \
                        os.path.join(ms_temp, input_label + "_y" + str(y+1) +\
                                     "_x" + str(x+1) + "_z" + str(z+1) + \
                                     ".png")
                    if not os.path.exists(output_file):
                        node += 1
                        dict_node_xyz[node] = \
                            [y, ymin, ymax, x, xmin, xmax, z]
                        if options.verbose:
                            print "node=",node," missing output file=",output_file
    return (node, dict_node_xyz)

# -----------------------------------------------------------------------

def extract_image_data(input_data, input_type, ymin, ymax,\
                       xmin, xmax, zmin, zmax, options):
    # Extract image data
    if input_type == "file":
        if re.search(".h5", input_data): 
            f  = h5py.File(input_data, 'r')
            key = f.keys()[0]
            data = numpy.transpose(f[key])
            data_shape = data.shape
            if len(data.shape) == 3:
                image_data = data[ymin:ymax, xmin:xmax, zmin:zmax]
            else:
                image_data = data[ymin:ymax, xmin:xmax, zmin:zmax, 1]
        elif re.search(".tif", input_data):
            if options.verbose:
                print "tiff input_data=", input_data
            data = numpy.transpose(tiff.imread(input_data))
            data_shape = data.shape
            if len(data.shape) == 3:
                image_data = data[ymin:ymax, xmin:xmax, zmin:zmax]
            else:
                image_data = data[ymin:ymax, xmin:xmax, zmin:zmax, 1]
            print "Input data shape=", image_data.shape
        else:
            sys.exit("\nUnsupported type of input file. Exit.")
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
        image_data = numpy.zeros((ymax-ymin,xmax-xmin), dtype="float")
        i = 0
        for file in files:
            im = misc.imread(file)
            if options.verbose:
                print "file=", file, " data_shape=", data_shape, " im.shape=", im.shape
            if len(im.shape) == 2:
                image_data = im[ymin:ymax, xmin:xmax]
            else:
                image_data = im[ymin:ymax, xmin:xmax,1]
            i = i + 1
    elif input_type == "DVID":
        connection = httplib.HTTPConnection("emdata2.int.janelia.org:80",\
                     timeout=5.0 )
        dvid_volume = voxels.VoxelsAccessor(connection, options.uuid, \
                                            options.dataset)
        image_data = numpy.asarray(dvid_volume[:, ymin:ymax, xmin:xmax, zmin:zmax])
        image_data = image_data.copy(order='C')
    return image_data

# -----------------------------------------------------------------------

def process_inputs(input_data, input_type, input_label, dict_node_xyz, options):

    node = int(options.node)
    if node > 0:
        y, ymin, ymax, x, xmin, xmax, z = dict_node_xyz[node]
        zmin = z
        zmax = z + 1
    else:
        ymin = 0
        ymax = int(options.ydim)
        xmin = 0
        xmax = int(options.xdim)
        zmin = max(0,               int(options.zmin))
        zmax = min(int(options.zdim,int(options.zmax)))
    if options.verbose:
        print "    node=", node, "y, ymin, ymax, x, xmin, xmax, z=",\
              [y, ymin, ymax, x, xmin, xmax, z]
        print "    ysize=", ymax - ymin

    # Extract image data
    image_data = extract_image_data(input_data, input_type, ymin, ymax,\
                                    xmin, xmax, zmin, zmax, options)
    # Extract output path
    if len(options.output_path) == 0:
        output_path = os.path.join(ms_temp, input_label + "_y" + str(y+1) +\
                                   "_x" + str(x+1) + "_z" + str(z+1) + ".png")
    else:
        output_path = options.output_path

    # Ourput data
    if options.verbose:
        print "Image_data shape=", image_data.shape
        print "Writing image data to file", output_path
    misc.imsave(output_path, image_data)   

# -----------------------------------------------------------------------

if __name__ == "__main__":

    usage = "\nUsage: %prog input_data input_type [ options ]\n"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = get_data_command_line_parser(parser)
    (options, args) = parser.parse_args()
    
    if len(args) == 2:
        input_data, input_type = args
        input_dim = [int(options.ydim), int(options.xdim), int(options.zdim)]
        input_label = "ms_DVID"
        if input_type in ["file", "directory"]:
            input_label = "ms_" + input_data.split('.')[0]
        num_nodes, dict_node_xyz = \
            map_node_to_xyz(input_dim, input_label, options)
        if options.verbose:
            print "num_nodes=", num_nodes
            print "input_dim=", input_dim
        if len(options.output_path) == 0 and len(options.node) == 0:
            print "\nPlease, specify the output name (with option -o)"
            sys.exit(2) 
        if num_nodes > 0:
            process_inputs(input_data, input_type, input_label, dict_node_xyz, options)
    else:
        print usage
        sys.exit(2)

