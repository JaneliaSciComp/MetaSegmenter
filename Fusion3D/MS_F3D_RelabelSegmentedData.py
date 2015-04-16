#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#
import os, sys, re, optparse
import numpy
import h5py
from PIL import Image 

ms_home = os.environ['MS_HOME']
ms_data = os.environ['MS_DATA']
ms_temp = os.environ['MS_TEMP']

# -----------------------------------------------------------------------

def relabel_command_line_parser(parser):
    parser.add_option("-D", "--debug",dest="debug", help="debugging; don't delete shell scripts", action="store_true", default=False)
    parser.add_option("-n", "--node",   dest="node", help="id of the cluster node to be used", metavar="node",  default=0)
    parser.add_option("-o", "--output_path",dest="output_path",help="output path",metavar="output_path",default="")
    parser.add_option("-v", "--verbose",action="store_true",dest="verbose",help="increase the verbosity level of output",default=False)
    parser.add_option("-z", "--zmin",dest="zmin",help="# of subsections in z direction", metavar="zmin", default=0)
    parser.add_option("-Z", "--zmax",dest="zmax",help="# of subsections in z direction", metavar="zmax", default=sys.maxint)
    return parser

# -----------------------------------------------------------------------

def get_image_data(input_data, input_type, options):
    image_data = []
    try:
        z = int(options.node) - 1
    except:
        z = options.node - 1
    if input_type == "file":
        input_path = os.path.join(ms_data, input_data)
        f  = h5py.File(input_path, 'r')
        key = f.keys()[0]
        image_data = numpy.transpose(f[key][z])
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
        image_data = numpy.zeros(data_shape, dtype="float")
        i = 0
        for file in files:
            im = misc.imread(file)
            image_data[:,:,i] = im[ymin:ymax, xmin:xmax,1]
            i = i + 1
    elif input_type == "DVID":
        connection = httplib.HTTPConnection("emdata2.int.janelia.org:80",\
                     timeout=5.0 )
        dvid_volume = voxels.VoxelsAccessor(connection, options.uuid, \
                                            options.dataset)
        image_data = numpy.asarray(dvid_volume[:, ymin:ymax, xmin:xmax, zmin])
        # Put data in C-contiguous order (last-index varies the fastest).
        image_data = image_data.copy(order='C')
    return image_data

# -----------------------------------------------------------------------

def relabel_image_data(image_data, input_label, options):
    print "image_data.dtype=", image_data.dtype
    z = int(options.zmin) + int(options.node) -1
    labels_file = os.path.join(ms_temp, input_label + "_labels" +\
                               "_z" + str(z+1) + ".txt")
    labels = numpy.loadtxt(labels_file, delimiter='\t')
    relabeled_image_data = numpy.zeros(image_data.shape, dtype=image_data.dtype)
    for i in range(1, len(labels)):
        relabeled_image_data[image_data == i] = labels[i]
    return relabeled_image_data

# -----------------------------------------------------------------------

def relabel_one_layer(input_data, input_type, options):  

    input_label    = "ms_DVID"
    if input_type == "directory":
        input_label = input_data
    else:
        input_label = input_data.split('.')[0]

    # Extract image data
    image_data = numpy.matrix(get_image_data(input_data, input_type, options))
    relabeled_image_data = relabel_image_data(image_data, input_label, options)
    
    # Store the relabeled data 
    z = int(options.zmin) + int(options.node) -1
    if len(options.output_path) == 0:
        output_path = os.path.join(ms_temp, input_label + "_relabeled" +\
                                   "_z" + str(z+1) + ".png")
    else:
        output_path = options.output_path
    img = Image.fromarray(relabeled_image_data).convert('RGB')
    if options.verbose:
        print "Saving the relabeled image in file: ", output_path
    img.save(output_path)

# -----------------------------------------------------------------------

if __name__ == "__main__":

    usage = "\nUsage: %prog input_data input_type [options (-h to list)]\n"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = relabel_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if len(args) == 2:
        input_data, input_type = args
        relabel_one_layer(input_data, input_type, options)
    else:
        print usage
        sys.exit(2)


