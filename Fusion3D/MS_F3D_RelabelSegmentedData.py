#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#
import os, sys, re, optparse
import numpy
import h5py
from PIL import Image 

import MS_LIB_Options

ms_home = os.environ['MS_HOME']
ms_data = os.environ['MS_DATA']
ms_temp = os.environ['MS_TEMP']

# -----------------------------------------------------------------------

def get_image_data(input_data, input_type, options):
    image_data = []
    try:
        z = int(options.zmin) + int(options.node) - 1
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
    labels_file = os.path.join(ms_temp, input_label + "_labels.txt")
    labels = numpy.array(numpy.loadtxt(labels_file, delimiter='\t'))
    labels1 = []
    for i in range(0, labels.shape[0]):
        if labels[i][0] == z:
            labels1.append([labels[i, 1], labels[i, 2]])
    relabeled_image_data = numpy.zeros(image_data.shape, dtype=image_data.dtype)
    for i in range(0, len(labels1)):
        print "Changing label ", int(labels1[i][0]), " for label ", int(labels1[i][1]), \
              " on area ", (image_data == labels1[i][0]).sum()
        relabeled_image_data[image_data == labels1[i][0]] = labels1[i][1]
#       sys.stdout.flush()
    return relabeled_image_data

# -----------------------------------------------------------------------

def relabel_one_layer(input_data, input_type, options):  

    input_label    = "ms3_DVID"
    if input_type == "directory":
        input_label = "ms3_" + input_data[4:]
    else:
        input_label = "ms3_" + input_data.split('.')[0][4:]

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
#   img = Image.fromarray(relabeled_image_data).convert('RGB')
    img = Image.fromarray(relabeled_image_data)
    if options.verbose:
        print "Saving the relabeled image in file: ", output_path
    img.save(output_path)

# -----------------------------------------------------------------------

if __name__ == "__main__":

    usage = "\nUsage: %prog input_data input_type [options (-h to list)]\n"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = MS_LIB_Options.RelabelSegmentedData_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if len(args) == 2:
        input_data, input_type = args
        relabel_one_layer(input_data, input_type, options)
    else:
        print usage
        sys.exit(2)


