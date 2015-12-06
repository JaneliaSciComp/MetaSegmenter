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
from scipy.ndimage.measurements import center_of_mass, find_objects

import MS_LIB_Dict
import MS_LIB_Options

ms_home = os.environ['MS_HOME']
ms_data = os.environ['MS_DATA']
ms_temp = os.environ['MS_TEMP']

# ----------------------------------------------------------------------

def compute_overlap(array_1, array_2, label1, label2, options):
    DL = float(options.dist_layers)  # by default distance between layers is 10 fold greater than a pixel size
    CS = float(options.critical_slope) # tan of the minimal angle between neural pili and z-layer

    # Determine the centers of mass of the two objects
    if options.debug:
        print "...computing the overlap for label1=", label1, " label2=", label2
        print "   shape1=", array_1.shape, " shape2=", array_2.shape
    center1 = numpy.int32(numpy.round(center_of_mass(array_1, labels=array_1, index=label1)))
    center2 = numpy.int32(numpy.round(center_of_mass(array_2, labels=array_2, index=label2)))
    y1_center = center1[0]
    x1_center = center1[1]
    y2_center = center2[0]
    x2_center = center2[1]

#   print "\ny1_center, x1_center, y2_center, x2_center=", y1_center, x1_center, y2_center, x2_center
    DC = numpy.sqrt(float(numpy.power(center1[0]-center2[0], 2) + \
                          numpy.power(center1[1]-center2[1], 2)))
    if options.debug:
        print "    DC=", DC
    
    if DL < DC * CS:
        return (0, 0, 0)

    # Determine the sizes of objects
    objects1 = find_objects(array_1, max_label=label1)
#   print "\nlen(objects1)=", len(objects1)
    loc1 = objects1[label1-1]
#   print "loc1=", loc1
    object1_shape = array_1[loc1].shape
    y1_min = loc1[0].start
    y1_max = loc1[0].stop
    x1_min = loc1[1].start
    x1_max = loc1[1].stop
#   print "y1_min,  y1_max, x1_min, x1_max=", y1_min,  y1_max, x1_min, x1_max

    objects2 = find_objects(array_2, max_label=label2)
#   print "len(objects2)=", len(objects2)
    loc2 = objects2[label2-1]
#   print "loc2=", loc2                          
    object2_shape = array_2[loc2].shape
#   print "\nobject1_shape=", object1_shape, " objects2_shape=", object2_shape
    y2_min = loc2[0].start
    y2_max = loc2[0].stop
    x2_min = loc2[1].start
    x2_max = loc2[1].stop
#   print "y2_min,  y2_max, x2_min, x2_max=", y2_min,  y2_max, x2_min, x2_max

    # Determine a min size of a parallelepoped that would
    # contain both the objects if their centers are placed at the center of the parallelepiped 
    ysize = max(y1_center - y1_min, y2_center -y2_min) + max(y1_max - y1_center, y2_max - y2_center)
    xsize = max(x1_center - x1_min, x2_center -x2_min) + max(x1_max - x1_center, x2_max - x2_center)
#   print "ysize=", ysize, " xsize=", xsize
    P1 = numpy.zeros([ysize, xsize])
    P2 = numpy.zeros([ysize, xsize])

    # Copy data from the original matrices to the parallelepoped1
#   print "P1.shape=", P1.shape, " P2.shape=", P2.shape
 
    P1_offset_y = max(y1_center - y1_min, y2_center -y2_min) - (y1_center - y1_min)  
    P2_offset_y = max(y1_center - y1_min, y2_center -y2_min) - (y2_center - y2_min)                    
    P1_offset_x = max(x1_center - x1_min, x2_center -x2_min) - (x1_center - x1_min)
    P2_offset_x = max(x1_center - x1_min, x2_center -x2_min) - (x2_center - x2_min)                    

#   print "P1_offset_y, P1_offset_x, P2_offset_y, P2_offset_x=", P1_offset_y, P1_offset_x, P2_offset_y, P2_offset_x

    P1[P1_offset_y:(P1_offset_y + y1_max - y1_min),
       P1_offset_x:(P1_offset_x + x1_max - x1_min)] = \
           array_1[y1_min : y1_max, x1_min : x1_max]
    P2[P2_offset_y:(P2_offset_y+object2_shape[0]),
       P2_offset_x:(P2_offset_x+object2_shape[1])] = \
           array_2[y2_min : y2_max, x2_min : x2_max]

    # Determin the overlap area
    overlap = ((P1 == label1) & (P2 == label2)).sum()
    area1   = (P1 == label1).sum()
    area2   = (P2 == label2).sum()
    projection_factor = DL / numpy.sqrt(numpy.power(DL,2) + numpy.power(DC,2))
    overlap = int(round(float(overlap) * projection_factor)) 
    print "final overlap=", overlap, "\n\n"
    return (overlap, area1, area2)

# ----------------------------------------------------------------------

def generate_overlap_and_area_matrices(image_fragment1, image_fragment2, options):
    fragment_matrix1 = numpy.array(image_fragment1)      
    fragment_matrix2 = numpy.array(image_fragment2)
    labels1 = numpy.uint32(numpy.unique(list(fragment_matrix1)))
    labels2 = numpy.uint32(numpy.unique(list(fragment_matrix2)))
    if options.verbose:
        print "len(labels1)=", len(labels1)
        print "len(labels2)=", len(labels2)

    overlap_matrix = numpy.zeros((len(labels1)+1, len(labels2)+1), dtype=int)
    area1_matrix   = numpy.zeros((len(labels1)+1, len(labels2)+1), dtype=int)
    area2_matrix   = numpy.zeros((len(labels1)+1, len(labels2)+1), dtype=int)
    if options.verbose:
        print "len(labels2)=", len(labels2)
    overlap_matrix[0 , 0] = 0
    overlap_matrix[1:, 0] = labels1
    overlap_matrix[0 ,1:] = labels2 
    area1_matrix[0 , 0] = 0
    area1_matrix[1:, 0] = labels1
    area1_matrix[0 ,1:] = labels2
    area2_matrix[0 , 0] = 0
    area2_matrix[1:, 0] = labels1
    area2_matrix[0 ,1:] = labels2
    for i in range(0, len(labels1)):
        for j in range(0, len(labels2)):
            if labels1[i] == 0 or labels2[j] == 0:
                continue
#           overlap_matrix[i+1, j+1] = (binary_class1 & binary_class2).sum()
            overlap_matrix[i+1, j+1], area1_matrix[i+1, j+1], area2_matrix[i+1, j+1] = \
                compute_overlap(fragment_matrix1, fragment_matrix2, \
                                labels1[i], labels2[j], options)

            if options.verbose and options.debug and overlap_matrix[i+1, j+1] > 0:
                print "i=", i, "/", len(labels1), " j=", j, "/", len(labels2), \
                      " overlap=", overlap_matrix[i+1, j+1], \
                      " area1=",     area1_matrix[i+1, j+1], \
                      " area2=",     area2_matrix[i+1, j+1]
    return (overlap_matrix, area1_matrix, area2_matrix)

# -----------------------------------------------------------------------

def read_input_data(input_data, input_type, ymin, ymax, xmin, xmax, \
                    zmin, zmax, options):
    if input_type == "file":
        input_path = os.path.join(ms_data, input_data)
        f  = h5py.File(input_path, 'r')
#       key = f.keys()[0]
#       image_data1 = numpy.transpose(f[key][zmin]) # read one chunk
#       image_data2 = numpy.transpose(f[key][zmax-1])
        dataset = numpy.transpose(f['/main'])
        data_shape = dataset.shape
        num_layers = dataset.shape[2]
        image_data1= dataset[:,:,zmin]
        image_data2= dataset[:,:,zmin+1]
        if options.verbose:
            print "num_layers=", num_layers
            print "numpy.max(image_data1)=", numpy.max(image_data1)
            print "numpy.max(image_data2)=", numpy.max(image_data2)
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
        image_data2 = numpy.asarray(dvid_volume[:, ymin:ymax, xmin:xmax, zmax-1])
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
        print "Storing the overlap matrix in file ", output_path
    numpy.savetxt(output_path, overlap_matrix, fmt='%10u', delimiter='\t')

# -----------------------------------------------------------------------

def output_area_matrix(area_matrix, input_label, x, y, z, options):
    # Extract output path
    if len(options.output_path) == 0:
        output_path = os.path.join(ms_temp, input_label + "_area" +\
                                   "_y" + str(y+1) +\
                                   "_x" + str(x+1) + "_z" + str(z+1) + ".txt")
    else:
        output_path = options.output_path

    # Ourput data
    if options.verbose:
        print "Storing the maxarea matrix in file ", output_path
    numpy.savetxt(output_path, area_matrix, fmt='%10u', delimiter='\t')

# -----------------------------------------------------------------------

def process_inputs(input_data, input_type, dict_node_xyz, options):
    if input_type == 'DVID':
        input_label = "ms3_DVID" 
    elif input_type == "directory":
        input_label = "ms3_" + input_data[4:] 
    else:
        input_label = "ms3_" + input_data.split('.')[0][4:]

    node = int(options.node)
    if node > 0:
        y, ymin, ymax, x, xmin, xmax, z = dict_node_xyz[node]
        zmin = z
        zmax = z+1
    else:
        ymin = 0
        ymax = int(options.ydim)
        xmin = 0
        xmax = int(options.xdim)
        zmin = int(options.zmin)
        zmax = int(options.zmax)

    if options.verbose:
        print "    node=", node, "y, ymin, ymax, x, xmin, xmax, z=",\
              [y, ymin, ymax, x, xmin, xmax, z]

    # Extract image data
    image_fragment1, image_fragment2 = \
        read_input_data(input_data, input_type, \
        ymin, ymax, xmin, xmax, zmin, zmax, options)

    if options.verbose:
        print "image_fragment1.shape=", image_fragment1.shape, \
             " image_fragment2.shape=", image_fragment2.shape
        print "\n...Generating overlap matrix..."
    overlap_matrix, area1_matrix, area2_matrix = \
        generate_overlap_and_area_matrices(image_fragment1, image_fragment2, \
                                             options) 
#   print "overlap_matrix=\n", overlap_matrix
    sparse_overlap_matrix = convert_to_sparse_format(overlap_matrix)
    output_overlap_matrix(sparse_overlap_matrix, input_label, \
                          x, y, zmin, options)

    sparse_area_matrix = convert_to_sparse_format2(overlap_matrix, \
                                                   area1_matrix, area2_matrix,)
    output_area_matrix(sparse_area_matrix, input_label, \
                          x, y, zmin, options)

# -----------------------------------------------------------------------

def convert_to_sparse_format(omatrix):
    sparse_matrix = [[ omatrix.shape[0]-1, omatrix.shape[1]-1, 0]]
    for i in range(1,omatrix.shape[0]):
        for j in range(1, omatrix.shape[1]):
            if omatrix[i][j] > 0: 
                sparse_matrix.append([omatrix[i,0],omatrix[0,j],omatrix[i,j]])
    print "sparse_matrix=\n", sparse_matrix
    return numpy.array(sparse_matrix)

# -----------------------------------------------------------------------

def convert_to_sparse_format2(omatrix, a1_matrix, a2_matrix):
    sparse_matrix = [[ omatrix.shape[0]-1, omatrix.shape[1]-1, 0, 0]]
    for i in range(1,omatrix.shape[0]):
        for j in range(1, omatrix.shape[1]):
            if omatrix[i][j] > 0:
                sparse_matrix.append([omatrix[i,0],   omatrix[0,j],\
                                      a1_matrix[i,j], a2_matrix[i,j]])
    print "sparse_matrix=\n", sparse_matrix
    return numpy.array(sparse_matrix)

# -----------------------------------------------------------------------

if __name__ == "__main__":

    usage = "\nUsage: %prog input_data input_type [ options ]\n"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = MS_LIB_Options.GenerateMatrices_command_line_parser(parser)
    (options, args) = parser.parse_args()
    
    if len(args) == 2:
        input_data, input_type = args
        input_dim = [int(options.ydim), int(options.xdim), int(options.zdim)]

        if input_type == 'DVID':
            input_label = "ms3_DVID"
        elif input_type == "directory":
            input_label = "ms3_" + os.path.basename(input_data)[4:]
        else:
            input_label = "ms3_" + os.path.basename(input_data).split('.')[0][4:]

        input_dim1 = [input_dim[0], input_dim[1], input_dim[2]-1];
        num_nodes, dict_node_xyz = \
            MS_LIB_Dict.map_node_to_xyz(input_dim1, input_label+"_fusion", ".txt", options)
        if options.verbose:
#           print "dict_node_xyz=", dict_node_xyz
            print "In MS_F3D_GenerateOverlapMatrix.py: num_nodes=", num_nodes
        if len(options.output_path) == 0 and len(options.node) == 0:
            print "\nPlease, specify the output name (with option -o)"
            sys.exit(2) 
        process_inputs(input_data, input_type, dict_node_xyz, options)
    else:
        print usage
        sys.exit(2)

