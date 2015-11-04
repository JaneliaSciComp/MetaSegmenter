#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#
import os, sys, re, optparse
import h5py
import numpy
from   zibopt import scip

import MS_LIB_Options

ms_home = os.environ['MS_HOME']
ms_data = os.environ['MS_DATA']
ms_temp = os.environ['MS_TEMP']

# -----------------------------------------------------------------------

def read_input_data(input_data, input_type, zmin, zmax, options):
    if input_type == "file":
        input_path = os.path.join(ms_data, input_data)
        if options.verbose:
            print "input_path=", input_path, "\nzmin=", zmin
        f  = h5py.File(input_path, 'r')
        key = f.keys()[0]
        image_data1 = numpy.transpose(f[key][zmin  ]) # read one chunk
        image_data2 = numpy.transpose(f[key][zmin+1])
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
        image_data2 = numpy.asarray(dvid_volume[:,:,:, zmax-1])
        # Put data in C-contiguous order (last-index varies the fastest).
        image_data1 = image_data1.copy(order='C')
        image_data2 = image_data2.copy(order='C')
    return (image_data1, image_data2)              

# -----------------------------------------------------------------------

def check_fusion_matrix(fmatrix):
    print "\nFusion matrix shape=", fmatrix.shape
    print "num_cols=", len(fmatrix[0])
    print "num_rows=", len(fmatrix[:,0])
    sum_col = [0,0,0]
    sum_row = [0,0,0]
    my_shape = fmatrix.shape
    # Explore rows    
    for i in range(1, my_shape[0]):
        row = fmatrix[i][1:]
        s = (row == 1).sum()
        if s <= 2:
            sum_col[s] +=1
        else:
            print "row #", i, " sum=", s
        if not len(row) == my_shape[1]-1:
            print "...row ", i, " contains ", len(row), " elements, not", my_shape[1]-1
    print "sum_col=", sum_col
    # Explore columns
    for j in range(1, my_shape[1]):
        col = fmatrix[1:,j]
        s = (col == 1).sum()
        if s <= 2:
            sum_row[s] +=1
        else:
            print "col #", j, " sum=", s
        if not len(col) == my_shape[0]-1:
            print "...col ", j, " contains ", len(col), " elements, not", my_shape[0]-1
    print "sum_row=", sum_row

# -----------------------------------------------------------------------

def check_sums(i, j, unified_fmatrix):
    if sum(unified_fmatrix[i,1:]) > 2:
        print "loss of constraint for row label", i
    if sum(unified_fmatrix[1:,j]) > 2:
        print "loss of constraint for col label", j

# -----------------------------------------------------------------------

def produce_unified_omatrix(input_data, input_type, options):  

    input_label    = "ms3_DVID"
    if input_type == "directory":
        input_label = "ms3_" + input_data[4:]
    else:
        input_label = "ms3_" + input_data.split('.')[0][4:]

    node = int(options.node)
    z    = node - 1
    if options.verbose:
        print "node=", node, " type=", type(node)

    # Extract image data
    image_data1, image_data2 = \
        read_input_data(input_data, input_type, \
        z, z+1, options)

    nx       = int(options.nx)
    ny       = int(options.ny)
    matrix1  = numpy.matrix(image_data1)
    matrix2  = numpy.matrix(image_data2)
    if options.verbose:
        print "matrix1.shape=", matrix1.shape, " matrix1.dtype=", matrix1.dtype
        print "matrix2.shape=", matrix2.shape, " matrix2.dtype=", matrix2.dtype
     
    labels1  = numpy.unique(list(matrix1))
    labels2  = numpy.unique(list(matrix2))
    print "nlabels1=", len(labels1), " nlabels2=", len(labels2)
    print "labels1=", labels1 
    print "labels2=", labels2

    # Overkap matrix
    unified_omatrix = numpy.zeros((len(labels1)+1,len(labels2)+1), dtype="int")
    print "unified_omatrix.shape=", unified_omatrix.shape
    unified_omatrix[0][1:] = labels2
    unified_omatrix[1:,0] = labels1

    output_opath = os.path.join(ms_temp, input_label + "_overlap" +\
                                "_z" + str(z+1) + ".txt")
    for y in range(0, ny):
        for x in range(0, nx):
            input_opath = os.path.join(ms_temp, input_label + "_overlap" +\
                                       "_y" + str(y+1) +\
                                       "_x" + str(x+1) + "_z" + str(z+1) + ".txt")
            input_apath = os.path.join(ms_temp, input_label + "_maxarea" +\
                                       "_y" + str(y+1) +\
                                       "_x" + str(x+1) + "_z" + str(z+1) + ".txt")
            try:    # inout matrix may be missing due to junk data
                if options.verbose:
                    print "input_opath=", input_opath
                omatrix = numpy.uint32(numpy.loadtxt(input_opath, delimiter='\t'))
                amatrix = numpy.uint32(numpy.loadtxt(input_apath, delimiter='\t'))
                if options.verbose:
                    print "\ninput_opath=", input_opath
                    print "    omatrix.shape=", omatrix.shape  
#                   print "    omatrix=", omatrix
                for k in range(1, omatrix.shape[0]):
                    li = int(omatrix[k][0]) # label of upper layer
                    lj = int(omatrix[k][1]) # label of lower layer
                    if options.debug:
                        print "   k=", k, " li=", li, " lj=", lj
                    overlap = omatrix[k][2]
                    maxarea = amatrix[k][2]
                    if options.debug:
                        print "   overlap=", overlap, " maxarea=", maxarea
                    i = list(labels1).index(int(li)) + 1 # row of a unified_omatrix which corresponds to li          
                    j = list(labels2).index(int(lj)) + 1 # row of a unified_omatrix which corresponds to lj          
                    if options.debug:
                        print "   i=", i, " j=", j
                    if overlap > max(unified_omatrix[i][j], int(options.overlap_area)) and \
                       overlap > max(unified_omatrix[i][j], float(options.overlap_fraction)*maxarea):
                       unified_omatrix[i][j] = overlap          
            except:
                print "Warning: y=", y, " x=", x, " k=", k,\
                      " skipped processing of the missing input matrix:", input_opath
    return unified_omatrix 

# -----------------------------------------------------------------------

def generate_fusion_matrix(overlap_matrix, options):
    solver = scip.solver()# quiet=False)
    matrix_shape = overlap_matrix.shape
    fusion_matrix = numpy.zeros(matrix_shape, dtype = numpy.uint32)
    # Setting the boundary elements, which are label ids
    for i in range(0, matrix_shape[0]):
        fusion_matrix[i][0] = overlap_matrix[i][0]
    for j in range(0, matrix_shape[1]):
        fusion_matrix[0][j] = overlap_matrix[0][j]

    # Define SCIP variables
    # Map index k of each SCIP variable to the coordinates of nonzero elements
    # of an overlap matrix
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
    # For higher efficiency, consider only those elements of Fusion matrix
    # which correspond to nonzero elements of an overlap matrix
    # 1) sum of elements of fusion_matrix in each row is <= 2
    for r in range(1, matrix_shape[0]):
        nz_linear_inds = []
        for c in range(1, matrix_shape[1]):
            if overlap_matrix[r][c] == 0:
                continue
            for key, value in dict_map.iteritems():
                if [r, c] == value:
                    nz_linear_inds.append(key)
        if len(nz_linear_inds) > 0:
            solver.constraint( sum(variables[k] for k in nz_linear_inds) <= 2 )# allow binary branching
            if options.verbose and options.debug:
                print "r=", r, " nz_linear_inds=", nz_linear_inds

    # 2) sum of elements of fusion_matrix in each column is <= 2
    for c in range(1, matrix_shape[1]):
        nz_linear_inds = []
        for r in range(1, matrix_shape[0]):
            if overlap_matrix[r][c] == 0:
                continue
            for key, value in dict_map.iteritems():
                if [r, c] == value:
                    nz_linear_inds.append(key)
        if len(nz_linear_inds) > 0:
            solver.constraint( sum(variables[k] for k in nz_linear_inds) <= 2 ) # allow binary branching
            if options.verbose and options.debug:
                print "c=", c, " nz_linear_inds=", nz_linear_inds

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

def convert_omatrix_to_sparse_format(matrix):
    sparse_matrix = [[matrix.shape[0], matrix.shape[1], 0]]                           
    for i in range(1,matrix.shape[0]):
        for j in range(1, matrix.shape[1]):
            if matrix[i][j] > 0: 
                sparse_matrix.append([matrix[i,0], matrix[0,j], matrix[i,j]])
    return numpy.array(sparse_matrix)

# -----------------------------------------------------------------------

def convert_fmatrix_to_sparse_format(matrix):
    sparse_matrix = [[matrix.shape[0], matrix.shape[1], 0]]
    for i in range(1,matrix.shape[0]):
        if sum(matrix[i,1:]) == 0: # region from upper layer does not map to anywhere
            sparse_matrix.append([matrix[i,0], -1, 0])
        for j in range(1, matrix.shape[1]):
            if matrix[i][j] > 0:
                sparse_matrix.append([matrix[i,0], matrix[0,j], matrix[i,j]])
    for j in range(1, matrix.shape[1]):
        if sum(matrix[1:,j]) == 0: # region from lower layer is not mapped to from anywhere
            sparse_matrix.append([-1, matrix[0,j], 0])
    return numpy.array(sparse_matrix)

# -----------------------------------------------------------------------

if __name__ == "__main__":

    usage = "\nUsage: %prog input_data input_type [options (-h to list)]\n"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = MS_LIB_Options.MergeMatrices_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if len(args) == 2:
        input_data, input_type = args
        input_label    = "ms3_DVID"
        if input_type == "directory":
            input_label = "ms3_" + input_data[4:]
        else:
            input_label = "ms3_" + input_data.split('.')[0][4:]

        unified_omatrix = produce_unified_omatrix(input_data,input_type,options)
        output_opath = os.path.join(ms_temp, input_label + "_overlap" +\
                                    "_z" + str(options.node) + ".txt")
        sparse_omatrix = convert_omatrix_to_sparse_format(unified_omatrix)
        numpy.savetxt(output_opath, sparse_omatrix, fmt='%10u', delimiter='\t')

        unified_fmatrix = generate_fusion_matrix(unified_omatrix, options)
        if options.debug:
            check_fusion_matrix(unified_fmatrix)
        sparse_fmatrix = convert_fmatrix_to_sparse_format(unified_fmatrix)
        if options.debug:
            print "sparse_fmatrix=", sparse_fmatrix
        output_fpath = os.path.join(ms_temp, input_label + "_fusion" +\
                                    "_z" + str(options.node) + ".txt")

        if options.verbose:
            print "Saving the sparse unified fusion matrix in file ", output_fpath
        numpy.savetxt(output_fpath, sparse_fmatrix, fmt='%10u', delimiter='\t')
 
    else:
        print usage
        sys.exit(2)


