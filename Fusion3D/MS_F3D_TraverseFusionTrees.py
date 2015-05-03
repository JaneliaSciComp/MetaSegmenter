#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#
import os, sys, re, optparse
import numpy

import MS_LIB_Options

ms_home = os.environ['MS_HOME']
ms_data = os.environ['MS_DATA']
ms_temp = os.environ['MS_TEMP']

sys.setrecursionlimit(2000)
recursion_depth = 0
labels  = []  # global variable: 1st index = layer; 2nd = region in the layer; value = final label to be assigned to the region

# -----------------------------------------------------------------------

# Finalu/new labels are stored in the list labels
# Trees have been employed for convenience of traversal
class Tree(object):
    def __init__(self, position, label):
        self.position              = position
        self.label                 = label
        self.children              = []         

    def traverse(self, dmap, umap, options):
        global recursion_depth
        global labels
        z = int(self.position[0]) # layer  id
        r = int(self.position[1]) # region id
        recursion_depth += 1
   
        # Generate children at lower layer
        if z < int(options.zmax)-1 and r in dmap[z].keys():
            self.children = dmap[z][r]
            if len(self.children) > 0:
                for c in self.children:
                    if c not in labels[z+1].keys():           
                        labels[z+1][c] = self.label
                        print "z=",z," r=",r," dmap=",dmap[z][r]," c=",c,\
                              " labels[",z+1,"][",c,"] now ",self.label
                        t = Tree([z+1,c], self.label)
                        t.traverse(dmap, umap, options) 

        # Genarate children at upper layer
        if z > 0 and r in umap[z].keys():
            self.children = umap[z][r]
            if len(self.children) > 0:
                for c in self.children:
                    if c not in labels[z-1].keys():         
                        labels[z-1][c] = self.label
                        print "z=",z," r=",r," umap=",umap[z][r]," c=",c,\
                              " labels[",z-1,"][",c,"] now ",self.label
                        t = Tree([z-1,c], self.label)
                        t.traverse(dmap, umap, options)

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

def initialize_labels_and_generate_maps(input_label, options):  

    global labels

    # Generate maps
    dmap   = [  ]     # list of "down map" hashes, or the maps of a region in a given layer to region(s) of adjacent lower layer 
    umap   = [{}]     # list of "up   map" hashes, or the maps of a region in a given layer to region(s) of adjacent upper layer 
    fmatrix= [  ]     # list of fusion matrices produced by the previous step of processing ( by MS_F3D_MergeMatrices.py )
    k = 0
    for z in range(int(options.zmin), int(options.zmax)):
        labels.append({})
        dmap.append({}) # initializing a hash at layer k
        umap.append({}) # initializing a hash at layer k
        if z < int(options.zmax)-1:
            input_fpath = os.path.join(ms_temp, input_label + "_fusion" +\
                                       "_z" + str(k+1) + ".txt")
            if options.verbose:
                print "Loading unified fusion matrix ", input_fpath
            fmatrix.append(numpy.loadtxt(input_fpath, delimiter='\t'))
            num_rows = int(fmatrix[k][0][0])-1 # 1st row in sparse fmatrix indicates shape of original fmatrix
            num_cols = int(fmatrix[k][0][1])-1
            dmap.append({})
            umap.append({})
            for r in range(1, num_rows+1):
                dmap[k  ][r] = []
            for c in range(1, num_cols+1):
                umap[k+1][c] = []
            my_shape = fmatrix[k].shape
            for i in range(1, my_shape[0]):
                r = int(fmatrix[k][i][0])
                c = int(fmatrix[k][i][1])
                dmap[k  ][r].append(c)             
                umap[k+1][c].append(r)
        else:
            # Process the last layer
            dmap.append({})
        k = k + 1
        sys.stdout.flush()
    return (dmap, umap)

# -----------------------------------------------------------------------

def produce_fusion_trees(dmap, umap, options):
    global labels
    global recursion_depth
    num_trees = 0;
    print "options.zmax=", options.zmax, " options.zmin=", options.zmin
    current_label = 1
    for z in range(int(options.zmin), int(options.zmax)):
        # Generate trees with root nodes at layer z
        if z < int(options.zmax)-1:
            for r in dmap[z].keys():       # r is the current label    
                if not r in labels[z].keys():         
                    print "Processing label ", r, " at layer ", z , " labels[z][r] was None, became=", r
                    labels[z][r] = current_label
                    recursion_depth = 0
                    t = Tree([z,r], current_label)     # start new tree
                    t.traverse(dmap, umap, options)
                    num_trees += 1
                    print "\n"
                    print "After traversal of tree #", num_trees, " starting at layer=", z, " recursion_depth=", recursion_depth
                    print "\n*************************************************"
                    current_label = current_label + 1
        else:
            for r in umap[z].keys():
                if not r in labels[z].keys():         
                    print "Processing label ", r, " at layer ", z , " labels[z][r] was None, became=", r
                    labels[z][r] = current_label
                    recursion_depth = 0
                    t = Tree([z,r], current_label)     # start new tree
                    t.traverse(dmap, umap, options)
                    num_trees += 1
                    print "\n"
                    print "After traversal of tree #", num_trees, " starting at layer=", z, " recursion_depth=", recursion_depth
                    print "\n*************************************************"
                    current_label = current_label + 1
    return num_trees

# -----------------------------------------------------------------------

def output_new_labels(input_label, options):
    global labels
    labels_list = []
    for z in range(0, len(labels)):
        keys = sorted(labels[z].keys())
        for i in range(0, len(keys)):
            k = keys[i]
            labels_list.append([z, k, labels[z][k]])
    output_path = os.path.join(ms_temp, input_label + "_labels.txt")
    numpy.savetxt(output_path, labels_list, fmt='%10u', delimiter='\t')
 
# -----------------------------------------------------------------------

if __name__ == "__main__":

    usage = "\nUsage: %prog input_data input_type [options (-h to list)]\n"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = MS_LIB_Options.TraverseFusionTrees_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if len(args) == 2:
        input_data, input_type = args
        input_label    = "ms3_DVID"
        if input_type == "directory":
            input_label = "ms3_" + input_data[4:]
        else:
            input_label = "ms3_" + input_data.split('.')[0][4:]
        if options.verbose:
            print "Generating maps ..."
        dmap, umap = initialize_labels_and_generate_maps(input_label, options)
        if options.verbose:
            print "Generating fusion trees ..."
        num_trees = produce_fusion_trees(dmap, umap, options)      
        output_new_labels(input_label, options)
    else:
        print usage
        sys.exit(2)
        

