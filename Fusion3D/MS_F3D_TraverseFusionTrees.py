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

recursion_depth = 0
labels  = []  # global variable: 1st index = layer; 2nd = region; value = new label

# -----------------------------------------------------------------------

# Finalu/new labels are stored in the list labels
# Trees have been employed for convenience of traversal
class Tree(object):
    def __init__(self, position, root_node, parent, label):
        self.root_node             = root_node
        self.position              = position
        self.label                 = label
        self.parent                = parent       
        self.children              = []         

    def traverse(self, dmap, umap, options):
        global recursion_depth
        global labels
        k = self.position[0] # layer  id
        r = self.position[1] # region id
        recursion_depth += 1
#       print "recursion_depth=", recursion_depth, " current label=", self.label
   
        # Generate children at lower layer
        if k < int(options.zmax)-1 and r in dmap[k].keys():
            self.children = dmap[k][r]
#           print "self.children lower layer=", self.children
            if not hasattr(self.children , '__array__'):
                self.children = [ self.children ]
            for c in self.children:
                if labels[k+1][c] == None:
                    print "labels[", k+1, "][", c, "] was ", labels[k+1][c], " now ", self.label
                    labels[k+1][c] = self.label
#                   print "k=", k, " r=", r, " picking child ", c, " from lower layer, ", k+1
                    t = Tree([k+1,c], self.root_node, self, self.label)
                    t.traverse(dmap, umap, options) 

        # Genarate children at upper layer
        if k > 0 and r in umap[k].keys():
            self.children = umap[k][r]
#           print "self.children upper layer=", self.children
            if not hasattr(self.children , '__array__'):
                self.children = [ self.children ]
            for c in self.children:
#               print "k=", k, " c=", c, " labels[k-1][c]=", labels[k-1][c]
                if labels[k-1][c] == None:
                    print "labels[", k-1, "][", c, "] was ", labels[k-1][c], " now ", self.label
                    labels[k-1][c] = self.label
#                   print "k=", k, " r=", r, " picking child ", c, " from upper layer, ", k-1
                    t = Tree([k-1,c], self.root_node, self, self.label)
                    t.traverse(dmap, umap, options)

# -----------------------------------------------------------------------

def initialize_labels_and_generate_maps(input_label, options):  

    global labels

    # Generate maps
    dmap   = [  ]
    umap   = [{}] 
    fmatrix= [  ]
    k = 0
    for z in range(int(options.zmin), int(options.zmax)):
        labels.append({})
        dmap.append({}) # initializing a hash at layer k
        umap.append({}) # initializing a hash at layer k
        if z < int(options.zmax):
            input_fpath = os.path.join(ms_temp, input_label + "_fusion" +\
                                       "_z" + str(k+1) + ".txt")
            if options.verbose:
                print "Loading unified fusion matrix ", input_fpath
            fmatrix.append(numpy.loadtxt(input_fpath, delimiter='\t'))
            print "Labels=", fmatrix[k][1:,0]
            my_shape = fmatrix[k].shape
            for i in range(1, my_shape[0]):
                # Initialize labels at layer k
                labels[k][i] = None
                dmap[k  ][i] = []
                for j in range(1, my_shape[1]):
                    if fmatrix[k][i][j] == 1:
                        dmap[k  ][i] =  j # map from current to next    layer 
            for j in range(1, my_shape[1]):
                # Initialize labels at layer k
                umap[k+1][j] = []
                for i in range(1, my_shape[0]):
                    if fmatrix[k][i][j] == 1:
                        umap[k+1][j] =  i # map from next    to current layer
        else:
            # Process the last layer
            dmap.append({})
            for j in range(1, my_shape[1]):
                # Initialize labels at last layer
                labels[k][j] = None
        k = k + 1
    return (dmap, umap, fmatrix)

# -----------------------------------------------------------------------

def produce_fusion_trees(dmap, umap, fmatrix, options):
    global labels
    global recursion_depth
    trees = []
    k = 0  # layer id
    print "options.zmax=", options.zmax, " options.zmin=", options.zmin
    for z in range(int(options.zmin), int(options.zmax)):
        # Generate trees with root nodes at layer z
        print "int(max(fmatrix[k][1:][0]))=", int(max(fmatrix[k][1:][0]))
        print " In produce fusion trees Labels=", fmatrix[k][1:,0]
        for r in range(1, int(max(fmatrix[k][1:,0]))):         
            if labels[k][r] == None:
                print "Processing label ", r, " at layer ", k , " labels[k][r] was=", labels[k][r], " became=", r
                labels[k][r] = r
                recursion_depth = 0
                t = Tree([k,r], [k,r], None, r)
                t.traverse(dmap, umap, options)
        k = k + 1
    return trees

# -----------------------------------------------------------------------

def output_new_labels(input_label, options):
    global labels
    print "len(labels)=", len(labels)
    for k in range(0, len(labels)):
        labels_list = [ 0 ]
        for i in range(1, max(labels[k].keys())):
            labels_list.append(int(labels[k][i]))
        output_path = os.path.join(ms_temp, input_label + "_labels" +\
                                       "_z" + str(k+1) + ".txt")
        print "k=", k, " saving labels ",  labels_list
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
        dmap, umap, fmatrix = initialize_labels_and_generate_maps(input_label, options)
        if options.verbose:
            print "Generating fusion trees ..."
        trees      = produce_fusion_trees(dmap, umap, fmatrix, options)      
        output_new_labels(input_label, options)
    else:
        print usage
        sys.exit(2)
        

