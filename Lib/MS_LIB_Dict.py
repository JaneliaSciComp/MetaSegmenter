# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

import os, sys

# --------------------------------------------------------------------

# Given node id, determine z-layer id, ymin, ymax, xmin and xmax

def map_node_to_xyz(input_dim, input_label, output_suffix, options):
    z_min   = max(0,                    int(options.zmin))
    z_max   = min(z_min + input_dim[2], int(options.zmax))
    z_range = range(z_min, z_max)
    y_range = range(0, int(options.ny))
    x_range = range(0, int(options.nx))
    if options.verbose:
        print "options.zmin=", options.zmin
        print "input_dim[2]=", input_dim[2], "options.zmax=", options.zmax
        print "x_range=", x_range, " y_range=", y_range, " z_range=", z_range
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
            if int(options.ny) > 1:
                ymax =   y_size + int(options.dy)
        elif y > 0 and y < max(y_range)-1:
            ymin =  y   *y_size - int(options.dy)
            ymax = (y+1)*y_size + int(options.dy)
        else:
            ymin =  y   *y_size - int(options.dy)
            ymax =  input_dim[0]
        for x in x_range:
            if x == 0:
                xmin =  0
                xmax =       x_size
                if int(options.nx) > 1:
                    xmax =   x_size + int(options.dx)
            elif x > 0 and x < max(x_range)-1:
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
                    if options.verbose:
                            print "node=",node
                else: # process only previously unprocessed data
                    if int(options.nx) > 1 and int(options.ny) > 1:
                        output_file = \
                            os.path.join(ms_temp, input_label + "_y" + str(y+1) +\
                                              "_x" + str(x+1) + "_z" + str(z+1) +\
                                              output_suffix)
                    elif int(options.nx) > 1 and int(options.ny) == 1:
                        output_file = \
                            os.path.join(ms_temp, input_label                   +\
                                              "_x" + str(x+1) + "_z" + str(z+1) +\
                                              output_suffix)
                    elif int(options.nx) == 1 and int(options.ny) > 1:
                        output_file = \
                            os.path.join(ms_temp, input_label + "_y" + str(y+1) +\
                                                                "_z" + str(z+1) +\
                                              output_suffix)
                    else:
                        output_file = \
                            os.path.join(ms_temp, input_label                   +\
                                                                "_z" + str(z+1) +\
                                              output_suffix)
                    if not os.path.exists(output_file):
                        node += 1
                        dict_node_xyz[node] = \
                            [y, ymin, ymax, x, xmin, xmax, z]
                        if options.verbose:
                            print "node=",node," missing output file=",output_file
    return (node, dict_node_xyz)

# ----------------------------------------------------------------------

def map_nodes_xmerge(xdim, ydim, zdim, options):
    z_range = range(max(0,   int(options.zmin)),\
                    min(zdim,int(options.zmax)))
    y_range = range(0, int(options.ny))

    x_size = int(round(float(xdim)/float(options.nx)))
    y_size = int(round(float(ydim)/float(options.ny)))
    if x_size <= 2*int(options.dx):
        print "\nx overlap ", options.dx, " exceeds half of x_size ", x_size
        sys.exit(2)
    dict_nodes_xmerge = {}
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
            ymax =  ydim
        for z in z_range:
            node += 1
            dict_nodes_xmerge[node] = [y, ymin, ymax, z]
    return (node, dict_nodes_xmerge)

# ----------------------------------------------------------------------

def map_nodes_ymerge(ydim, zdim, options):
    z_range = range(max(0,   int(options.zmin)),\
                    min(zdim,int(options.zmax)))

    y_size = int(round(float(ydim)/float(options.ny)))
    if y_size <= 2*int(options.dy):
        print "\ny overlap ", options.dy, " exceeds half of y_size ", y_size
        sys.exit()
    dict_nodes_ymerge = {}
    node = 0
    for z in z_range:
        node += 1
        dict_nodes_ymerge[node] = z
    return (node, dict_nodes_ymerge)

# ----------------------------------------------------------------------

def map_nodes_zmerge(zdim, options):

    dict_nodes_zmerge = {}
    dict_nodes_zmerge[1] = "BW.png"
    dict_nodes_zmerge[2] = "Seg.png"
#   dict_nodes_zmerge[3] = "RGB.png"
    return (3, dict_nodes_zmerge)


