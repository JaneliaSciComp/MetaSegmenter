#! /usr/local/python-2.7.6/bin/python     
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#
# File MS_S2D_MergeImageFragments.py
#
# ------------------------- imports -------------------------

import os
import sys, re 
import optparse
import h5py
import numpy

import MS_LIB_Dict
import MS_LIB_Options

ms_home = os.environ['MS_HOME']
ms_data = os.environ['MS_DATA']
ms_temp = os.environ['MS_TEMP']

# -----------------------------------------------------------------------

def compile_code(options):
    print "\nCompiling MS_S2D_MergeImageFragments.m ..."
    input_path = os.path.join(ms_home,"Segmentation2D", \
                              "MS_S2D_MergeImageFragments.m")
    command = "mcc -m -N -p signal -p optim -p stats -p shared -p images " +\
              " -R -singleCompThread " +\
              input_path + " -o " + options.executable + \
              " -d " +  os.path.join(ms_home, "Bin")
    if options.verbose:
        command += " -v "
        print "command= ", command, "\n"
    os.system(command)
    os.system("chmod +x " + os.path.join(ms_home, "Bin", \
                                         options.executable))

# -----------------------------------------------------------------------

def process_inputs_xmerge(input_label, xdim, ydim, node, dict_nodes_xmerge, options):

    node        = int(node)
    y           = dict_nodes_xmerge[node][0]
    ymin        = dict_nodes_xmerge[node][1]
    ymax        = dict_nodes_xmerge[node][2]
    z           = dict_nodes_xmerge[node][3]

    if options.ny > 1:
        output_path = os.path.join(ms_temp, \
            input_label + "_y" + str(y+1) + "_z" + str(z+1) + "_BW.png")
    else:
        output_path = os.path.join(ms_temp, \
            input_label +                   "_z" + str(z+1) + "_BW.png")
    executable_path = os.path.join(ms_home, "Bin", options.executable)
    command_xmerge = executable_path + " x " + str(options.nx)      + " " +\
                     str(options.dx) + "   " + str(options.verbose) + " " +\
                     str(ydim)       + "   " + str(xdim)            + " " +\
                     output_path 
    command_rm     = "rm -f " 
    for x in range(0, int(options.nx)):
        input_path = os.path.join(ms_temp,\
                     input_label + "_y" + str(y+1) + "_x" + str(x+1) + \
                                   "_z" + str(z+1) + "_BW.png")
        command_xmerge += " " + input_path
        command_rm     += " " + input_path
    os.system(command_xmerge)

    if options.verbose:
        print "In MS_S2D_MergeImageFragments.py/process_inputs_xmerge: "
        print "command_xmerge=", command_xmerge
        print "command_rm=", command_rm
 
    if not options.debug:
        os.system(command_rm)

# ----------------------------------------------------------------------

def process_inputs_ymerge(input_label, xdim, ydim, node, dict_nodes_ymerge, options):

    node        = int(node)
    z           = dict_nodes_ymerge[node]

    # Merge fragments in y-direction
    bw_output_path = os.path.join(ms_temp, input_label + "_z" + str(z+1) + "_BW.png")
    merge_executable_path = os.path.join(ms_home, "Bin", options.executable)
    command_ymerge = merge_executable_path  + " y " + str(options.ny)     + " " +\
                     str(options.dy)  + "   " + str(int(options.verbose)) + " " +\
                     str(ydim)        + "   " + str(xdim)                 + " " +\
                     bw_output_path
    command_rm     = "rm -f "
    for y in range(0, int(options.ny)):
        if int(options.ny) > 1:
            input_path = os.path.join(ms_temp,\
                         input_label + "_y" + str(y+1) + "_z" + str(z+1) + \
                         "_BW.png")
        else:
            input_path = os.path.join(ms_temp,\
                         input_label + "_z" + str(z+1) + "_BW.png")
        command_ymerge += " " + input_path
        command_rm     += " " + input_path
    if options.verbose:
        print "In MS_S2D_MergeImageFragments.py/process_inputs_ymerge: "
        print "command_ymerge=", command_ymerge
        print "command_rm=", command_rm

    os.system(command_ymerge)
    if not options.debug:
        os.system(command_rm)

# ----------------------------------------------------------------------

def process_inputs_zmerge(zdim, input_label, node, dict_nodes_zmerge, options):
    node       = int(node)
    match_str0 = dict_nodes_zmerge[node] 
    match_str  = input_label + '*' + match_str0 + "$"

    command_rm     = "rm -f "
    z_range = range(max(0, int(options.zmin)), min(zdim, int(options.zmax)))
    i = 0
    for z in z_range:
        input_file = os.path.join(ms_temp, input_label + "_z" + str(z+1) + "_Seg.h5")
        command_rm += " " + input_file
    
        # Read input file
        if re.search(".h5", input_file):
            fin = h5py.File(input_file, 'r')
            key = fin.keys()[0]
            image_data = numpy.squeeze(fin[key])
        else:
            # Using PIL
            image_data = numpy.asarray(Image.open(input_file))
            print "Using PIL: max_label=", numpy.max(image_data)
        if i == 0:
            stack_shape = [len(z_range), image_data.shape[0], image_data.shape[1]]
            stack_data = numpy.zeros(stack_shape, dtype = numpy.uint64)
        stack_data[i,:,:] = numpy.uint64(image_data[:,:])
        i = i + 1

    output_path = os.path.join(ms_data, \
        input_label + "_" + match_str0.split(".")[0] + ".h5")
    fout  = h5py.File(output_path, 'w')
    fout.create_dataset('main', stack_data.shape, data=stack_data)

    if options.verbose:
        print "In MS_S2D_MergeImageFragments.py/process_inputs_zmerge: "
        print "output_path=", output_path
        print "command_rm=", command_rm

    if not options.debug:
        os.system(command_rm)

# ----------------------------------------------------------------------

def parse_option(option):
    if re.search("-", option) and re.search(",", option):
        sys.exit("Please, use either dash or comma in specifying " + option + "\n")
    if re.search("-", option):
        splt = option.split("-")
        items = range(int(splt[0]), int(splt[1])+1)
    else:
        if re.search(",", option):
            items  = option.split(",")
        else:
            items  = [option]
    return items

# -----------------------------------------------------------------------

if __name__ == "__main__":
   
    usage = "\nUsage: %prog input_label merge_dir xdim ydim zdim node [options (-h to list)]\n"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = MS_LIB_Options.MergeImageFragments_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if options.compile:
        compile_code(options)
    if len(args) == 6:      
        input_label, mdir, xdim, ydim, zdim, node = args[0:6]
        if mdir == 'x':
            num_nodes, dict_nodes_xmerge = MS_LIB_Dict.map_nodes_xmerge(xdim, ydim, zdim, options)
            process_inputs_xmerge(input_label, xdim, ydim, node, dict_nodes_xmerge, options)
        elif mdir == 'y':
            num_nodes, dict_nodes_ymerge = MS_LIB_Dict.map_nodes_ymerge(ydim, zdim, options)
            process_inputs_ymerge(input_label, xdim, ydim, node, dict_nodes_ymerge, options)
        elif mdir == 'z':
            num_nodes, dict_nodes_zmerge = MS_LIB_Dict.map_nodes_zmerge(zdim, options)
            process_inputs_zmerge(zdim, input_label, node, dict_nodes_zmerge, options)
        else:
            print "Unsapported merge direction ", mdir
            sys.exit(2)
    elif options.compile:
        sys.exit(2)
    else:
        print usage
        sys.exit(2)
