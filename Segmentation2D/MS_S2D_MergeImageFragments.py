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

import MS_Dict
import MS_Options

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
    if options.verbose:
        print "command_xmerge=", command_xmerge
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
        print "command_rm=", command_rm
 
    if not options.debug:
        os.system(command_rm)

# ----------------------------------------------------------------------

def process_inputs_ymerge(input_label, xdim, ydim, node, dict_nodes_ymerge, options):

    node        = int(node)
    z           = dict_nodes_ymerge[node]

    output_path = os.path.join(ms_temp, input_label + "_z" + str(z+1) + "_BW.png")
    executable_path = os.path.join(ms_home, "Bin", options.executable)
    command_ymerge = executable_path  + " y " + str(options.ny)           + " " +\
                     str(options.dy)  + "   " + str(int(options.verbose)) + " " +\
                     str(ydim)        + "   " + str(xdim)                 + " " +\
                     output_path
    command_rm     = "rm -f "
    for y in range(0, int(options.ny)):
        input_path = os.path.join(ms_temp,\
                     input_label + "_y" + str(y+1) + "_z" + str(z+1) + \
                     "_BW.png")
        command_ymerge += " " + input_path
        command_rm     += " " + input_path
    os.system(command_ymerge)

    if options.verbose:
        print "In MS_S2D_MergeImageFragments.py/process_inputs_ymerge: "
        print "command_rm=", command_rm

    if not options.debug:
        os.system(command_rm)

# ----------------------------------------------------------------------

def process_inputs_zmerge(zdim, input_label, node, dict_nodes_zmerge, options):
    node      = int(node)
    match_str0 = dict_nodes_zmerge[node] 
    match_str  = input_label + '*' + match_str0 + "$"

    output_path = os.path.join(ms_data, \
        input_label + "_" + match_str0.split(".")[0] + ".h5")
    executable_path = os.path.join(ms_home, "Utilities", \
                                   "MS_UT_CreateH5Stack.py")
    command_zmerge = executable_path    + " " + ms_temp + " -t data " + \
                     " -m " + match_str + " -u _y -o " + output_path +\
                     " -z " + str(options.zmin) + " -Z " + str(options.zmax)
    if options.verbose:
        print "\ncommand_zmerge=", command_zmerge
    command_rm     = "rm -f "
    z_range = range(max(0, int(options.zmin)), min(zdim, int(options.zmax)))
    for z in z_range:
        if match_str0 in ["_Seg.png", "_RGB.png"]:
            input_file = input_label + "_z" + str(z+1) + "_" + match_str0
        else:
            input_file = input_label + "_z" + str(z+1) + "_" + match_str0
    input_path = os.path.join(ms_temp, input_file)
    command_rm     += " " + input_path
    os.system(command_zmerge)

    if options.verbose:
        print "In MS_S2D_MergeImageFragments.py/process_inputs_zmerge: "
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
    parser = MS_Options.MergeImageFragments_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if options.compile:
        compile_code(options)
    if len(args) == 6:      
        input_label, mdir, xdim, ydim, zdim, node = args[0:6]
        if mdir == 'x':
            num_nodes, dict_nodes_xmerge = MS_Dict.map_nodes_xmerge(xdim, ydim, zdim, options)
            process_inputs_xmerge(input_label, xdim, ydim, node, dict_nodes_xmerge, options)
        elif mdir == 'y':
            num_nodes, dict_nodes_ymerge = MS_Dict.map_nodes_ymerge(ydim, zdim, options)
            process_inputs_ymerge(input_label, xdim, ydim, node, dict_nodes_ymerge, options)
        elif mdir == 'z':
            num_nodes, dict_nodes_zmerge = MS_Dict.map_nodes_zmerge(zdim, options)
            process_inputs_zmerge(zdim, input_label, node, dict_nodes_zmerge, options)
        else:
            print "Unsapported merge direction ", mdir
            sys.exit(2)
    elif options.compile:
        sys.exit(2)
    else:
        print usage
        sys.exit(2)
