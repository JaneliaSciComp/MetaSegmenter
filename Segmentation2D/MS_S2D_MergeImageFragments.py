#! /usr/local/python-2.7.8/bin/python     
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#
# File MS_S2D_MergeImageFragments.py
#
# ------------------------- imports -------------------------

import os
import sys, re 
import optparse

ms_home = os.environ['MS_HOME']
ms_data = os.environ['MS_DATA']
ms_temp = os.environ['MS_TEMP']

# -----------------------------------------------------------------------

def merge_command_line_parser(parser):
    parser.add_option("-c", "--compile",action="store_true",dest="compile",help="compile Matlab code",metavar="compile",default=False)
    parser.add_option("-d", "--direction", dest="dir", help="merge direction: x or y", metavar="dir", default="")
    parser.add_option("-D", "--debug",dest="debug",help="debugging mode; don't delete shell scripts",action="store_true",default=False)
    parser.add_option("-e", "--executable", dest="executable", help="name of matlab executable", metavar="executable", default="MS_S2D_MergeImageFragments")
    parser.add_option("-o", "--out_folder", dest="out_folder", help="where to place executable", metavar="out_folder", default=ms_home + "/Bin")
    parser.add_option("-s", "--submission_command", dest="submission_command", help="source, qsub or qsub_debug", metavar="submission_command", default="qsub")
    parser.add_option("-v", "--verbose",action="store_true",dest="verbose",help="increase the verbosity level of output",default=False)
    parser.add_option("-X", "--nx",  dest="nx",  help="# of subsections in x direction", metavar="nx", default=1)
    parser.add_option("-Y", "--ny",  dest="ny",  help="# of subsections in y direction", metavar="ny", default=1)
    parser.add_option("-x", "--dx",  dest="dx",  help="# of scans for image overlap in x direction", metavar="dx", default=50)
    parser.add_option("-y", "--dy",  dest="dy",  help="# of scans for image overlap in y direction", metavar="dy", default=50)
    parser.add_option("-z", "--zmin",dest="zmin",help="# of subsections in z direction", metavar="zmin", default=0)
    parser.add_option("-Z", "--zmax",dest="zmax",help="# of subsections in z direction", metavar="zmax", default=sys.maxint)      
    return parser

# -----------------------------------------------------------------------

def compile_code(options):
    print "\nCompiling MS_S2D_MergeImageFragments.m ..."
    input_path = os.path.join(ms_home,"2DSegmentation", \
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

def process_inputs_xmerge(input_label, node, dict_nodes_xmerge, options):

    node        = int(node)
    y           = dict_nodes_xmerge[node][0]
    ymin        = dict_nodes_xmerge[node][1]
    ymax        = dict_nodes_xmerge[node][2]
    z           = dict_nodes_xmerge[node][3]

    output_path = os.path.join(ms_temp, \
        input_label + "_y" + str(y+1) + "_z" + str(z+1) + "_BW.png")
    executable_path = os.path.join(ms_home, "Bin", options.executable)
    command_xmerge = executable_path + " x " + str(options.nx)      + " " +\
                     str(options.dx) + "   " + str(options.verbose) + " " +\
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

def process_inputs_ymerge(input_label, node, dict_nodes_ymerge, options):

    node        = int(node)
    z           = dict_nodes_ymerge[node]

    output_path = os.path.join(ms_temp, input_label + "_z" + str(z+1) + "_BW.png")
    executable_path = os.path.join(ms_home, "Bin", options.executable)
    command_ymerge = executable_path  + " y " + str(options.ny)      + " " +\
                     str(options.dy)  + "   " + str(int(options.verbose)) + " " +\
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
    match_str  = match_str0 + "$"

    output_path = os.path.join(ms_data, \
        input_label + "_" + match_str0.split(".")[0] + ".h5")
    executable_path = os.path.join(ms_home, "Utilities", \
                                   "MS_UT_CreateH5Stack.py")
    command_zmerge = executable_path    + " " + ms_temp + " -t data " + \
                     " -c -m " + match_str + " -u _y -o " + output_path +\
                     " -z " + str(options.zmin) + " -Z " + str(options.zmax)
    if options.verbose:
        print "\ncommand_zmerge=", command_zmerge
    command_rm     = "rm -f "
    z_range = range(max(0, int(options.zmin)), min(zdim, int(options.zmax)))
    for z in z_range:
        if match_str0 in ["_Seg.png", "_RGB.png"]:
            input_file = input_label + "_z" + str(z+1) + match_str0
        else:
            input_file = input_label + "_z" + str(z+1) + match_str0
    input_path = os.path.join(ms_temp, input_file)
#   command_zmerge += " " + input_path
    command_rm     += " " + input_path
    os.system(command_zmerge)

    if options.verbose:
        print "In MS_S2D_MergeImageFragments.py/process_inputs_zmerge: "
        print "command_rm=", command_rm

    if not options.debug:
        os.system(command_rm)

# ----------------------------------------------------------------------

def map_nodes_xmerge(xdim, ydim, zdim, options):
    print "options.zmin=", options.zmin, " options.zmax=", options.zmax, " zdim=", zdim
    z_range = range(max(0,   int(options.zmin)),\
                    min(zdim,int(options.zmax)))
    print "z_range=", z_range
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
    dict_nodes_zmerge[3] = "RGB.png"
    return (3, dict_nodes_zmerge)

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
   
    usage = "\nUsage: %prog input_label merge_dir ydim xdim zdim node [options (-h to list)]\n"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = merge_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if options.compile:
        compile_code(options)
    if len(args) == 6:      
        input_label, mdir, xdim, ydim, zdim, node = args[0:6]
        if mdir == 'x':
            num_nodes, dict_nodes_xmerge = map_nodes_xmerge(xdim, ydim, zdim, options)
            process_inputs_xmerge(input_label, node, dict_nodes_xmerge, options)
        elif mdir == 'y':
            num_nodes, dict_nodes_ymerge = map_nodes_ymerge(ydim, zdim, options)
            process_inputs_ymerge(input_label, node, dict_nodes_ymerge, options)
        elif mdir == 'z':
            num_nodes, dict_nodes_zmerge = map_nodes_zmerge(zdim, options)
            process_inputs_zmerge(zdim, input_label, node, dict_nodes_zmerge, options)
        else:
            print "Unsapported merge direction ", mdir
            sys.exit(2)
    elif options.compile:
        sys.exit(2)
    else:
        print usage
        sys.exit(2)
