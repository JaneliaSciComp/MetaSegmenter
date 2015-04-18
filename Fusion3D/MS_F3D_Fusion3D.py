#! /usr/local/python-2.7.6/bin/python     
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

# Purpose: submit 3D fusion processing steps to cluster,
#          then post-process results

# ------------------------- imports -------------------------

import os, glob, gc
import shutil, commands
import sys, re, optparse
import numpy
import h5py
from scipy import misc
import imghdr
import matplotlib
from PIL import Image
import ntpath

ms_home = os.environ['MS_HOME'] # where source code is located
ms_data = os.environ['MS_DATA'] # dir for initial input data and final output 
ms_temp = os.environ['MS_TEMP'] # dir for intermediate data files

# ----------------------------------------------------------------------

def fusion3d_command_line_parser(parser):
    parser.add_option("-A", "--project",dest="project_code",help="code to be used with qsub",metavar="project_code", default="flyTEM")
    parser.add_option("-D", "--debug",dest="debug",help="don't delete intermediate outputs", action="store_true", default=False)
    parser.add_option("-e", "--executable",dest="executable",help="executable",metavar="executable",default="MS_S2D_Segmentation2D")
    parser.add_option("-l", "--nlen",   dest="nlen", help="# of subsections for processing a fragment in y (length) direction",metavar="nlen", default="1")
    parser.add_option("-m", "--maxsize",dest="msize", help="# of subsections for processing a fragment in y (length) direction",metavar="nlen",default=sys.maxint)
    parser.add_option("-n", "--node",   dest="node", help="id of the cluster node to be used", metavar="node",  default=0)
    parser.add_option("-o", "--output_folder",dest="output_folder",help="output folder",metavar="output_folder",default=ms_data)
    parser.add_option("-p", "--processing_step_beg",dest="processing_step_beg",help="start processing from matr(1), merge(2), trav(3), relab(4) or epil(5)",metavar="processing_step_beg",default=1)
    parser.add_option("-P", "--processing_step_end",dest="processing_step_end",help="end processing with matr(1), merge(2), trav(3), relab(4) or epil(5)",metavar="processing_step_end",default=5)
    parser.add_option("-s", "--sub",dest="submission_command", help="source, qsub or qsub_debug", metavar="submission_command", default="qsub")
    parser.add_option("-S", "--slots",dest="num_slots", help="# of cluster slots per 1 job (default=1)", metavar="num_slots", default=1)
    parser.add_option("-v", "--verbose",action="store_true",dest="verbose",help="increase the verbosity level of output",default=False)
    parser.add_option("-w", "--nwid",   dest="nwid", help="# of subsections for processing a fragment in x (width) direction", metavar="nwid",  default="1")
    parser.add_option("-X", "--nx",  dest="nx",  help="# of image fragments in x direction", metavar="nx", default=1)
    parser.add_option("-Y", "--ny",  dest="ny",  help="# of image fragments in y direction", metavar="ny", default=1)
    parser.add_option("-x", "--dx",  dest="dx",  help="# of scans for fragment overlap in x direction", metavar="dx", default=50)
    parser.add_option("-y", "--dy",  dest="dy",  help="# of scans for fragment overlap in y direction", metavar="dy", default=50)
    parser.add_option("-z", "--zmin",dest="zmin",help="min z-layer to be processed", metavar="zmin", default=0)
    parser.add_option("-Z", "--zmax",dest="zmax",help="max z-layer to be processed", metavar="zmax", default=sys.maxint)
    
    return parser

# -----------------------------------------------------------------------

def is_int(my_sting):     
    try: 
        int(my_sting)
        return True
    except ValueError:
        return False

# -----------------------------------------------------------------------------

def is_DVID_subvolume(input_data, options):
    is_subvolume = True
    if not input_data[0] == "[" or not input_data[len(input_data)-1] == "]":
        is_subvolume = False
        print "Input does not start with '[' or does not end with ']'"
        return is_subvolume
    ranges = input_data[1:len(input_data)-1].split(",")
    if not len(ranges) == 3:
        is_subvolume = False
        print "Input does not contain 3 ranges"
        return is_subvolume
    for range in ranges:
        coords = range.split(":")
        if not len(coords) == 2:
            is_subvolume = False
            print "Input coordinates", coords, " are missing ':'"
            return is_subvolume
        if not is_int(coords[0]) or \
           not is_int(coords[1]) or \
           int(coords[0]) > int(coords[1]):
            is_subvolume = False
            print "Input coordinates", coords, " are in bad order"
            return is_subvolume
    return is_subvolume

# -----------------------------------------------------------------------------

def get_input_type(input_data, options):
    input_type = ""
    if   os.path.isfile(input_data) or \
         os.path.isfile(os.path.join(ms_data,input_data)) or \
         os.path.isfile(os.path.join(ms_temp,input_data)):
        input_type = "file"
    elif os.path.isdir(input_data) or \
         os.path.isdir(os.path.join(ms_data,input_data)) or \
         os.path.isdir(os.path.join(ms_temp,input_data)):
        input_type = "directory"
    elif is_DVID_subvolume(input_data, options):
        input_type = "DVID"
    return input_type

# ----------------------------------------------------------------------

def get_data_dimensions(input_data, input_type, options):
    input_dim = [0,0,0]
    if input_type == "file":
        try:
            f  = h5py.File(os.path.join(ms_data,input_data), 'r')
            key = f.keys()[0]
            data = numpy.transpose(f[key])
            input_dim = data.shape[0:3]
        except:
            print "Input file",input_data,"is not an HDF5 stack\n"
            sys.exit(2)
    elif input_type == "directory":
        files = os.listdir(os.path.join(ms_data,input_data))
        num_files = 0
        image_size = []
        for ifile in files:
            myfile = os.path.join(os.path.join(ms_data,input_data), ifile)
            try:
                image_shape = Image.open(myfile).size
            except:
                myfile = os.path.join(os.path.join(ms_temp,input_data), ifile)
            try:
                suffix3 = myfile[(len(myfile)-3):len(myfile)]
                suffix4 = myfile[(len(myfile)-4):len(myfile)]
#               print "    suffix3=", suffix3, " suffix4=", suffix4
#               print "suffix3 in ['png','tif']:", suffix3 in ['png','tif']
                if  suffix3 in ['png', 'tif' ] or \
                    suffix4 in ['tiff','jpeg'] or \
                    imghdr.what(myfile) in ['png','tiff','jpeg']:
                    num_files = num_files + 1
                    if num_files == 1:
#                       print "Reading first file", myfile," ..."
#                       image_shape = misc.imread(myfile).shape
                        try:
                            # Using PIL
                            image_shape = Image.open(myfile).size           
                        except:
                            try:
                                # Using scipy
                                image_shape = misc.imread(myfile).shape
                            except:
                                try:
                                    # Using Matplotlib
                                    image_shape = matplotlib.image.imread(myfile, format = 'png').shape
                                except:
                                    print "Unable to read image file", myfile
                                    sys.exit(2)
                    input_dim[0:2] = image_shape[0:2]
            except:
                continue
        input_dim[2] = min(num_files, int(options.zmax)-int(options.zmin))
    elif input_type == "DVID":
        ranges = input_data[1:len(input_data)-1].split(",")
        for i in range(0,3):
            coords = ranges[i].split(":")
            input_dim[i] = int(coords[1]) - int(coords[0]) 
    return input_dim

# ----------------------------------------------------------------------

# Given node id, determine z-layer id, ymin, ymax, xmin and xmax

def map_node_to_xyz(input_dim, options):
    z_range = range(max(0,           int(options.zmin)),\
                    min(input_dim[2],int(options.zmax)))
    y_range = range(0, int(options.ny))
    x_range = range(0, int(options.nx))

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
            ymax =       y_size + int(options.dy)
        elif y > 0 and y < max(y_range):
            ymin =  y   *y_size - int(options.dy)
            ymax = (y+1)*y_size + int(options.dy) 
        else:
            ymin =  y   *y_size - int(options.dy)
            ymax =  input_dim[0]
        for x in x_range:
            if x == 0:
                xmin =  0
                xmax =       x_size + int(options.dx)
            elif x > 0 and x < max(x_range):
                xmin =  x   *x_size - int(options.dx)
                xmax = (x+1)*x_size + int(options.dx)
            else:
                xmin =  x   *x_size - int(options.dx)
                xmax =  input_dim[1]
            for z in z_range:  
                node += 1
                dict_node_xyz[node] = \
                    [y, ymin, ymax, x, xmin, xmax, z]
    return (node, dict_node_xyz)

# ----------------------------------------------------------------------

# Genarate overgal matrix between segmented fragments   

def create_generate_matrices_job(outfolderpath, \
                            input_data, input_type, input_dim, \
                            num_nodes, input_label, options):
    generate_matrices_script_path = \
        os.path.join(outfolderpath, "GenerateMatrices_script_ms3." +\
                     input_label + ".sh")
    executable_matrices = os.path.join(ms_home, "Fusion3D",\
                                       "MS_F3D_GenerateMatrices.py")
    if options.verbose:
        print "Number of generate matrices nodes=", num_nodes
    scr = open(generate_matrices_script_path, 'wt')
    scr.write("#!/usr/bash\n")
    scr.write("#$ -t 1-" + str(num_nodes) + "\n")
    command_matrices=    "python " + executable_matrices + \
        "    '"+ input_data + "' " + " file -n $SGE_TASK_ID " + \
        " -X " + str(options.nx)   + " -x " + str(options.dx) + \
        " -Y " + str(options.ny)   + " -y " + str(options.dy) + \
        " -z " + str(options.zmin) + " -Z " + str(options.zmax) +\
        " -L " + str(input_dim[0]) + " -W " + str(input_dim[1]) +\
        " -H " + str(input_dim[2])
    if options.verbose:
        command_matrices += " -v "
    if options.debug:
        command_matrices += " -D "
    command_matrices    += "\n"
    scr.write(command_matrices)   
    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", \
                  generate_matrices_script_path]))
    scr.write("\n")
    scr.close()
    os.chmod(generate_matrices_script_path, 0777)
    return generate_matrices_script_path

# -----------------------------------------------------------------------------

def create_merge_job(outfolderpath, input_data, input_type,\
                       input_dim, input_label, options):
    zmin = max(0,            int(options.zmin))
    zmax = min(input_dim[2], int(options.zmax))
    num_nodes = zmax - zmin
    if options.verbose:
        print "Number of merge nodes=", num_nodes

    merge_script_path = \
        os.path.join(outfolderpath, "Merge_fusion_script_ms3." +\
                     input_label + ".sh")
    executable_path = os.path.join(ms_home,"Fusion3D", \
                                   "MS_F3D_MergeFusionMatrices.py")
    scr = open(merge_script_path, 'wt')
    scr.write("#!/usr/bash\n")
    scr.write("#$ -t 1-" + str(num_nodes) + "\n")
    command_merge = executable_path \
                      + " " + input_data + " " + input_type \
                      + " -n " + " $SGE_TASK_ID " \
                      + " -X " + str(options.nx)  \
                      + " -Y " + str(options.ny)  \
                      + " -z " + str(zmin)        \
                      + " -Z " + str(zmax)           
    if options.verbose:
        command_merge += " -v "
    if options.debug:
        command_merge += " -D"
    command_merge += "\n"
    scr.write(command_merge) 
    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", \
                  merge_script_path]))
    scr.close()
    os.chmod(merge_script_path, 0777)
    return merge_script_path

# -----------------------------------------------------------------------------

def create_traverse_job(outfolderpath, input_data, input_type,\
                      input_dim, input_label, options):
    zmin = max(0,            int(options.zmin))
    zmax = min(input_dim[2], int(options.zmax))
    if options.verbose:
        print "Number of traverse  nodes=", 1
    traverse_script_path = \
        os.path.join(outfolderpath, "Traverse_script_ms3." +\
                     input_label + ".sh")
    executable_path = os.path.join(ms_home,"Fusion3D", \
                                   "MS_F3D_TraverseFusionTrees.py")
    scr = open(traverse_script_path,'wt')
    scr.write("#!/usr/bash\n")
    command_traverse = executable_path \
                   + "    " + input_data + "   "  + input_type \
                   + " -z " + str(zmin)  + " -Z " + str(zmax)
    if options.verbose:
        command_traverse += " -v "
    if options.debug:
        command_traverse += " -D"
    command_traverse += "\n"
    scr.write(command_traverse)
    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", \
                  traverse_script_path]))
    scr.close()
    os.chmod(traverse_script_path, 0777)
    return traverse_script_path

# -----------------------------------------------------------------------------

def create_relabel_job(outfolderpath, input_data, input_type,\
                       input_dim, input_label, options):
    relabel_script_path = \
        os.path.join(outfolderpath, "Relabel_script_ms3." +\
                     input_label + ".sh")
    executable_path = os.path.join(ms_home,"Fusion3D", \
                                   "MS_F3D_RelabelSegmentedData.py")
    zmin = max(0,            int(options.zmin))
    zmax = min(input_dim[2], int(options.zmax))
    num_nodes = zmax - zmin
    scr = open(relabel_script_path,'wt')
    scr.write("#!/usr/bash\n")
    scr.write("#$ -t 1-" + str(num_nodes) + "\n")
    command_relabel = executable_path \
                   + " -n " + " $SGE_TASK_ID " \
                   + "    " + input_data + "    " + input_type \
                   + " -z " + str(zmin)  + " -Z " + str(zmax)
    if options.verbose:
        command_relabel += " -v "
    if options.debug:
        command_relabel += " -D"
    command_relabel += "\n"
    scr.write(command_relabel)
    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", \
                  relabel_script_path]))
    scr.close()
    os.chmod(relabel_script_path, 0777)
    return relabel_script_path

 # -----------------------------------------------------------------------------

def create_epilog_job(outfolderpath, input_dim, input_label, options):
    epilog_script_path = \
        os.path.join(ms_temp, "Epilog_script_" +\
                     input_label + ".sh")
    executable_path = os.path.join(ms_home, "Utilities", \
                                   "MS_UT_CreateH5Stack.py")
    scr = open(epilog_script_path, 'wt')
    print "\ninput_label=", input_label
    output_path = os.path.join(ms_data, input_label + ".h5")
    scr.write("#!/usr/bash\n")
    command_epilog = executable_path    + " " + ms_temp + " -t data " + \
                     " -m _relabeled -u _y -o " + output_path
    if options.verbose:
        print "command_epilog=", command_epilog
    if options.verbose:
        command_epilog += " -v "
    command_epilog += "\n"
    scr.write(command_epilog)
    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", \
                  epilog_script_path]))
    scr.close()
    os.chmod(epilog_script_path, 0777)
    return epilog_script_path

# -----------------------------------------------------------------------------

def submit_job(script_path, base_command, jobid_in, options):
    prog_name = ""
    if re.search("Matrices",  script_path):
        prog_name = "ms3.matr" 
    elif re.search("Merge_fusion",  script_path):
        prog_name = "ms3.mergfus"
    elif re.search("Traverse",  script_path):
        prog_name = "ms3.travrs"
    elif re.search("Relabel",  script_path):
        prog_name = "ms3.relabl"
    elif re.search("Epilog", script_path):
        prog_name = "ms3.epilog"
    jobid_out = 0
    if options.submission_command == "qsub":
        command = base_command + " -V -N " + prog_name +\
                                 " -o /dev/null -e /dev/null "
    elif re.search("qsub", options.submission_command):
        command = base_command + " -V -N " + prog_name
    else:
        command = base_command
    if re.search("qsub", options.submission_command):
        if jobid_in > 0:
            command += " -hold_jid " + str(jobid_in)
        command += " " + script_path
        if not re.search("epilog", script_path):
            res = commands.getstatusoutput(command)
            jobid_out = (res[1].split()[2]).split(".")[0]
        else:
            os.system(command)
            jobid_out = 0
    elif not options.submission_command == "none":
        command += " " + script_path
        os.system(command)
    if  options.verbose and \
        not options.submission_command in ["none", "source"]:
        print "\nSubmit job command=", command, " res=", res
        print "jobid_out=", jobid_out
    return jobid_out

# -----------------------------------------------------------------------------

def submit_all_jobs(generate_matrices_script_path,\
                    merge_script_path, traverse_script_path,\
                    relabel_script_path, epilog_script_path, options):
    qsub = os.path.join(os.environ['SGE_ROOT'], "bin", "lx-amd64", "qsub")       
    if re.search("qsub", options.submission_command):
        base_command = qsub
    else:
        base_command = options.submission_command
    if  len(options.project_code) > 0 and \
        re.search("qsub", options.submission_command):
        base_command += "  -A " + options.project_code

    jid1 = 0
    if int(options.processing_step_beg) == 1 and \
       int(options.processing_step_end) >= 1:
        jid1 = submit_job(generate_matrices_script_path, base_command, 0, options)

    jid2 = 0
    if int(options.processing_step_beg) <= 2 and \
       int(options.processing_step_end) >= 2:
        jid2 = submit_job(merge_script_path,base_command,jid1, options)

    jid3 = 0
    if int(options.processing_step_beg) <= 3 and \
       int(options.processing_step_end) >= 3:
        jid3 = submit_job(traverse_script_path,base_command,jid2, options)

    jid4 = 0
    if int(options.processing_step_beg) <= 4 and \
       int(options.processing_step_end) >= 4:
        jid4 = submit_job(relabel_script_path, base_command,jid3, options)

    if int(options.processing_step_beg) <= 5 and \
       int(options.processing_step_end) >= 5:
        submit_job(epilog_script_path, base_command, jid4, options)

    return

# -----------------------------------------------------------------------------

# Create and submit all cluster jobs
def process_inputs(input_data, input_type, input_dim, \
                   num_nodes, input_label, options):
   
    if int(options.processing_step_beg) == 1:
        generate_matrices_script_path = \
            create_generate_matrices_job(ms_temp, \
                                input_data, input_type, input_dim, \
                                num_nodes, input_label, options)
    else:
        generate_matrices_script_path = ""

    if int(options.processing_step_beg) <= 2:
        merge_script_path = \
            create_merge_job(ms_temp, input_data,\
                               input_type, input_dim, input_label, options)
    else:
        merge_script_path = ""

    if int(options.processing_step_beg) <= 3:
        traverse_script_path = \
            create_traverse_job(ms_temp, input_data,\
                               input_type, input_dim, input_label, options)
    else:
        traverse_script_path = ""

    if int(options.processing_step_beg) <= 4:
        relabel_script_path = \
            create_relabel_job(ms_temp, input_data,\
                               input_type, input_dim, input_label, options)
    else:
        relabel_script_path = ""

    epilog_script_path = \
        create_epilog_job(ms_data, input_dim, input_label, \
                          options)
    if options.verbose:
        print "generate_matrices_script_path=", generate_matrices_script_path
        print "merge_script_path=",      merge_script_path
        print "traverse_script_path=",traverse_script_path
        print "relabel_script_path=",  relabel_script_path
        print "epilog_script_path=",    epilog_script_path
    submit_all_jobs(generate_matrices_script_path,\
                      merge_script_path,\
                      traverse_script_path,\
                      relabel_script_path,\
                      epilog_script_path, options)

# -----------------------------------------------------------------------------

if __name__ == "__main__":
   
    usage = "Usage: \n\
    %prog input_data [options (-h to list)] \nwhere\n\
    input_data=\n\
    1) path to image stack file in HDF5 format, or \n\
    2) path to folder containing image files, or \n\
    3) coordinates of DVID subvolume: [xmin:xmax,ymin:ymax,zmin:zmax]"   

#   Steps of processing:
#   1) MS_F3D_GenerateMatrices.py
#   2) MS_F3D_MergeFusionMatrices.py
#   3) MS_F3D_TraverseFusionTrees.py
#   4) MS_F3D_RelabelSegmentedData.py
#   5) produce ouput (in HDF5 or other format)

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = fusion3d_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if len(args) == 1:
        input_data = args[0]
        input_type = get_input_type(input_data, options)
        if input_type == 'DVID':
            input_label = "ms3_DVID"
        elif input_type == "directory":
            input_label = "ms3_" + ntpath.split(input_data)[1][4:]
        else:
            input_label = "ms3_" + ntpath.split(input_data)[1].split('.')[0][4:]
        if input_type in ["file", "directory", "DVID"]:
            if options.verbose and int(options.node) == 0:
                print "Input type=", input_type
        else:
            print "\nIncorrectly specified input data", input_data, "\n"
            parser.print_usage()
            sys.exit(2)         
        input_dim = get_data_dimensions(input_data, input_type, options)
        num_nodes, dict_node_xyz = \
            map_node_to_xyz(input_dim, options)
        if options.verbose:
            print "Input data dimensions: ", input_dim, " num_nodes=", num_nodes
            print "\ninput_label=", input_label   
            print "\nProcessing input data ..."
        process_inputs(input_data, input_type, input_dim, \
                num_nodes, input_label, options)
    else:
        parser.print_usage()
        sys.exit(2)
