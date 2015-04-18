#! /usr/local/python-2.7.6/bin/python     
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

# Purpose: submit Matlab's 2D segmentation job(s) to cluster,
#          then post-process results

# ------------------------- imports -------------------------

import os, glob, gc
import shutil, commands
import tifffile as tiff
import sys, re, optparse
import numpy
import h5py
from scipy import misc
import matplotlib
from PIL import Image

ms_home = os.environ['MS_HOME'] # where source code is located
ms_data = os.environ['MS_DATA'] # dir for initial input data and final output 
ms_temp = os.environ['MS_TEMP'] # dir for intermediate data files

# ----------------------------------------------------------------------

def segmentation_command_line_parser(parser):
    parser.add_option("-A", "--project",dest="project_code",help="code to be used with qsub",metavar="project_code", default="flyTEM")
    parser.add_option("-c", "--compile",action="store_true",dest="compile",help="compile Matlab file",metavar="compile",default=False)
    parser.add_option("-C", "--compile_all",action="store_true",dest="compile_all", help="compile all Matlab code", metavar="compile_all", default=False)
    parser.add_option("-D", "--debug",dest="debug",help="don't delete intermediate outputs", action="store_true", default=False)
    parser.add_option("-e", "--executable",dest="executable",help="executable",metavar="executable",default="MS_S2D_Segmentation2D")
    parser.add_option("-f", "--fracBlack", dest="fracBlack", help="fracBlack for detecting neurons (default=automatic)",metavar="fracBlack", default=None)
    parser.add_option("-F", "--fracBlack2",dest="fracBlack2",help="fracBlack for detect. dark str.",metavar="fracBlack2",default=None)
    parser.add_option("-l", "--nlen",   dest="nlen", help="# of subsections for processing a fragment in y (length) direction",metavar="nlen", default="1")
    parser.add_option("-m", "--maxsize",dest="msize", help="# of subsections for processing a fragment in y (length) direction",metavar="nlen",default=sys.maxint)
    parser.add_option("-n", "--node",   dest="node", help="id of the cluster node to be used", metavar="node",  default=0)
    parser.add_option("-o", "--output_folder",dest="output_folder",help="output folder",metavar="output_folder",default=ms_data)
    parser.add_option("-p", "--processing_start",dest="processing_start",help="start processing from segm(=1),xmerg(=2),ymerg(=3),or epilog(=4) ",metavar="processing_start",default=1)
    parser.add_option("-P", "--processing_end",dest="processing_end",help="complete processing at step segm(=1),xmerg(=2),ymerg(=3),or epilog(=4) ",metavar="processing_end",default=4)
    parser.add_option("-s", "--sub",dest="submission_command", help="source, qsub or qsub_debug", metavar="submission_command", default="qsub")
    parser.add_option("-S", "--slots",dest="num_slots", help="# of cluster slots per 1 job (default=1)", metavar="num_slots", default=1)
    parser.add_option("-U", "--unprocessed",action="store_true",dest="unprocessed", help="reprocess only the data for which an output file does not exist",default=False)
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

def compile_code(options):
    print "\nCompiling MS_S2D_Segmentation2D.m ..."
    input_path = os.path.join(ms_home,"Segmentation2D", \
                              "MS_S2D_Segmentation2D.m")
    command = "mcc -m -N -p signal -p optim -p stats -p images " + \
              "-p shared -p curvefit -R -singleCompThread " +\
              input_path + " -o " + options.executable + \
              " -d " + os.path.join(ms_home, "Bin")
    if options.verbose:
        command += " -v "
        print "command= ", command, "\n"
    os.system(command)
    os.system("chmod +x " + os.path.join(ms_home, "Bin", \
                                         options.executable))
    if options.compile_all:
        print "\nCompiling MS_S2D_MergeImageFragments.m ..."
        input_path = os.path.join(ms_home,"Segmentation2D", \
                                  "MS_S2D_MergeImageFragments.m")
        command = "mcc -m -N -p signal -p optim -p stats -p shared -p images " +\
                  " -R -singleCompThread " +\
                  input_path + " -o MS_S2D_MergeImageFragments " + \
                  " -d " +  os.path.join(ms_home, "Bin")
        if options.verbose:
            command += " -v "
            print "command= ", command, "\n"
        os.system(command)
        command_chmod = "chmod +x " + os.path.join(ms_home, "Bin", \
                                         "MS_S2D_MergeImageFragments")
        if options.verbose:
            print "command_chmod= ", command_chmod, "\n"
        os.system(command_chmod) 

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
        if re.search(".h5", input_data):
            try:
                f  = h5py.File(os.path.join(ms_data,input_data), 'r')
                key = f.keys()[0]
                image_data = numpy.transpose(f[key])
                input_dim = image_data.shape[0:3]
            except:
                print "Input file",input_data,"is not an HDF5 stack\n"
                sys.exit(2)
        elif re.search(".tif", input_data):
            try:
                image_data = numpy.transpose(tiff.imread(input_data))               
                input_dim = image_data.shape[0:3]
            except:
                print "Input file",input_data,"is not a tiff stack\n"
                sys.exit(2)
        else:
            print "Unsupported input file\n"                
            sys.exit(2)
    elif input_type == "directory":
        files = sorted(os.listdir(os.path.join(ms_data,input_data)))
        num_files = 0
        image_size = []
        for ifile in files:
            myfile = os.path.join(os.path.join(ms_data,input_data), ifile)
            try:
                image_shape = Image.open(myfile).size
            except:
                myfile = os.path.join(os.path.join(ms_temp,input_data), ifile)
            try:
                if re.search(".png",  myfile) or \
                   re.search(".tif",  myfile) or \
                   re.search(".jpeg", myfile):
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
                    input_dim[0:2] = [image_shape[1],image_shape[0]] # transpose
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

def map_node_to_xyz(input_dim, input_label,options):
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
                if not options.unprocessed:
                    node += 1
                    dict_node_xyz[node] = \
                        [y, ymin, ymax, x, xmin, xmax, z]
                else: # process only previously unprocessed data
                    output_file = \
                        os.path.join(ms_temp, input_label + "_y" + str(y+1) +\
                                     "_x" + str(x+1) + "_z" + str(z+1) + \
                                     "_BW.png")
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
    dict_node_xmerge = {}
    node = 0
    for y in y_range:
        if y == 0:
            ymin = 0
            ymax = y_size + int(options.dy)
        elif y > 0 and y < max(y_range):
            ymin = (y - 1)*y_size - int(options.dy)
            ymax = (y - 1)*y_size + int(options.dy)
        else:
            ymin = (y - 1)*y_size - int(options.dy)
            ymax = input_dim[0]
        for z in z_range:
            node += 1
            dict_node_xmerge[node] = [y, ymin, ymax, z]
    return (node, dict_node_xmerge)

# ----------------------------------------------------------------------

def map_nodes_ymerge(ydim, zdim, options):
    z_range = range(max(0,   int(options.zmin)),\
                    min(zdim,int(options.zmax)))

    y_size = int(round(float(ydim)/float(options.ny)))
    if y_size <= 2*int(options.dy):
        print "\ny overlap ", options.dy, " exceeds half of y_size ", y_size
        sys.exit()
    dict_node_ymerge = {}
    node = 0
    for z in z_range:
        node += 1
        dict_node_ymerge[node] = z
    return (node, dict_node_ymerge)

# ----------------------------------------------------------------------

# Extract fragment data and segment the fragment

def create_data_extraction_and_segmentation_job(outfolderpath, \
                            input_data, input_type, input_dim, \
                            num_nodes, input_label, options):
    print "input_dim[2]=", input_dim[2], " options.zmin=", options.zmin, " options.zmax=", options.zmax
    options.zmin = max(0,            int(options.zmin))
    options.zmax = min(input_dim[2], int(options.zmax))
    fragment_segmentation_script_path = \
        os.path.join(outfolderpath, "Segmentation_script." +\
                     input_label + ".sh")
    executable_extr = os.path.join(ms_home, "Segmentation2D", \
                                   "MS_S2D_ExtractInputData.py")
    executable_segm = os.path.join(ms_home, "Segmentation2D", \
                                   "MS_S2D_Segmentation2D.py")
    if options.verbose:
        print "Number of segm. nodes=", num_nodes

    my_fracBlack = str(options.fracBlack)
    if my_fracBlack == "None":
        my_fracBlack = "0"
    my_fracBlack2 = str(options.fracBlack2)
    if my_fracBlack2 == "None":
        my_fracBlack2 = "0"

    scr = open(fragment_segmentation_script_path, 'wt')
    scr.write("#!/usr/bash\n")
    scr.write("#$ -t 1-" + str(num_nodes) + "\n")
    command_extr_data = \
        "python " + executable_extr + " '" + input_data + "' " + \
        input_type + " -n $SGE_TASK_ID " + \
        " -X " + str(options.nx)   + " -x " + str(options.dx) + \
        " -Y " + str(options.ny)   + " -y " + str(options.dy) + \
        " -L " + str(input_dim[0]) + " -W " + str(input_dim[1]) + \
        " -H " + str(input_dim[2])
    command_segm_data = \
        "python " + executable_segm  + " '" + input_data + "' " +\
        " -o '"+ ms_temp  + "' -n $SGE_TASK_ID  "   +\
        " -X " + str(options.nx)   + " -x " + str(options.dx) + \
        " -Y " + str(options.ny)   + " -y " + str(options.dx) + \
        " -z " + str(options.zmin) + " -Z " + str(options.zmax) +\
        " -l " + options.nlen      + " -w " + options.nwid + \
        " -f " + my_fracBlack      + " -F " + my_fracBlack2           
    if options.unprocessed:
        command_extr_data += " -U "
        command_segm_data += " -U "
    if options.debug:
        command_extr_data += " -D "
        command_segm_data += " -D "
    if options.verbose:
        command_extr_data += " -v "
        command_segm_data += " -v "
    command_extr_data += "\n"
    command_segm_data += "\n"
    scr.write(command_extr_data)
    scr.write(command_segm_data)
    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", \
                  fragment_segmentation_script_path]))
    scr.write("\n")
    scr.close()
    os.chmod(fragment_segmentation_script_path, 0777)
    return fragment_segmentation_script_path

# ----------------------------------------------------------------------

def create_merging_job(mdir, outfolderpath, input_data, input_type,\
                       input_dim, input_label, options):
    options.zmin = max(0,            int(options.zmin))
    options.zmax = min(input_dim[2], int(options.zmax))
    if mdir == 'x':
        num_nodes, dict_node_xmerge = \
            map_nodes_xmerge(input_dim[1], input_dim[0], input_dim[2], options)
        if options.verbose:
            print "Number of x-merge nodes=", num_nodes
    elif mdir == 'y':
        num_nodes, dict_node_ymerge = \
            map_nodes_ymerge(input_dim[0], input_dim[2], options)
        if options.verbose:
            print "Number of y-merge nodes=", num_nodes
    elif mdir == 'z':
        num_nodes = 3
        if options.verbose:
            print "Number of z-merge nodes=", num_nodes
    else:
        print "Unsupported merge direction:", mdir
        sys.exit(2)

    if mdir in ['x', 'y']:
        merge_script_path = \
            os.path.join(outfolderpath, "Merge_in_" + mdir + "_dir_script." +\
                         input_label + ".sh")
    else:
        merge_script_path = \
            os.path.join(outfolderpath, "Epilog_script." +\
                         input_label + ".sh")
    executable_path = os.path.join(ms_home,"Segmentation2D", \
                                   "MS_S2D_MergeImageFragments.py")
    scr = open(merge_script_path, 'wt')
    scr.write("#!/usr/bash\n")
    scr.write("#$ -t 1-" + str(num_nodes) + "\n")
    command_merge = executable_path + " " + input_label + " " + mdir \
                      + " " + str(input_dim[0]) + " " + str(input_dim[1]) \
                      + " " + str(input_dim[2]) + " $SGE_TASK_ID " \
                      + " -X " + str(options.nx)   \
                      + " -Y " + str(options.ny)   \
                      + " -z " + str(options.zmin) \
                      + " -Z " + str(options.zmax) 
    if options.verbose:
        command_merge += " -v "
    if options.debug:
        command_merge += " -D"
    command_merge += "\n"
    if options.verbose:
        print "command_merge=", command_merge
    scr.write(command_merge) 
    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", \
                  merge_script_path]))
    scr.close()
    os.chmod(merge_script_path, 0777)
    return merge_script_path

# -----------------------------------------------------------------------------

def submit_segmentation_job(segm_script_path, base_command, options):
    prog_name = "ms2.segm"
    jobid1 = 0
    if options.submission_command == "qsub":
        command = base_command + " -V -N " + prog_name +\
                                 " -o /dev/null -e /dev/null "
    else:
        command = base_command + " -V -N " + prog_name
    if options.num_slots > 1:
        command += " -pe batch " + str(options.num_slots) 
    if re.search("qsub", options.submission_command):
        command += " " + segm_script_path   
        res = commands.getstatusoutput(command)
        jobid1 = (res[1].split()[2]).split(".")[0]
    elif not options.submission_command == "none":
        os.system(command)
    if options.verbose and not options.submission_command == "none":
        print "Submit segmentation job command=", command, "\nres=", res

    # Submit epilog script
    return jobid1            

#------------------------------------------------------------------------------

def submit_xmerge_job(xmerge_script_path, base_command, jobid1, options):
    prog_name = "ms2.xmrg" 
    jobid2 = 0
    if options.submission_command == "qsub":
        command = base_command + " -V -N " + prog_name +\
                                 " -o /dev/null -e /dev/null "
    else:
        command = base_command + " -V -N " + prog_name
    if re.search("qsub", options.submission_command):
        if int(jobid1) > 0:
            command += " -hold_jid " + str(jobid1)
        command += " " + xmerge_script_path
        res2 = commands.getstatusoutput(command)
        jobid2 = (res2[1].split()[2]).split(".")[0]
    elif not options.submission_command == "none":
        command += " " + xmerge_script_path
        os.system(command)
    if options.verbose and not options.submission_command == "none":
        print "\nSubmit x-merge job command=", command, " res2=", res2
        print "jobid2=", jobid2
    return jobid2

# -----------------------------------------------------------------------------

def submit_ymerge_job(ymerge_script_path, base_command, jobid2, options):
    prog_name = "ms2.ymrg"  
    jobid3 = 0
    if options.submission_command == "qsub":
        command = base_command + " -V -N " + prog_name +\
                                 " -o /dev/null -e /dev/null "
    else:
        command = base_command + " -V -N " + prog_name
    if re.search("qsub", options.submission_command):
        if int(jobid2) > 0:
            command += " -hold_jid " + str(jobid2)
        command += " " + ymerge_script_path
        res3 = commands.getstatusoutput(command)
        jobid3 = (res3[1].split()[2]).split(".")[0]
            
    elif not options.submission_command == "none":
        command += " " + ymerge_script_path
        os.system(command)
    if options.verbose and not options.submission_command == "none":
        print "\nSubmit y-merge job command=", command, " res3=", res3
        print "jobid3=", jobid3
    return jobid3

# -----------------------------------------------------------------------------

def submit_epilog_job(epilog_script_path, base_command, jobid3, options):
    prog_name = "ms2.epilog"
    command = base_command + " -V -N " + prog_name
    if options.submission_command == "qsub":
        command += " -o /dev/null -e /dev/null "
    if int(jobid3) > 0:
        command += " -hold_jid " + str(jobid3)
    command += " " +  epilog_script_path
    if not options.submission_command == "none":
        os.system(command)
    if options.verbose and not options.submission_command == "none":
        print "\nSubmit epilog job command=", command
    return 

# -----------------------------------------------------------------------------

def submit_array_jobs(segm_script_path,  \
                      xmerge_script_path,\
                      ymerge_script_path,\
                      epilog_script_path, options):
    qsub = os.path.join(os.environ['SGE_ROOT'], "bin", "lx-amd64", "qsub")       
    if re.search("qsub", options.submission_command):
        base_command = qsub
    else:
        base_command = options.submission_command
    if len(options.project_code) > 0:
        base_command += "  -A " + options.project_code

    jid1 = 0
    if int(options.processing_start) == 1 and int(options.processing_end) >= 1:
        jid1 = submit_segmentation_job(segm_script_path, base_command, options)

    jid2 = 0
    if int(options.processing_start) <= 2 and int(options.processing_end) >= 2:
        jid2 = submit_xmerge_job(xmerge_script_path,base_command,jid1, options)

    jid3 = 0
    if int(options.processing_start) <= 3 and int(options.processing_end) >= 3:
        jid3 = submit_ymerge_job(ymerge_script_path,base_command,jid2, options)

    if int(options.processing_start) <= 4 and int(options.processing_end) >= 4:
        submit_epilog_job(epilog_script_path, base_command, jid3, options)

    return

# -----------------------------------------------------------------------------

# Create and submit all cluster jobs
def process_input_data_high_level(input_data, input_type, input_dim, \
                                  num_nodes, input_label, options):
    outfolderpath = ms_temp
   
    if int(options.processing_start) == 1:
        segmentation_script_path = \
            create_data_extraction_and_segmentation_job(outfolderpath, \
                                input_data, input_type, input_dim, \
                                num_nodes, input_label, options)
    else:
        segmentation_script_path = ""

    if int(options.processing_start) <= 2:
        fragment_xmerge_script_path = \
            create_merging_job('x', outfolderpath, input_data,\
                               input_type, input_dim, input_label, options)
    else:
        fragment_xmerge_script_path = ""

    if int(options.processing_start) <= 3:
        fragment_ymerge_script_path = \
            create_merging_job('y', outfolderpath, input_data,\
                               input_type, input_dim, input_label, options)
    else:
        fragment_ymerge_script_path = ""

    if int(options.processing_start) <= 4:
        epilog_script_path = \
            create_merging_job('z', outfolderpath, input_data, \
                               input_type, input_dim, input_label, options)
    else:
        epilog_script_path = ""

    if options.verbose:
        print "segmentation_script_path=", segmentation_script_path
        print "fragment_xmerge_script_path=", fragment_xmerge_script_path
        print "fragment_ymerge_script_path=", fragment_ymerge_script_path
        print "epilog_script_path=", epilog_script_path
    submit_array_jobs(segmentation_script_path,\
                      fragment_xmerge_script_path,\
                      fragment_ymerge_script_path,\
                      epilog_script_path, options)

# -----------------------------------------------------------------------------

# Process individual fragment segmentation job                 
def process_input_data_low_level(dict_node_xyz, input_label,options):
    outfolderpath = ms_temp     
    executable_path = os.path.join(ms_home, "Bin", options.executable)
    node = int(options.node)
    y    = dict_node_xyz[node][0]
    ymin = dict_node_xyz[node][1]
    ymax = dict_node_xyz[node][2]
    x    = dict_node_xyz[node][3]
    xmin = dict_node_xyz[node][4]
    xmax = dict_node_xyz[node][5]
    z    = dict_node_xyz[node][6]
    input_file = os.path.join(ms_temp, input_label + "_y" + str(y+1) +\
                              "_x" + str(x+1) + "_z" + str(z+1) + ".png")
    output_file = os.path.join(ms_temp, input_label + "_y" + str(y+1) +\
                              "_x" + str(x+1) + "_z" + str(z+1) + "_BW.png")
    my_fracBlack = str(options.fracBlack)
    if my_fracBlack == "None":
        my_fracBlack = "0"
    my_fracBlack2 = str(options.fracBlack2)
    if my_fracBlack2 == "None":
        my_fracBlack2 = "0"
    command_segm = executable_path + " " + input_file + " " + \
                   my_fracBlack    + " " + my_fracBlack2 + " " +\
                   " verbose "     + str(int(options.verbose)) + " " + \
                   " ny      "     + options.nlen              + " " + \
                   " nx      "     + options.nwid              + " " + \
                   " outBW   "     + output_file
    command_rm   = "rm -f " + input_file
    if int(options.msize) < sys.maxint:             
        command_segm +=  ' maxSize '    + options.msize   
    if options.verbose:
        print "Input data file=", input_file
        data_shape = misc.imread(input_file).shape
        print "Input data shape=", data_shape
        print "low level command_segm=", command_segm
    os.system(command_segm)
    if not options.debug:
        os.system(command_rm)

# -----------------------------------------------------------------------------

if __name__ == "__main__":
   
    usage = "Usage: \n\
    %prog input_data [options (-h to list)] \nwhere\n\
    input_data=\n\
    1) path to image stack file in HDF5 format, or \n\
    2) path to folder containing image files, or \n\
    3) coordinates of DVID subvolume: [xmin:xmax,ymin:ymax,zmin:zmax]"   

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = segmentation_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if options.compile or options.compile_all:
        compile_code(options)
    if len(args) == 1:
        input_data = args[0]
        input_type = get_input_type(input_data, options)
        input_label = "ms2_DVID"
        if input_type in ["file", "directory"]:
            input_label = "ms2_" + input_data
        if input_type in ["file", "directory", "DVID"]:
            if options.verbose and int(options.node) == 0:
                print "Input type=", input_type
        else:
            print "\nIncorrectly specified input data", input_data, "\n"
            parser.print_usage()
            sys.exit(2)         
        input_dim = get_data_dimensions(input_data, input_type, options)
        if options.verbose:
            print "Input data dimensions: ", input_dim
        num_nodes, dict_node_xyz = \
            map_node_to_xyz(input_dim, input_label, options)
        if options.verbose:
            print "num_nodes=", num_nodes
            print "\nProcessing input data ..."
        if int(options.node) == 0 and num_nodes > 0:
            process_input_data_high_level(input_data, input_type, input_dim, \
                num_nodes, input_label, options)
        elif num_nodes > 0:
            process_input_data_low_level(dict_node_xyz, input_label, options)
    elif options.compile or options.compile_all:
        sys.exit(2)
    else:
        parser.print_usage()
        sys.exit(2)
