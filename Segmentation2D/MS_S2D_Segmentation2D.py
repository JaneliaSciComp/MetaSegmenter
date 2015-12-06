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

import MS_LIB_Dict 
import MS_LIB_Options
import MS_LIB_IO

ms_home = os.environ['MS_HOME'] # where source code is located
ms_data = os.environ['MS_DATA'] # dir for initial input data and final output 
ms_temp = os.environ['MS_TEMP'] # dir for intermediate data files

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

# Extract fragment data and segment the fragment

def create_data_extraction_and_segmentation_job(outfolderpath, \
                            input_data, input_type, input_dim, \
                            num_nodes, input_label, options):
    print "input_dim[2]=", input_dim[2], " options.zmin=", options.zmin, " options.zmax=", options.zmax

    zmin = max(0,                   int(options.zmin))
    zmax = min(zmin + input_dim[2], int(options.zmax))

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
        " -z " + str(zmin)         + " -Z " + str(zmax) +\
        " -L " + str(input_dim[0]) + " -W " + str(input_dim[1]) + \
        " -H " + str(input_dim[2])
    command_segm_data = \
        "python " + executable_segm  + " '" + input_data + "' " +\
        " -o '"+ ms_temp  + "' -n $SGE_TASK_ID  "   +\
        " -X " + str(options.nx)   + " -x " + str(options.dx) + \
        " -Y " + str(options.ny)   + " -y " + str(options.dx) + \
        " -z " + str(zmin)         + " -Z " + str(zmax) +\
        " -l " + str(options.nlen) + " -w " + str(options.nwid) + \
        " -f " + my_fracBlack      + " -F " + my_fracBlack2 +\
        " -r " + str(options.resize_scale) +\
        " -b " + str(options.no_dark)          
    if len(options.memb_prob) > 0:
        memb_prob_dir_path = os.path.join(ms_data, options.memb_prob)
        command_segm_data += " -M '" + memb_prob_dir_path + "'"
        command_extr_data += " -M '" + memb_prob_dir_path + "'"
    if len(options.mito_prob) > 0:
        mito_prob_dir_path = os.path.join(ms_data, options.mito_prob)
        command_segm_data += " -m '" + mito_prob_dir_path + "'"
        command_extr_data += " -m '" + mito_prob_dir_path + "'"
    if options.unprocessed:
        command_extr_data += " -U "
        command_segm_data += " -U "
    if options.debug:
        command_extr_data += " -D "
        command_segm_data += " -D "
    if options.verbose:
        command_extr_data += " -v "
        command_segm_data += " -v "
    if options.output_thresholds:
        command_segm_data += " -T "
    if options.use_memb_probs:
        command_segm_data += " -u "
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

    zmin = max(0,            int(options.zmin))
    zmax = min(input_dim[2], int(options.zmax))
    if mdir == 'x':
        num_nodes, dict_node_xmerge = \
            MS_LIB_Dict.map_nodes_xmerge(input_dim[1], input_dim[0], input_dim[2], options)
        if options.verbose:
            print "Number of x-merge nodes=", num_nodes
    elif mdir == 'y':
        num_nodes, dict_node_ymerge = \
            MS_LIB_Dict.map_nodes_ymerge(input_dim[0], input_dim[2], options)
        if options.verbose:
            print "Number of y-merge nodes=", num_nodes
    elif mdir == 'z':
        num_nodes = 2
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
                      + " -z " + str(zmin) \
                      + " -Z " + str(zmax) \
                      + " -i " + str(options.uint) 
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
        res = ''
    if options.verbose and not options.submission_command == "none":
        print "Submit segmentation job command=", command, "\nres=", res, "\n"
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
        print "\nSubmit x-merge job command=", command, " res2=", res2, "\n"
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
        if options.num_slots > 1:
            command += " -pe batch " + str(options.num_slots)
        command += " " + ymerge_script_path
        if options.verbose:
            print "Submitting ymerge command=", command
        res3 = commands.getstatusoutput(command)
        jobid3 = (res3[1].split()[2]).split(".")[0]
    elif not options.submission_command == "none":
        command += " " + ymerge_script_path
        os.system(command)
        res3 = ''
    if options.verbose and not options.submission_command == "none":
        print "\nSubmit y-merge job command=", command, " res3=", res3, "\n"
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
    if options.nx > 1 and int(options.processing_start) <= 2 and int(options.processing_end) >= 2:
        jid2 = submit_xmerge_job(xmerge_script_path,base_command,jid1, options)

    if options.nx == 1:
        jid2 = jid1

    jid3 = 0
    # if options.ny == 1, still submit ymerge_job, but for only producing Seg/RGB files
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

    if options.nx > 1 and int(options.processing_start) <= 2:
        fragment_xmerge_script_path = \
            create_merging_job('x', outfolderpath, input_data,\
                               input_type, input_dim, input_label, options)
    else:
        fragment_xmerge_script_path = ""

    # if options.ny == 1, still create ymerge_job, but for only producing Seg/RGB files
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

def get_iput_path_list(input_data, options):
    # Create a list of files
    files = []
    i = 0
    dir_path = os.path.join(ms_data,input_data)
    files_list = sorted(os.listdir(dir_path))
    print "\ninput_data=", input_data
    print "files_list=", files_list
    files = []
    for i in range(0, len(files_list)):
        file = files_list[i]
        file_path = os.path.join(dir_path, file)
        if i in range(int(options.zmin), int(options.zmax)):
            files.append(os.path.join(input_data, file_path))
        i = i+1
    return files

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
    if options.verbose:
        print "options.nx=", options.nx, " options.ny=", options.ny, \
              "ymin=", ymin, " ymax=", ymax, \
              "xmin=", xmin, " xmax=", xmax, \
              "z=", z
    if   int(options.nx) > 1 and int(options.ny) > 1:
        input_path = os.path.join(ms_temp, input_label + "_y" + str(y+1) +\
                                  "_x" + str(x+1) + "_z" + str(z+1) + ".png")
        output_path = os.path.join(ms_temp, input_label + "_y" + str(y+1) +\
                                  "_x" + str(x+1) + "_z" + str(z+1) + "_BW.png")
    elif int(options.nx) == 1 and int(options.ny) > 1:
        input_path = os.path.join(ms_temp, input_label + "_y" + str(y+1) +\
                                                    "_z" + str(z+1) + ".png")
        output_path = os.path.join(ms_temp, input_label + "_y" + str(y+1) +\
                                                    "_z" + str(z+1) + "_BW.png")
    elif int(options.nx) > 1 and int(options.ny) == 1:
        input_path = os.path.join(ms_temp, input_label +                  \
                                  "_x" + str(x+1) + "_z" + str(z+1) + ".png")
        output_path = os.path.join(ms_temp, input_label +                  \
                                  "_x" + str(x+1) + "_z" + str(z+1) + "_BW.png")
    else:
        input_path = os.path.join(ms_temp, input_label +                  \
                                                    "_z" + str(z+1) + ".png")
        output_path = os.path.join(ms_temp, input_label +                  \
                                                    "_z" + str(z+1) + "_BW.png")
    if options.output_thresholds:
        if  input_type == "directory":
            file_paths = get_iput_path_list(input_data, options)
            print "z=", z, " file_paths[z]=",  file_paths[z]
            input_path = file_paths[z]
        input_file = os.path.basename(input_path)
        input_file_suffix = "." + input_file.split(".")[len(input_file.split("."))-1]
        output_file = input_file.replace(input_file_suffix, ".thr")
        output_path = os.path.join(ms_data, output_file)

    if options.verbose:
        print "input_path=", input_path, "\noutput_path=", output_path

    if len(options.memb_prob) > 0:
        input_memb_prob_path = os.path.join(ms_temp, input_label +          \
                                            "_z" + str(z+1) + "_memb_prob.h5")
        if options.verbose:
            print "memb_prob_path =", input_memb_prob_path   

    if len(options.mito_prob) > 0:
        input_mito_prob_path = os.path.join(ms_temp, input_label +          \
                                            "_z" + str(z+1) + "_mito_prob.h5")
        if options.verbose:
            print "mito_prob_path =", input_mito_prob_path
      
    my_fracBlack = str(options.fracBlack)
    if my_fracBlack == "None":
        my_fracBlack = "0"
    my_fracBlack2 = str(options.fracBlack2)
    if my_fracBlack2 == "None":
        my_fracBlack2 = "0"
    command_segm = executable_path + " " + input_path + " " + \
                   my_fracBlack    + " " + my_fracBlack2 + " " +\
                   " verbose "     + str(int(options.verbose)) + " " + \
                   " ny      "     + options.nlen              + " " + \
                   " nx      "     + options.nwid              + " " + \
                   " outBW   "     + output_path               + " " + \
                   " resize  "     + str(options.resize_scale)
    if int(options.no_dark) > 0:
        command_segm += " noDark  1 "
    if len(options.memb_prob) > 0:
        command_segm += " membPr " + input_memb_prob_path
    if len(options.mito_prob) > 0:
        command_segm += " mitoPr " + input_mito_prob_path
    if options.output_thresholds:
        command_segm += " outThr " + output_path
    if options.use_memb_probs:
        command_segm += " useMembPr 1 "

    command_rm   = "rm -f " + input_path
    if int(options.msize) < sys.maxint:             
        command_segm +=  ' maxSize '    + options.msize   
    if options.verbose:
        print "Input data file=", input_path
        data_shape = misc.imread(input_path).shape
        print "Input data shape=", data_shape
        print "low level command_segm=", command_segm
        print "Writing segmented image to file", output_path 
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
    parser = MS_LIB_Options.Segmentation2D_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if options.compile or options.compile_all:
        compile_code(options)
    if len(args) == 1:
        input_data = args[0]
        input_type = MS_LIB_IO.get_input_type(input_data, options)
        input_label = "ms2_DVID"
        if input_type in ["file", "directory"]:
            input_label = "ms2_" + input_data.split('.')[0]
        if input_type in ["file", "directory", "DVID"]:
            if options.verbose and int(options.node) == 0:
                print "Input type=", input_type, " input_data=", input_data
        else:
            print "\nIncorrectly specified input data", input_data, "\n"
            parser.print_usage()
            sys.exit(2)         
        input_dim = MS_LIB_IO.get_data_dimensions(input_data, input_type, options)
        if options.verbose:
            print "Input data dimensions: ", input_dim
        if options.output_thresholds:
            options.ny = 1
            options.nx = 1
            options.dy = 0
            options.dx = 0
        num_nodes, dict_node_xyz = \
            MS_LIB_Dict.map_node_to_xyz(input_dim, input_label, "_BW.png", options)
        if options.verbose:
            print "\nint(options.node)=", int(options.node), " num_nodes=", num_nodes
            print "options.nlen=", options.nlen
#           print "dict_node_xyz =", dict_node_xyz
        if int(options.node) == 0 and num_nodes > 0:
            print "\nProcessing input data high level..."
            process_input_data_high_level(input_data, input_type, input_dim, \
                num_nodes, input_label, options)
        elif num_nodes > 0:
            print "\nProcessing input data low  level..."
            process_input_data_low_level(dict_node_xyz, input_label, options)
        else:
            sys.exit("\nCannot process, since num_nodes == 0")
    elif options.compile or options.compile_all:
        sys.exit(2)
    else:
        parser.print_usage()
        sys.exit(2)
