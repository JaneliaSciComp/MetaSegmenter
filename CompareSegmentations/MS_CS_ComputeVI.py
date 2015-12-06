#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

# Purpose: compare segmentation stacks using the VI criterion
# as defined in:
# Marina Meila, "Comparing clusterings.an information based distance"
# Journal of Multivariate Analysis 98 (2007) 873 - 895
#

import os, sys, re 
import h5py
import numpy
import json 
import ntpath
import optparse
import commands

import MS_LIB_Dict
import MS_LIB_Options

ms_home = os.environ['MS_HOME'] # where source code is located
ms_data = os.environ['MS_DATA'] # dir for initial input data and final output
ms_temp = os.environ['MS_TEMP'] # dir for intermediate data files

# -------------------------------------------------------------------------------

def read_input(stack_file_name):
    if re.search(".h5", stack_file_name):
        f = h5py.File(stack_file_name, 'r')
        key = f.keys()[0]
        data = numpy.uint64(numpy.transpose(f[key]))
    elif re.search(".tif", stack_file_name):
        data = numpy.uint64(numpy.transpose(tiff.imread(stack_file_name)))
    else:
        sys.exit("Unrecognized stack type")
    return data

# -------------------------------------------------------------------------------

def create_compute_counts_job(num_layers, stack1, stack2, options):
    num_nodes, dict_node_z = \
            MS_LIB_Dict.map_node_to_z(num_layers, options)

    label1 = ntpath.basename(stack1).split('.')[0]
    label2 = ntpath.basename(stack1).split('.')[0]

    compute_counts_script_path = \
        os.path.join(ms_temp, "VI_ComputeCounts_script_" +\
                     label1 + "_" + label2 + ".sh")
    executable_counts = os.path.join(ms_home, "CompareSegmentations",\
                                     "MS_CS_ComputeVI.py")
    if options.verbose:
        print "Number of compute counts nodes=", num_nodes
    scr = open(compute_counts_script_path, 'wt')
    scr.write("#!/usr/bash\n")
    scr.write("#$ -t 1-" + str(num_nodes) + "\n")
    command_counts =  "python " + executable_counts + \
        "    '"+ stack1 + "' " + " '"+ stack2 + "' " + \
        " -n $SGE_TASK_ID " + \
        " -z " + str(options.zmin)         + " -Z " + str(options.zmax) 
    if options.verbose:
        command_counts += " -v "
    if options.debug:
        command_counts += " -D "
    command_counts       += "\n"
    scr.write(command_counts)          
    if options.verbose:
        print "In create_compute_counts_job: command=", command_counts   
    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", \
                  compute_counts_script_path]))
    scr.write("\n")
    scr.close()
    os.chmod(compute_counts_script_path, 0777)
    return compute_counts_script_path

# -------------------------------------------------------------------------------

def create_epilog_job(ms_data, label1, label2, options):
    epilog_script_path = \
        os.path.join(ms_temp, "VI_Epilog_script_" +\
                     label1 + "_" + label2 + ".sh")
    scr = open(epilog_script_path, 'wt')

    if len(options.output_file) > 0:
        output_path = options.output_file
    else:
        output_path = os.path.join(ms_data, "VI" + str(options.dim) + \
                                   "D_" + label1 + "_" + label2 + \
                                   "_z" + str(options.zmin) + "-" +str(options.zmax) + ".txt")

    scr.write("#!/usr/bash\n")
    command_epilog = ""
    if int(options.dim) == 3:
        executable_epilog = os.path.join(ms_home, "CompareSegmentations",\
                                         "MS_CS_CountsToVI3D.py")
    elif int(options.dim) == 2:
        executable_epilog = os.path.join(ms_home, "CompareSegmentations",\
                                         "MS_CS_CountsToVI2D.py")
    else:
        sys.exit("Incorrectly specified VI dim")
    command_epilog = executable_epilog + " " + label1  + " " + label2 + \
                      " " + str(options.zmin) + " " + str(options.zmax) + " " + output_path
    if options.verbose:
        print "command_epilog=", command_epilog
    command_epilog += "\n"
    scr.write(command_epilog)
    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", \
                  epilog_script_path]))
    scr.close()
    os.chmod(epilog_script_path, 0777)
    return epilog_script_path

# -------------------------------------------------------------------------------

def get_counts(z, data1, data2, options):
    N1  = {}
    N2  = {}
    N12 = {}

    # Extract counts
    im1 = data1[:,:,z]
    im2 = data2[:,:,z]
    print "Creating list of labels1 ..."
    labels1 = numpy.unique(im1)            
    print " ... len(labels1)=", len(labels1)
    print "Creating list of labels2 ..."
    labels2 = numpy.unique(im2)      
    print " ... len(labels2)=", len(labels2)
    print "Computing counts N1 ..."
    for k1 in labels1:
        count = (im1 == k1).sum()
        if not k1 in N1.keys():
            N1[str(k1)]  = count
        else:
            N1[str(k1)] += count
    print "Computing counts N2 ..."
    for k2 in labels2:
        count = (im2 == k2).sum()
        if not k2 in N2.keys():
            N2[str(k2)]  = count
        else:
            N2[str(k2)] += count
    print "Computing counts N12 ..."
    for k1 in labels1:
        for k2 in labels2:
            k12 = str(k1) + "_" + str(k2)
            count = ((im1 == k1) & (im2 == k2)).sum()
            if not k12 in N12.keys():
                N12[k12]  = count
            else:
                N12[k12] += count
    if options.verbose:
        print "len(labels1)=", len(labels1), " len(keys1)=", len(N1.keys()), \
             " len(labels2)=", len(labels2), " len(keys2)=", len(N2.keys())
    if options.verbose:
        print "\nN1.values=", N1.values()
        print "\nN2.values=", N2.values()
        print "\nN12.values=", N12.values()

    return (N1, N2, N12)

# -------------------------------------------------------------------------------

def output_counts(z, label1, label2, N1, N2, N12, options):
    output_path1  = os.path.join(ms_temp, "VI_" + label1 + "_" + label2 + "_N1_counts_z"  + str(z) + ".json")
    with open(output_path1, 'w') as fp:
        json.dump(N1, fp)

    output_path2  = os.path.join(ms_temp, "VI_" + label1 + "_" + label2 + "_N2_counts_z"  + str(z) + ".json")
    with open(output_path2, 'w') as fp:
        json.dump(N2, fp)

    output_path12 = os.path.join(ms_temp, "VI_" + label1 + "_" + label2 + "_N12_counts_z" + str(z) + ".json")
    with open(output_path12, 'w') as fp:
        json.dump(N12, fp)

# -------------------------------------------------------------------------------

def submit_job(script_path, base_command, jobid_in, options):
    prog_name = ""
    print "script_path=", script_path
    if re.search("Counts",  script_path):
        prog_name = "vi.counts"
    elif re.search("Epilog", script_path):
        prog_name = "vi.epilog"
    jobid_out = 0
    if options.submission_command == "qsub":
        command = base_command + " -V -N " + prog_name +\
                                 " -o /dev/null -e /dev/null "
    elif re.search("qsub", options.submission_command):
        command = base_command + " -V -N " + prog_name
    else:
        command = base_command
    res = ""
    if re.search("qsub", options.submission_command):
        if jobid_in > 0:
            command += " -hold_jid " + str(jobid_in)
        command += " " + script_path
        if not re.search("Epilog", script_path):
            res = commands.getstatusoutput(command)
            jobid_out = (res[1].split()[2]).split(".")[0]
        else:
            os.system(command)
            jobid_out = 0
    elif not options.submission_command == "none":
        command += " " + script_path
        os.system(command)
    if  options.verbose:
        print "\nSubmit job command=", command, " res=", res
        print "jobid_out=", jobid_out
    return jobid_out

# -------------------------------------------------------------------------------

def submit_all_jobs(compute_counts_script_path, epilog_script_path, options):
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
        jid1 = submit_job(compute_counts_script_path, base_command, 0, options)

    if int(options.processing_step_beg) <= 2 and \
       int(options.processing_step_end) >= 2:
        submit_job(epilog_script_path, base_command, jid1, options)

    return

# -------------------------------------------------------------------------------

def process_inputs_high_level(num_layers, stack1, stack2, options):

    label1 = ntpath.basename(stack1).split('.')[0]
    label2 = ntpath.basename(stack2).split('.')[0]

    print "options.processing_step_beg=", options.processing_step_beg
    if int(options.processing_step_beg) == 1:
        compute_counts_script_path = \
            create_compute_counts_job(num_layers, stack1, stack2, options)
    else:
        compute_counts_script_path = ""

    if int(options.processing_step_end) == 2:
        epilog_script_path = \
            create_epilog_job(ms_data, label1, label2, options) 
    else:
        epilog_script_path = ""

    if options.verbose:
        print "compute_counts_script_path=", compute_counts_script_path    
        print "epilog_script_path=",    epilog_script_path
    submit_all_jobs(compute_counts_script_path, epilog_script_path, options)

# -----------------------------------------------------------------------------

def process_inputs_low_level(num_layers, stack1, stack2, options):
    num_nodes, dict_node_z = \
            MS_LIB_Dict.map_node_to_z(num_layers, options)

    outfolderpath = ms_temp
    node = int(options.node)
    z    = dict_node_z[node]
    if options.verbose:
        print "z=", z

    # Read input data
    f1 = h5py.File(stack1, 'r')
    key = f1.keys()[0]
    data1 = numpy.transpose(f1[key])
#   data1 = numpy.squeeze(data[:, :, z])
    
    f2 = h5py.File(stack2, 'r')
    key = f2.keys()[0]
    data2 = numpy.transpose(f2[key])
#   data2 = numpy.squeeze(data[:, :, z])

    # Compute counts
    N1, N2, N12 = get_counts(z, data1, data2, options)

    print "data1.shape=", data1.shape, " data2.shape=", data2.shape
    # Output results
    label1 = ntpath.basename(stack1).split('.')[0]
    label2 = ntpath.basename(stack2).split('.')[0]
    output_counts(z, label1, label2, N1, N2, N12, options)

# -----------------------------------------------------------------------------

if __name__ == "__main__":

    usage = "Usage: \n\
    %prog seg_stack1 seg_stack2 [ node ] [ options (-h to list) ] "

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = MS_LIB_Options.ComputeVI_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if len(args) == 2:             
        stack_file_name1 = args[0]
        stack_file_name2 = args[1]
    else:
        parser.print_usage()
        sys.exit(2)

    if options.verbose:
        print "stack_file_name1=", stack_file_name1
        print "stack_file_name2=", stack_file_name2

    f1 = h5py.File(stack_file_name1, 'r')
    key = f1.keys()[0]
    num_layers1 = int(f1[key].shape[0])

    f2 = h5py.File(stack_file_name2, 'r')
    key = f2.keys()[0]
    num_layers2 = int(f2[key].shape[0])

    if not num_layers1 == num_layers2:
        print "num_layers1=", num_layers1, " num_layers2=", num_layers2
        sys.exit('Number of layers in the input files must be the same')
    else:
        options.zmin = max(int(options.zmin), 0)
        options.zmax = min(int(options.zmax), num_layers1)
    
    stack1_path = os.path.join(ms_data, stack_file_name1)
    stack2_path = os.path.join(ms_data, stack_file_name2)

    if int(options.node) == 0:
        process_inputs_high_level(num_layers1, stack1_path, stack2_path, options)
    else:
        process_inputs_low_level( num_layers1, stack1_path, stack2_path, options)


