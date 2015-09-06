#! /usr/local/python-2.7.6/bin/python     
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

# Purpose: submit Matlab's 2D segmentation job(s) to cluster,
#          then post-process results

# ------------------------- imports -------------------------

import os 
import sys, re, optparse
import numpy
import h5py
from PIL import Image
import mimetypes
import ntpath

import MS_LIB_Dict 
import MS_LIB_Options
import MS_LIB_IO

ms_home = os.environ['MS_HOME'] # where source code is located
ms_data = os.environ['MS_DATA'] # dir for initial input data and final output 
ms_temp = os.environ['MS_TEMP'] # dir for intermediate data files
rh_home = os.environ['RHOANA_HOME'] # where source code is located

# -----------------------------------------------------------------------

def produce_segmentation(input_data, probs_file_list, options):
    segmentation_script = os.path.join(rh_home, "rhoana", \
        "Segment", "segment.py")
    for i in range(0, len(probs_file_list)):
        probs_file = ntpath.basename(probs_file_list[i])
        output_file = probs_file.split(".")[0] + ".h5"
        if re.search("mitochondria", options.classifier):
            output_dir = os.path.join(ms_data, "segmentation_2D_bw_mitochondria_" \
                                      + options.output_tag)
        else:
            output_dir = os.path.join(ms_data, "segmentation_2D_bw_membranes_" \
                                      + options.output_tag)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        command = "python " + segmentation_script + " " + probs_file_list[i] +\
                  " " + os.path.join(output_dir, output_file)
        if options.verbose:
            print "\nSegmentation command=", command, "\n"
        os.system(command)
    return 

# -----------------------------------------------------------------------------

def compute_segmented_object_probabilities(input_data, file_list, options):
    classify_image_executable = options.executable
    classifier_file = options.classifier
    output_file_list = []
    for i in range(0, len(file_list)):
        input_file_path  = file_list[i]
        print "input_file_path=", input_file_path
        input_file       = ntpath.basename(input_file_path)     
        output_file      = input_file.split(".")[0] + ".h5"
        input_tag = input_data
        if input_data[:4] == "raw_":
            input_tag = input_data[4:]
        output_dir    = os.path.join(ms_data, options.seg_type + "_prob_" + input_tag)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        output_file_path = os.path.join(output_dir, output_file)
        output_file_list.append(output_file_path)
        command = classify_image_executable + " " + input_file_path + " " + \
                  classifier_file + " " + output_file_path     
        if options.verbose:
            print "...command=", command
        os.system(command) 
    return output_file_list

# -----------------------------------------------------------------------------

def extract_file_list_from_directory(input_dir):
    file_list = []
    files = os.listdir(input_dir)
    for file in files:
        file_path = os.path.join(ms_data, input_dir, file)
        if mimetypes.guess_type(file_path)[0] in \
        ['image/png', 'image/jpeg', 'image/tiff']:
            file_list.append(file_path)
    if len(file_list) > 0:
        file_list = sorted(file_list)
    return file_list

# -----------------------------------------------------------------------------

def extract_file_list_from_fileoffiles(input_data):
    file_list = []
    infile = open(input_data, mode='r')
    while True:
        line = infile.readline().strip()
        if os.path.isfile(line) and mimetypes.guess_type(input_data)[0] in \
        ['image/png', 'image/jpeg', 'image/tiff']:
            file_list.append(line)
    if len(file_list) > 0:
        file_list = sorted(file_list)
    return file_list

# -----------------------------------------------------------------------------

# Process individual fragment segmentation job                 
def process_input_data(input_data, input_type, input_label, options):
    input_file_list = []
    if input_type == "file" and mimetypes.guess_type(input_data)[0] in \
        ['image/png', 'image/jpeg', 'image/tiff']:
        input_file_list = [ os.path.join(ms_data, input_data) ]
    elif input_type == "fileoffiles":
        input_file_list = extract_file_list_from_fileoffiles(input_data)
    elif input_type == "directory":
        input_file_list = extract_file_list_from_directory(input_data)
    else:
        print "\nUnrecognizable input type"
        sys.exit()
    
    # Filter input list
    filtered_input_file_list = []
    for i in range(0, len(input_file_list)):
        if i >= int(options.zmin) and i <= int(options.zmax):
            filtered_input_file_list.append(input_file_list[i]) 
    if options.verbose: 
        print "input_file_list=", filtered_input_file_list 

    if len(input_file_list) > 0:
        if int(options.processing_start) == 1:
            probs_file_list = compute_segmented_object_probabilities(input_data, \
                                  filtered_input_file_list, options)

        if int(options.processing_end) == 2:
            produce_segmentation(input_data, probs_file_list, options)
    return

# -----------------------------------------------------------------------------

if __name__ == "__main__":
   
    usage = "Usage: \n\
    %prog input_data [options (-h to list)] \nwhere\n\
    input_data=\n\
    1) path to image stack file in HDF5 format, or \n\
    2) path to folder containing image files, or \n\
    3) coordinates of DVID subvolume: [xmin:xmax,ymin:ymax,zmin:zmax]"   

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = MS_LIB_Options.RhoanaSegmentation2D_command_line_parser(parser)
    (options, args) = parser.parse_args()

    # Input arguments:
    # 1) input data (=file in HDF5 format, or file of files, or directory, or DVID entry)
    # 2) 
    # Input options:
    # 1) Classifier file (if different from default
    # 2) Number of segmentations

    if len(args) == 1:
        input_data = args[0]
        input_type = MS_LIB_IO.get_input_type(input_data, options)
        input_label = "rh2_DVID"
        if input_type in ["file", "fileoffiles", "directory"]:
            input_label = "rh2_" + input_data.split('.')[0]
        if input_type in ["file", "fileoffiles", "directory", "DVID"]:
            if options.verbose:
                print "Input type=", input_type
        else:
            print "\nIncorrectly specified input data", input_data, "\n"
            parser.print_usage()
            sys.exit(2)         
        if options.verbose:
            print "\nProcessing input data ..."
        if len(options.output_tag) == 0:
            options.output_tag = input_data
        process_input_data(input_data, input_type, input_label, options)
    else:
        parser.print_usage()
        sys.exit(2)
