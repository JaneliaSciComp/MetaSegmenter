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
        output_file = probs_file.split(".")[0] + ".png"
        output_dir = os.path.join(ms_data, "seg_" + options.output_tag)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        command = "python " + segmentation_script + " " + probs_file_list[i] +\
                  " " + os.path.join(ms_data, output_file)
        os.system(command)
    return 

# -----------------------------------------------------------------------------

def invert_probabilities(h5_infile, h5_outfile):
    print "h5_outfile=", h5_outfile
    f   = h5py.File(h5_infile,  'r')
    f2  = h5py.File(h5_outfile, 'w')
    keys = f.keys()
    for k in keys:
        if not k == "probabilities":
            f2.create_dataset(k,      data = f[k])
        else: 
            probs  = f["probabilities"]
            inverse_probs = 1. - numpy.array(probs)
            f2.create_dataset("probabilities", data = inverse_probs)
    return 

# -----------------------------------------------------------------------------

def compute_membrane_probabilities(input_data, file_list, options):
    classify_image_executable = options.executable
    classifier_file = options.classifier
    output_file_list = []
    for i in range(0, len(file_list)):
        input_file_path  = file_list[i]
        print "input_file_path=", input_file_path
        input_file       = ntpath.basename(input_file_path)     
        output_file      = input_file.split(".")[0] + ".h5"
        output_dir       = os.path.join(ms_data, "probs_"     + input_data)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        if input_data[0:4] == "raw_":
            output_dir   = os.path.join(ms_data, "probs_"     + input_data[4:])
            output_dir2  = os.path.join(ms_data, "probs_inv_" + input_data[4:])
        output_file_path = os.path.join(output_dir, output_file)
        output_file_list.append(output_file_path)
        command = classify_image_executable + " " + input_file_path + " " + \
                  classifier_file + " " + output_file_path     
        if options.verbose:
            print "...command=", command
        os.system(command) 
        # Inverse probabilities
        if options.inverse_probabilities:
            output_dir2      = os.path.join(ms_data, "probs_inv_" + input_data)
            if input_data[0:4] == "raw_":
               output_dir2  = os.path.join(ms_data, "probs_inv_" + input_data[4:])
            if not os.path.isdir(output_dir2):
                os.mkdir(output_dir2)
            output_file_path2= os.path.join(output_dir2, output_file)
            invert_probabilities(output_file_path, output_file_path2)
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

def compute_features(data_tag, options):
    file_list = sorted(os.listdir(options.raw_images))
    z = int(options.zmin)
    Z = int(options.zmax)+1
    file_list = file_list[z:Z]
    executable = os.path.join(rh_home, "rhoana", 
        "ClassifyMembranes", "compute_features")
    for i in range(0, len(file_list)):
        file_h5 = file_list[i].split(".")[0] + ".h5"
        command = executable + " " + \
                  os.path.join(ms_data, options.raw_images, file_list[i]) + " " +\
                  os.path.join(ms_data, options.features,   file_h5)       
        if options.verbose:
            print "\nComputing features for ", file_list[i]
        os.system(command)
    return 

# -----------------------------------------------------------------------------

def produce_training_images(data_tag, options):
    file_list = sorted(os.listdir(options.raw_images))
    z = int(options.zmin)
    Z = int(options.zmax)+1
    file_list = file_list[z:Z]
    for i in range(0, len(file_list)):
        if options.verbose:
            print "\nCreating training images for ", file_list[i]
        raw_image          = os.path.join(options.raw_images, file_list[i])
        membranes_bw       = os.path.join(options.membrane_labels, file_list[i])
        mitochondria_bw    = os.path.join(options.mitochondria_labels, file_list[i])
        membranes_training = os.path.join(ms_data,
                         "membrane_training_"     + data_tag, file_list[i])
        mitochondria_training = os.path.join(ms_data,
                         "mitochondria_training_" + data_tag, file_list[i])
        command_train_memb = os.path.join(ms_home, "Utilities", \
                             "MS_UT_CreateTrainingImage.py") + " " + raw_image +\
                             " " +  membranes_bw + " " +  membranes_training
        if options.verbose:
            print "\n...Running the command:", command_train_memb
        os.system(command_train_memb)
        command_train_mito = os.path.join(ms_home, "Utilities", \
                             "MS_UT_CreateTrainingImage.py") + " " + raw_image +\
                             " " +  mitochondria_bw + " " +  mitochondria_training
        if options.verbose:
            print "\n...Running the command:", command_train_mito
        os.system(command_train_mito)
    return

# -----------------------------------------------------------------------------

def perform_classifier_training(data_tag, options):
    file_list = sorted(os.listdir(options.raw_images))
    z = int(options.zmin)
    Z = int(options.zmax)+1
    file_list = file_list[z:Z]
    executable = os.path.join(rh_home, "rhoana", "ClassifyMembranes", "train_gb")
    if options.classifier_type == "RF":
        executable = os.path.join(rh_home, "rhoana", "ClassifyMembranes", \
                                  "train_randomforest")
    command_train_memb = executable
    command_train_mito = executable
    for i in range(0, len(file_list)):
        file    = file_list[i]
        file_h5 = file.split(".")[0] + ".h5"
        train_memb = os.path.join(ms_data,
                         "membrane_training_" + data_tag, file)
        train_mito = os.path.join(ms_data,
                         "mitochondria_training_" + data_tag, file)
        features   = os.path.join(ms_data,
                         "features_" + data_tag, file_h5)
        command_train_memb += " " + train_memb + " " + features
        command_train_mito += " " + train_mito + " " + features

    command_train_memb += " " + os.path.join(ms_home, "RhoanaSegmentation", \
                                options.classifier_type + "_classifier_membrates_" +\
                                data_tag + ".txt")
    if options.verbose:
        print "\n...Training ", options.classifier_type, " classifier on membrane examples:"
        print "command_train_memb=", command_train_memb
    os.system(command_train_memb)

    command_train_mito += " " + os.path.join(ms_home, "RhoanaSegmentation", \
                                options.classifier_type + "_classifier_mitochondria_" +\
                                data_tag + ".txt")
    if options.verbose:
        print "\n...Training ", options.classifier_type, " classifier on mitochondria examples"
        print "command_train_mito=", command_train_mito
    os.system(command_train_mito)

    return

# -----------------------------------------------------------------------------

# Process individual fragment segmentation job                 
def train_classifier(data_tag, options):
    default_tag = "150324_pedunculus"
    if not data_tag == default_tag:
        options.raw_images          = os.path.join(ms_data,"raw_" + data_tag)
        options.membrate_labels     = os.path.join(ms_data,\
                                      "membrane_labels_" + data_tag)
        options.mitochondria_labels = os.path.join(ms_data,\
                                      "mitochondria_labels_" + data_tag)
    if options.processing_start <= 1:
        compute_features(data_tag, options)

    if options.processing_start <= 2 and options.processing_end  >= 2:
        produce_training_images(data_tag, options)

    if options.processing_start <= 3 and options.processing_end  >= 3:
        perform_classifier_training(data_tag, options)
    
    return

# -----------------------------------------------------------------------------

if __name__ == "__main__":
   
    usage = "Usage: %prog <data_tag> [options (-h to list)]"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = MS_LIB_Options.RhoanaTrainClassifier_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if len(args) == 1:
        data_tag = args[0]
        train_classifier(data_tag, options)
    else:
        parser.print_usage()
        sys.exit(2)
