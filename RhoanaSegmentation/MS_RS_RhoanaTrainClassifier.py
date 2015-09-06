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

def compute_features(data_tag, options):
    file_list = sorted(os.listdir(options.raw_dir))
    z = int(options.zmin)
    Z = int(options.zmax)+1
    file_list = file_list[z:Z]
    executable = os.path.join(rh_home, "rhoana", 
        "ClassifyMembranes", "compute_features")
    for i in range(0, len(file_list)):
        file_h5 = file_list[i].split(".")[0] + ".h5"
        command = executable + " " + \
                  os.path.join(ms_data, options.raw_dir, file_list[i]) + " " +\
                  os.path.join(ms_data, options.features_dir, file_h5)       
        if options.verbose:
            print "\nComputing features for ", file_list[i]
        os.system(command)
    return 

# -----------------------------------------------------------------------------

def produce_training_images(data_tag, options):
    file_list = sorted(os.listdir(options.raw_dir))
    z = int(options.zmin)
    Z = int(options.zmax)+1
    file_list = file_list[z:Z]
    for i in range(0, len(file_list)):
        if options.verbose:
            print "\nCreating training images for ", file_list[i]
        raw_file        = os.path.join(options.raw_dir, file_list[i])
        training_file   = os.path.join(options.training_dir, file_list[i])
        pos_labels_file = os.path.join(options.pos_labels_dir, \
                                       file_list[i])
        command_train   = os.path.join(ms_home, "Utilities", \
                             "MS_UT_CreateTrainingImage.py") + " " + raw_file +\
                             " " + pos_labels_file + " " + training_file      
        if len(options.neg_labels_dir) > 0:
            neg_labels_file = os.path.join(options.neg_labels_dir, \
                                       file_list[i])
            command_train += " " + neg_labels_file
        if options.verbose:
            print "\n...Running the command:", command_train
        os.system(command_train)
    return

# -----------------------------------------------------------------------------

def perform_classifier_training(data_tag, options):
    file_list = sorted(os.listdir(options.raw_dir))
    z = int(options.zmin)
    Z = int(options.zmax)+1
    file_list = file_list[z:Z]
    executable = os.path.join(rh_home, "rhoana", "ClassifyMembranes", "train_gb")
    if options.classifier_type == "RF":
        executable = os.path.join(rh_home, "rhoana", "ClassifyMembranes", \
                                  "train_randomforest")
    command_train = executable
    for i in range(0, len(file_list)):
        file    = file_list[i]
        file_h5 = file.split(".")[0] + ".h5"
        training_file = os.path.join(ms_data, options.training_dir, file)
        features_file = os.path.join(ms_data, options.features_dir, file_h5)
        command_train += " " + training_file + " " + features_file
    classifier_file =  os.path.join(ms_home, "RhoanaSegmentation", \
                       options.classifier_type + "_classifier_" + options.training_type +\
                       "_" + data_tag + ".txt")
    command_train += " " + classifier_file  
    if options.verbose:
        print "\n...Training ", options.classifier_type, " classifier on membrane examples:"
        print "   Output classifier file=", classifier_file
        print "   Command_train=", command_train
    os.system(command_train)

    return

# -----------------------------------------------------------------------------

# Process individual fragment segmentation job                 
def train_classifier(data_tag, options):
    if len(options.raw_dir) == 0:
        options.raw_dir     = os.path.join(ms_data,"raw_" + data_tag)

    if len(options.training_dir) == 0:
        options.training_dir = os.path.join(ms_data, options.training_type + \
                                         "_training_" + data_tag)
    if len(options.features_dir) == 0:
        options.features_dir = os.path.join(ms_data, "features_" + data_tag)  

    if int(options.processing_start) <= 1:
        compute_features(data_tag, options)

    if int(options.processing_start) <= 2 and int(options.processing_end)  >= 2:
        produce_training_images(data_tag, options)

    if int(options.processing_start) <= 3 and int(options.processing_end) >= 3:
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
        if not options.training_type in ["memb_vs_rest", "mito_vs_rest", \
                                         "mitomemb_vs_rest", "memb_vs_mito"]:
            print "\nIncorrectly specified training type"
            sys.exit(2)
        train_classifier(data_tag, options)
    else:
        parser.print_usage()
        sys.exit(2)
