#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#
# Purpose: perform 2D segmentation using GALA
#          then post-process results
#
# ------------------------- imports -------------------------

# imports
import os, sys, re, optparse
from gala import imio, classify, features, agglo, evaluate as ev
from six.moves import map
import numpy as np
import h5py
from PIL import Image
from skimage.color import label2rgb
import ntpath

import MS_LIB_Dict
import MS_LIB_Options
import MS_LIB_IO

ms_home = os.environ['MS_HOME']   # where source code is located
ms_data = os.environ['MS_DATA']   # dir for initial input data and final output
ms_temp = os.environ['MS_TEMP']   # dir for intermediate data files
gala_home          = os.environ['GALA_HOME']    # where source code is located
ilastik_home       = os.environ['ILASTIK_HOME'] # where source code is located
ilastik_home2      = os.environ['ILASTIK_HOME2']
neuroproof_2d_home = os.environ['NEUROPROOF_2D_HOME']

# -----------------------------------------------------------------------------

def extract_file_list_from_directory(input_dir, ind_min, ind_max):
    file_list = []
    print "input_dir=", input_dir
    files = sorted(os.listdir(input_dir))
    for i in range(0, len(files)):
        if i >= ind_min and i <= min(ind_max, len(files)):
            file_path = os.path.join(ms_data, input_dir, files[i])
            file_list.append(file_path)
    return file_list

# -----------------------------------------------------------------------------

def create_dir(dir_path):
    if os.path.isdir(dir_path):
        command = "rm -rf " + dir_path + "/*"
        print "\nRunning the commmand: " , command, "\n"
        os.system(command)
    else:
        os.mkdir(dir_path)

# -----------------------------------------------------------------------------

def remove_dir(dir_path):
    command = "rm -rf " + dir_path 
    print "\nRunning the commmand: " , command, "\n"
    os.system(command)

# -----------------------------------------------------------------------------

def create_mosaic_file(image_file_list, mosaic_file):
    list_len = len(image_file_list)
    for i in range(0, list_len):
        data1 = np.asarray(Image.open(image_file_list[i]))
        if i == 0:
            shape1 = data1.shape
            mosaic_data = np.zeros([shape1[0]*list_len, shape1[1]], dtype = data1.dtype)
        mosaic_data[(i*shape1[0]):((i+1)*shape1[0]), :] = data1
    img = Image.fromarray(mosaic_data)
    img.save(mosaic_file)                           
    print "image_file_list=", image_file_list, " mosaic_file=", mosaic_file

# -----------------------------------------------------------------------------

def create_ilastik_trainig_file(raw_file_list, gt_file_list):
    create_dir("tmp_raw_mosaic_dir")
    create_dir("tmp_gt_mosaic_dir")
    create_mosaic_file(raw_file_list, os.path.join("tmp_raw_mosaic_dir", "mosaic_raw.png"))
    create_mosaic_file( gt_file_list, os.path.join("tmp_gt_mosaic_dir",  "mosaic_gt.png"))

    # Create only one raw and one labels file
    command = "MS_UT_Im2H5.py " + os.path.join("tmp_gt_mosaic_dir", "mosaic_gt.png") + \
              " groundtruth.h5"
    print "\nRunning the commmand: " , command, "\n"
    os.system(command)

    # Generate Ilastik project file
    training_ilp = "training.ilp"
    python_script = os.path.join(ilastik_home, "ilastik", "bin", "train_headless.py")
#   python_script = os.path.join(ilastik_home2, "src", "ilastik", "ilastik", "bin", "train_headless.py")
    command = "python " + python_script +\
              " " + training_ilp + " 'tmp_raw_mosaic_dir/*.png' groundtruth.h5/stack"
    print "\nRunning the commmand: " , command, "\n"
    os.system(command)

#   remove_dir("tmp_raw_mosaic_dir")
#   remove_dir("tmp_gt_mosaic_dir")

    return training_ilp         

# -----------------------------------------------------------------------------

def populate_one_file(training_ilp, raw_file, target_dir):
    # Populate temporary dirs
    create_dir("tmp_raw")        
    command = "cp " + raw_file + " " + " tmp_raw"
    print "\nRunning the commmand: " , command, "\n"
    os.system(command)

    # Remove old files, if they exist
    if not os.path.isdir('seg_data'):
        print "\nRunning the commmand: mkdir seg_data"
        os.mkdir('seg_data')
    else:
        print "\nRunning the commmand: rm -f seg_data/*"
        os.system("rm -rf seg_data/*")
    if os.path.isfile('STACKED_prediction.h5'):
        command = "rm -f STACKED_prediction.h5"
        print "\nRunning the commmand: " , command, "\n"
        os.system(command)

    # Produce the probabilities and oversegmentation files from the inputs
    if not os.path.islink('gala-segmentation-pipeline'):
        os.symlink('gala-segmentation-pipeline', os.path.join(gala_home, 'gala',\
                   'bin', 'gala-segmentation-pipeline'))
    command = "python gala-segmentation-pipeline " + \
              " . -I 'tmp_raw/*.png' --ilp-file " + training_ilp + " " + \
              " --enable-gen-supervoxels --enable-gen-agglomeration " + \
              " --enable-gen-pixel --seed-size 5 --segmentation-thresholds 0.0"
    print "\nRunning the commmand: " , command, "\n"
    os.system(command)
#   remove_dir("tmp_raw")

    # Populate the image file
    base_name = ntpath.basename(raw_file).split(".")[0]
    suffix    = ntpath.basename(raw_file).split(".")[1]
    img_file  = os.path.join(target_dir, "img." + base_name[1:] + "." + suffix)
    command = "cp " + raw_file + " " + img_file
    print "\nRunning the commmand: ", command, "\n"
    os.system(command)

    # Populate the predictions file             
    base_name = ntpath.basename(raw_file).split(".")[0]
    print "target_dir=", target_dir, " basename=",  base_name
    
    probs_h5 = os.path.join(target_dir, "img." + base_name[1:] + "_prediction.h5") 
    command = "MS_UT_ProbsNF2GALA.py STACKED_prediction.h5 " + probs_h5
    print "\nRunning the commmand: ", command, "\n"
    os.system(command)

    # Populate watersfed file
    base_name = ntpath.basename(raw_file).split(".")[0]
    watershed_h5 = os.path.join(target_dir, "img." + base_name[1:] + "_watershed.h5")
    command = "python " + os.path.join(neuroproof_2d_home, "generate_watershed_dan.py") +\
              " " + probs_h5 + " " + watershed_h5 + " 0"
    print "\nRunning the commmand: ", command, "\n"
    os.system(command)

# -----------------------------------------------------------------------------

def populate_training_dir(training_raw, groundtruth_labels_dir, options):
    # Create training dir
    create_dir(options.training_dir)

    # Create a list of training files
    raw_file_list = extract_file_list_from_directory(training_raw, \
                                 int(options.tmin), int(options.tmax))
    groundtruth_labels_file_list =  extract_file_list_from_directory(groundtruth_labels_dir,\
                                 int(options.tmin), int(options.tmax))
    if not len(raw_file_list) == len(groundtruth_labels_file_list):
        sys.exit("# files in raw and groundtruth training lists are not the same")
  
    print "groundtruth_labels_file_list=", groundtruth_labels_file_list 
    # Create ilastik training file
    training_ilp = create_ilastik_trainig_file(raw_file_list, groundtruth_labels_file_list)
 
    # Populate training_dir
    print "options.training_dir=", options.training_dir
    for i in range(0, len(raw_file_list)):
        populate_one_file(training_ilp, raw_file_list[i], options.training_dir) 
        base_name = ntpath.basename(raw_file_list[i]).split(".")[0]
        groundtruth_seg_file = os.path.join(options.training_dir, "img." + base_name[1:] + "_groundtruth.h5")
        prediction_file      = os.path.join(options.training_dir, "img." + base_name[1:] + "_prediction.h5")
        command = "python " + os.path.join(neuroproof_2d_home, "generate_seg_groundtruth.py") \
                            + " " + groundtruth_labels_file_list[i] + " " + prediction_file \
                            + " " + groundtruth_seg_file
        print "\nRunning the commmand: ", command, "\n"
        os.system(command)
    return training_ilp

# -----------------------------------------------------------------------------

def populate_production_dir(training_ilp, production_raw, options):
    # Create production dir
    create_dir(options.production_dir)

    # Create a list of training files
    raw_file_list = extract_file_list_from_directory(production_raw, \
                                 int(options.zmin), int(options.zmax))

    # Populate production die
    for i in range(0, len(raw_file_list)):
        populate_one_file(training_ilp, raw_file_list[i], options.production_dir)
    return 

# -----------------------------------------------------------------------------

def learn_classifier(classifier_h5, training_dir, options):
    command = os.path.join(neuroproof_2d_home, "NeuroProof_stack_learn") + \
              " -train_dir " + training_dir + "  -mito_thd 0.35  " + \
              " -classifier " + classifier_h5 + " -iteration 1 -strategy 2 "
    print "\nRunning the commmand: " , command, "\n"
    os.system(command)
    return 

# -----------------------------------------------------------------------------

def perform_segmentation(output_dir, classifier_h5, production_dir, options):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    command = os.path.join(neuroproof_2d_home, "NeuroProof_stack") + \
              " -image_dir " + production_dir + " -mito_thd 0.35  " + \
              " -classifier " + classifier_h5 + " -algorithm 1 -threshold 0.2"
    print "\nRunning the commmand: " , command, "\n"
    os.system(command)
    return 

# -----------------------------------------------------------------------------

def create_segmentation_stack(output_dir, options):
    result_files = []
    print "output_dir=", output_dir
    for file in os.listdir(output_dir):
        if re.search("_m.h5", file):
            result_files.append(os.path.join(output_dir, file))
    result_files = sorted(result_files)
    for i in range(0, len(result_files)):
        f  = h5py.File(result_files[i], 'r')
        data1 = f["stack"]                        
        if i == 0:
            data_stack = np.zeros([len(result_files), data1.shape[1], data1.shape[2]])
            print "data1.shape=", data1.shape
        data_stack[i,:,:] = data1
    if len(result_files) > 0:
        f  = h5py.File(os.path.join(output_dir, "resulting_stack.h5"), 'w') 
        f.create_dataset('stack', data_stack.shape, data = data_stack)

# -----------------------------------------------------------------------------

def process_data(raw_dir, groundtruth_labels_dir, options):
    training_ilp  = "training.ilp" # ilastik training file
    classifier_h5 = "np_classifier.h5"
    output_dir    = os.path.join("production_dir", "output")

    if int(options.processing_start) == 1 and int(options.processing_end) >= 1:
        training_ilp = populate_training_dir(raw_dir, \
                                    groundtruth_labels_dir, options)

    if int(options.processing_start) <= 2 and int(options.processing_end) >= 2:
        populate_production_dir(training_ilp, raw_dir, options)

    if int(options.processing_start) <= 3 and int(options.processing_end) >= 3:
        learn_classifier(classifier_h5, options.training_dir, options)

    if int(options.processing_start) <= 4 and int(options.processing_end) >= 4:
       perform_segmentation(output_dir, classifier_h5, options.production_dir, options)
  
    if int(options.processing_start) <= 5 and int(options.processing_end) >= 5:
        create_segmentation_stack(output_dir, options)

# -----------------------------------------------------------------------------

if __name__ == "__main__":

    usage = "Usage: \n\
    %prog <raw_train_data> <gt_seg_labels> <raw_data> [options (-h to list)]"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = MS_LIB_Options.NeuroProofSegmentation2D_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if len(args) == 2:
        raw_dir   = args[0] 
        gt_labels_dir  = args[1]
        
        process_data(raw_dir, gt_labels_dir, options)
    else:
        parser.print_usage()
        sys.exit(2)
 
