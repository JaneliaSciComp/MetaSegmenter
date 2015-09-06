#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#
# Purpose: perform 2D segmentation using GALA
#          then post-process results
#
# ------------------------- imports -------------------------

#from __future__ import absolute_import
#from __future__ import print_function
# imports
import os, sys, re, optparse
from gala import imio, classify, features, agglo, evaluate as ev
from six.moves import map
import numpy as np

import MS_LIB_Dict
import MS_LIB_Options
import MS_LIB_IO

from PIL import Image

ms_home = os.environ['MS_HOME']   # where source code is located
ms_data = os.environ['MS_DATA']   # dir for initial input data and final output
ms_temp = os.environ['MS_TEMP']   # dir for intermediate data files
gala_home    = os.environ['GALA_HOME']    # where source code is located
ilastik_home = os.environ['ILASTIK_HOME'] # where source code is located

# -----------------------------------------------------------------------------

def perform_segmentation(training_raw, groundtruth_annot_dir, \
                         groundtruth_seg_dir, production_raw, \
                         options):
    raw_files = os.listdir(os.path.join(ms_data, training_raw))
    labels    = os.listdir(os.path.join(ms_data, groundtruth_annot_dir))

    # Assume only one raw and one labels file for now
    command = "MS_UT_CreateH5LabelsForIlastik.py " + groundtruth_annot_dir + \
              " groundtruth.h5"
    print "\nRunning the commmand: " , command
    os.system(command)

    # Generate Ilastik project file
    command = "python " + os.path.join(ilastik_home, "ilastik", "bin", "train_headless.py") + \
              " training.ilp '" + training_raw +  "/*.png' groundtruth.h5/main"
    print "\nRunning the commmand: " , command
    os.system(command)

    # Produce the probabilities and oversegmentation files from the inputs
    if not os.path.isdir('seg_data'):
        print "\nRunning the commmand: mkdir seg_data"
        os.mkdir('seg_data')
    else:
        print "\nRunning the commmand: rm -f seg_data/*"
        os.system("rm -rf seg_data/*")
    if os.path.isfile('STACKED_prediction.h5'):
        command = "rm -f STACKED_prediction.h5"
        print "\nRunning the commmand: " , command
        os.system(command)
    if not os.path.islink('gala-segmentation-pipeline'):
        os.symlink('gala-segmentation-pipeline', os.path.join(gala_home, 'gala',\
                   'bin', 'gala-segmentation-pipeline'))
    command = "python gala-segmentation-pipeline " + \
              " . -I '" + training_raw + "/*.png' --ilp-file training.ilp " + \
              " --enable-gen-supervoxels --enable-gen-agglomeration " + \
              " --enable-gen-pixel --seed-size 5   --segmentation-thresholds 0.0"
    print "\nRunning the commmand: " , command
    os.system(command)

    # Put probabilities in the right format:
    command = "MS_UT_ProbsNF2GALA.py STACKED_prediction.h5 probabilities_training.h5"
    print "\nRunning the commmand: ", command
    os.system(command)

    # Rename the oversegmentation training file
    for file in os.listdir("seg_data"):
        if re.search(".h5", file):
            oversegmentation_file = os.path.join("seg_data", file)
    command = "mv " + oversegmentation_file + " oversegmentation_training.h5"
    print "\nRunning the commmand: ", command
    os.system(command)

    # Create the groundtruuth segmentation labels file
    command = "MS_UT_CreateH5Groundtruth.py " + groundtruth_seg_dir + " seg_labels_file.h5"
    print "\nRunning the commmand: " , command
    os.system(command)

    # Read in training data
    print "\nReading training data into gala..."
    gt_train, pr_train, ws_train = (map(imio.read_h5_stack,
                                    ['seg_labels_file.h5', 'probabilities_training.h5',
                                     'oversegmentation_training.h5']))

    # create a feature manager
    fm = features.moments.Manager()
    fh = features.histogram.Manager()
    fc = features.base.Composite(children=[fm, fh])

    # create graph and obtain a training dataset
    g_train = agglo.Rag(ws_train, pr_train, feature_manager=fc)
    (X, y, w, merges) = g_train.learn_agglomerate(gt_train, fc)[0]
    y = y[:, 0] # gala has 3 truth labeling schemes, pick the first one
    print((X.shape, y.shape)) # standard scikit-learn input format

    # train a classifier, scikit-learn syntax
    rf = classify.DefaultRandomForest().fit(X, y)
    # a policy is the composition of a feature map and a classifier
    learned_policy = agglo.classifier_probability(fc, rf)

    if not os.path.isdir('seg_data'):
        print "\nRunning the commmand: mkdir seg_data"
        os.mkdir('seg_data')
    else:
        command = "rm -f seg_data/*"
        print "\nRunning the commmand: ", command
        os.system(command)
    if os.path.isfile('STACKED_prediction.h5'):
        command = "rm -f STACKED_prediction.h5"
        print "\nRunning the commmand: " , command
        os.system(command)

    # Produce the probabilities and oversegmentation files from the inputs
    command = "python " + os.path.join(gala_home, "gala", "bin", "gala-segmentation-pipeline") + \
              " . -I '" + productiobn_raw + "/*.png' --ilp-file training.ilp " + \
              " --enable-gen-supervoxels --enable-gen-agglomeration " + \
              " --enable-gen-pixel --seed-size 5   --segmentation-thresholds 0.0"
    print "\nRunning the commmand: ", command
    os.system(command)

    # Put probabilities in the right format:
    if os.path.isdir("seg_data"):
        print "\nRunning the commmand: rm -rf seg_data"
        os.system("rm -rf seg_data")
    command = "MS_UT_ProbsNF2GALA.py STACKED_prediction.h5 probabilities_production.h5"
    print "\nRunning the commmand: ", command
    os.system(command)

    # Rename the oversegmentation production file
    for file in os.listdir("seg_data"):
        if re.search(".h5", file):
            oversegmentation_file = os.path.join("seg_data", file)
    command = "mv " + oversegmentation_file + " oversegmentation_production.h5"

    # get the production data and make a RAG with the trained policy
    pr_test, ws_test = (map(imio.read_h5_stack,
                        ['probabilities_production.h5', 'oversegmentation_production.h5']))
    g_test = agglo.Rag(ws_test, pr_test, learned_policy, feature_manager=fc)
    g_test.agglomerate(0.5) # best expected segmentation
    seg_production = g_test.get_segmentation()

    # Output the segmentation file
    img = Image.fromarray(seg_production)
    img.save('segmentation_production.h5')


# -----------------------------------------------------------------------------

def output_results(ws_test, gt_test, seg_test1):
    results = np.vstack((
        ev.split_vi(ws_test, gt_test),
        ev.split_vi(seg_test1, gt_test),
        ))
    print(results)

# -----------------------------------------------------------------------------

if __name__ == "__main__":

    usage = "Usage: \n\
    %prog <raw_train_data> <gt_labels> <raw_data> [options (-h to list)]"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = MS_LIB_Options.GalaSegmentation2D_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if len(args) == 4:
        training_raw    = args[0] 
        gt_annot_labels = args[1]
        gt_seg_labels   = args[2]
        production_raw  = args[3] 
        perform_segmentation(training_raw, gt_annot_labels, gt_seg_labels, \
                             production_raw, options)
    else:
        parser.print_usage()
        sys.exit(2)
 
