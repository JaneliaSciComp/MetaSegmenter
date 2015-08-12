#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
import os, sys

ms_home = os.environ['MS_HOME'] # where source code is located
ms_data = os.environ['MS_DATA'] # dir for initial input data and final output
ms_temp = os.environ['MS_TEMP'] # dir for intermediate data files
rh_home = os.environ['RHOANA_HOME'] # where source code is located

# ----------------------------------------------------------------------

def Segmentation2D_command_line_parser(parser):
    parser.add_option("-A", "--project",dest="project_code",help="code to be used with qsub",metavar="project_code", default="flyTEM")
    parser.add_option("-c", "--compile",action="store_true",dest="compile",help="compile Matlab file",metavar="compile",default=False)
    parser.add_option("-C", "--compile_all",action="store_true",dest="compile_all", help="compile all Matlab code", metavar="compile_all", default=False)
    parser.add_option("-D", "--debug",dest="debug",help="don't delete intermediate outputs", action="store_true", default=False)
    parser.add_option("-e", "--executable",dest="executable",help="executable",metavar="executable",default="MS_S2D_Segmentation2D")
    parser.add_option("-f", "--fracBlack", dest="fracBlack", help="fracBlack for detecting neurons (default=automatic)",metavar="fracBlack", default=None)
    parser.add_option("-F", "--fracBlack2",dest="fracBlack2",help="fracBlack for detect. dark str.",metavar="fracBlack2",default=None)
    parser.add_option("-i","--uint_type",dest="uint",help="8,16,32 or 64",metavar="uint",default="16")
    parser.add_option("-l", "--nlen",   dest="nlen", help="# of subsections for processing a fragment in y (length) direction",metavar="nlen", default="1")
    parser.add_option("-m", "--maxsize",dest="msize", help="# of subsections for processing a fragment in y (length) direction",metavar="nlen",default=sys.maxint)
    parser.add_option("-n", "--node",   dest="node", help="id of the cluster node to be used", metavar="node",  default=0)
    parser.add_option("-o", "--output_folder",dest="output_folder",help="output folder",metavar="output_folder",default=ms_data)
    parser.add_option("-p", "--processing_start",dest="processing_start",help="start processing from segm(=1),xmerg(=2),ymerg(=3),or epilog(=4) ",metavar="processing_start",default=1)
    parser.add_option("-P", "--processing_end",dest="processing_end",help="complete processing at step segm(=1),xmerg(=2),ymerg(=3),or epilog(=4) ",metavar="processing_end",default=4)
    parser.add_option("-r", "--resize",dest="resize_scale", help="scale for resizing BW image", metavar="resize_scale", default=1)
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

def ExtractInputData_command_line_parser(parser):
    parser.add_option("-d", "--dataset", dest="dataset", help="dataset in DVID store",default='test1')
    parser.add_option("-D", "--debug",dest="debug", help="debugging; don't delete shell scripts", action="store_true", default=False)
    parser.add_option("-i", "--uuid",dest="uuid",help="uuid for DVID store",default='fe7')
    parser.add_option("-L", "--length", dest="ydim", help="y-size (length) of an image", metavar="ydim", default="")
    parser.add_option("-H", "--height", dest="zdim", help="z-size (height, or # layers) of an image stack", metavar="zdim", default="")
    parser.add_option("-n", "--node",   dest="node", help="id of the cluster node to be used", metavar="node",  default="")
    parser.add_option("-o", "--output_path",dest="output_path",help="output path",metavar="output_path",default="")
    parser.add_option("-r", "--uri",dest="uri",help="uri for DVID store",default='/api/repo/fe7/instance')
    parser.add_option("-U", "--unprocessed",action="store_true",dest="unprocessed", help="reprocess only the data for which an output file does not exist",default=False)
    parser.add_option("-v", "--verbose",action="store_true",dest="verbose",help="increase the verbosity level of output",default=False)
    parser.add_option("-W", "--width", dest="xdim", help="x-size (width) of an image", metavar="xdim", default="")
    parser.add_option("-X", "--nx",  dest="nx",  help="# of subsections in x direction", metavar="nx", default=1)
    parser.add_option("-Y", "--ny",  dest="ny",  help="# of subsections in y direction", metavar="ny", default=1)
    parser.add_option("-x", "--dx",  dest="dx",  help="# of scans for image overlap in x direction", metavar="dx", default=50)
    parser.add_option("-y", "--dy",  dest="dy",  help="# of scans for image overlap in y direction", metavar="dy", default=50)
    parser.add_option("-z", "--zmin",dest="zmin",help="# of subsections in z direction", metavar="zmin", default=0)
    parser.add_option("-Z", "--zmax",dest="zmax",help="# of subsections in z direction", metavar="zmax", default=sys.maxint)
    return parser

# -----------------------------------------------------------------------

def MergeImageFragments_command_line_parser(parser):
    parser.add_option("-c", "--compile",action="store_true",dest="compile",help="compile Matlab code",metavar="compile",default=False)
    parser.add_option("-d", "--direction", dest="dir", help="merge direction: x or y", metavar="dir", default="")
    parser.add_option("-D", "--debug",dest="debug",help="debugging mode; don't delete shell scripts",action="store_true",default=False)
    parser.add_option("-e", "--executable", dest="executable", help="name of matlab executable", metavar="executable", default="MS_S2D_MergeImageFragments")
    parser.add_option("-i","--uint_type",dest="uint",help="8,16,32 or 64",metavar="uint",default="16")
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

# ----------------------------------------------------------------------

def Fusion3D_command_line_parser(parser):
    parser.add_option("-A", "--project",dest="project_code",help="code to be used with qsub",metavar="project_code", default="flyTEM")
    parser.add_option("-a", "--overlap_area",dest="overlap_area",help="min overlap_area",metavar="overlap_area",default="200")
    parser.add_option("-D", "--debug",dest="debug",help="don't delete intermediate outputs", action="store_true", default=False)
    parser.add_option("-e", "--executable",dest="executable",help="executable",metavar="executable",default="MS_S2D_Segmentation2D")
    parser.add_option("-f", "--overlap_fraction",dest="overlap_fraction",help="overlap_fraction",metavar="overlap_fraction",default="0.5")
    parser.add_option("-i","--uint_type",dest="uint",help="8,16,32 or 64",metavar="uint",default="64")
    parser.add_option("-l", "--nlen",   dest="nlen", help="# of subsections for processing a fragment in y (length) direction",metavar="nlen", default="1")
    parser.add_option("-m", "--maxsize",dest="msize", help="# of subsections for processing a fragment in y (length) direction",metavar="nlen",default=sys.maxint)
    parser.add_option("-n", "--node",   dest="node", help="id of the cluster node to be used", metavar="node",  default=0)
    parser.add_option("-o", "--output_folder",dest="output_folder",help="output folder",metavar="output_folder",default=ms_data)
    parser.add_option("-p", "--processing_step_beg",dest="processing_step_beg",help="start processing from matr(1), merge(2), trav(3), relab(4) or epil(5)",metavar="processing_step_beg",default=1)
    parser.add_option("-P", "--processing_step_end",dest="processing_step_end",help="end processing with matr(1), merge(2), trav(3), relab(4) or epil(5)",metavar="processing_step_end",default=5)
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

def GenerateMatrices_command_line_parser(parser):
    parser.add_option("-a", "--overlap_area",dest="overlap_area",help="min overlap_area",metavar="overlap_area",default="200")
    parser.add_option("-d", "--dataset", dest="dataset", help="dataset in DVID store",default='test1')
    parser.add_option("-D", "--debug",dest="debug", help="debugging; don't delete shell scripts", action="store_true", default=False)
    parser.add_option("-f", "--overlap_fraction",dest="overlap_fraction",help="overlap_fraction",metavar="overlap_fraction",default="0.3")
    parser.add_option("-i", "--uuid",dest="uuid",help="uuid for DVID store",default='fe7')
    parser.add_option("-L", "--length", dest="ydim", help="y-size (length) of an image", metavar="ydim", default="")
    parser.add_option("-H", "--height", dest="zdim", help="z-size (height, or # layers) of an image stack", metavar="zdim", default="")
    parser.add_option("-n", "--node",   dest="node", help="id of the cluster node to be used", metavar="node",  default=0)
    parser.add_option("-o", "--output_path",dest="output_path",help="output path",metavar="output_path",default="")
    parser.add_option("-r", "--uri",dest="uri",help="uri for DVID store",default='/api/repo/fe7/instance')
    parser.add_option("-U", "--unprocessed",action="store_true",dest="unprocessed", help="reprocess only the data for which an output file does not exist",default=False)
    parser.add_option("-v", "--verbose",action="store_true",dest="verbose",help="increase the verbosity level of output",default=False)
    parser.add_option("-W", "--width", dest="xdim", help="x-size (width) of an image", metavar="xdim", default="")
    parser.add_option("-X", "--nx",  dest="nx",  help="# of subsections in x direction (default=1)", metavar="nx", default=1)
    parser.add_option("-Y", "--ny",  dest="ny",  help="# of subsections in y direction (default=1)", metavar="ny", default=1)
    parser.add_option("-x", "--dx",  dest="dx",  help="# of scans for image overlap in x direction", metavar="dx", default=50)
    parser.add_option("-y", "--dy",  dest="dy",  help="# of scans for image overlap in y direction", metavar="dy", default=50)
    parser.add_option("-z", "--zmin",dest="zmin",help="# of subsections in z direction", metavar="zmin", default=0)
    parser.add_option("-Z", "--zmax",dest="zmax",help="# of subsections in z direction", metavar="zmax", default=sys.maxint)
    return parser

# -----------------------------------------------------------------------

def MergeMatrices_command_line_parser(parser):
    parser.add_option("-a", "--overlap_area",dest="overlap_area",help="min overlap_area",metavar="overlap_area",default="200")
    parser.add_option("-D", "--debug",dest="debug", help="debugging; don't delete shell scripts", action="store_true", default=False)
    parser.add_option("-n", "--node",   dest="node", help="id of the cluster node to be used", metavar="node",  default=0)
    parser.add_option("-o", "--output_path",dest="output_path",help="output path",metavar="output_path",default="")
    parser.add_option("-v", "--verbose",action="store_true",dest="verbose",help="increase the verbosity level of output",default=False)
    parser.add_option("-X", "--nx",  dest="nx",  help="# of subsections in x direction", metavar="nx", default=1)
    parser.add_option("-Y", "--ny",  dest="ny",  help="# of subsections in y direction", metavar="ny", default=1)
    parser.add_option("-z", "--zmin",dest="zmin",help="# of subsections in z direction", metavar="zmin", default=0)
    parser.add_option("-Z", "--zmax",dest="zmax",help="# of subsections in z direction", metavar="zmax", default=sys.maxint)
    return parser

# -----------------------------------------------------------------------

def TraverseFusionTrees_command_line_parser(parser):
    parser.add_option("-D", "--debug",dest="debug", help="debugging; don't delete shell scripts", action="store_true", default=False)
    parser.add_option("-o", "--output_path",dest="output_path",help="output path",metavar="output_path",default="")
    parser.add_option("-v", "--verbose",action="store_true",dest="verbose",help="increase the verbosity level of output",default=False)
    parser.add_option("-X", "--nx",  dest="nx",  help="# of subsections in x direction", metavar="nx", default=1)
    parser.add_option("-Y", "--ny",  dest="ny",  help="# of subsections in y direction", metavar="ny", default=1)
    parser.add_option("-z", "--zmin",dest="zmin",help="# of subsections in z direction", metavar="zmin", default=0)
    parser.add_option("-Z", "--zmax",dest="zmax",help="# of subsections in z direction", metavar="zmax", default=sys.maxint)
    return parser

# -----------------------------------------------------------------------

def RelabelSegmentedData_command_line_parser(parser):
    parser.add_option("-D", "--debug",dest="debug", help="debugging; don't delete shell scripts", action="store_true", default=False)
    parser.add_option("-n", "--node",   dest="node", help="id of the cluster node to be used", metavar="node",  default=0)
    parser.add_option("-o", "--output_path",dest="output_path",help="output path",metavar="output_path",default="")
    parser.add_option("-v", "--verbose",action="store_true",dest="verbose",help="increase the verbosity level of output",default=False)
    parser.add_option("-z", "--zmin",dest="zmin",help="# of subsections in z direction", metavar="zmin", default=0)
    parser.add_option("-Z", "--zmax",dest="zmax",help="# of subsections in z direction", metavar="zmax", default=sys.maxint)
    return parser

# ----------------------------------------------------------------------

def CreateH5Stack_command_line_parser(parser):
    parser.add_option("-a","--alpha", dest="alpha", help="smoothing coefficient",  metavar="alpha", default="0.01")
    parser.add_option("-c","--chunked",action="store_true",dest="chunked",help="each layer is chunk",metavar="chunked",default=False)
    parser.add_option("-i","--uint_type",dest="uint",help="8,16,32 or 64",metavar="uint",default="8")
    parser.add_option("-m","--match_str", dest="match_string",help="match str. for input file names", metavar="mstr",default="")
    parser.add_option("-o","--output_name",dest="output_name",help="name of the output HDF5 file", metavar="out", default="output.h5")
    parser.add_option("-t","--type", dest="output_type", help="output type (='data','labels' or 'mask')", metavar="ot", default="")
    parser.add_option("-v","--verbose",action="store_true",dest="verbose",help="increase verbosity of output", default=False)
    parser.add_option("-u","--unmatch_str",dest="unmatch_string",help="unmatch str. for input file names",metavar="unmstr",default="")
    parser.add_option("-z", "--zmin",dest="zmin",help="min z-layer to be processed", metavar="zmin", default=0)
    parser.add_option("-Z", "--zmax",dest="zmax",help="max z-layer to be processed", metavar="zmax", default=sys.maxint)
    return parser

# ----------------------------------------------------------------------

def RhoanaSegmentation2D_command_line_parser(parser):
    parser.add_option("-A", "--project",dest="project_code",help="code to be used with qsub",metavar="project_code", default="flyTEM")
    parser.add_option("-c", "--classifier",dest="classifier",help="classifier file",metavar="classifier", \
                                           default=os.path.join(rh_home, "rhoana", "ClassifyMembranes", "GB_classifier.txt"))
    parser.add_option("-D", "--debug",dest="debug",help="don't delete intermediate outputs", action="store_true", default=False)
    parser.add_option("-e", "--executable",dest="executable",help="executable",metavar="executable",\
                                           default=os.path.join(rh_home, "rhoana", "ClassifyMembranes", "classify_image"))    
    parser.add_option("-i", "--inverse_probabilities", action="store_true",dest="inverse_probabilities",help="output membrane probabilities", default=False)
    parser.add_option("-m", "--memprob",dest="memprobs",help="membrane probabilities dir name",metavar="membprobs",default="membranes")
    parser.add_option("-o", "--output_folder",dest="output_folder",help="output folder",metavar="output_folder",default=ms_data)
    parser.add_option("-p", "--processing_start",dest="processing_start",help="start processing from probs(=1) or segm(=2)",\
                      metavar="processing_start",default=1)
    parser.add_option("-P", "--processing_end",dest="processing_end",help="complete processing at step probs(=1) or segm(=2) ",\
                      metavar="processing_end",default=2)
    parser.add_option("-t", "--output_tag",dest="output_tag",help="output tag, default=input_data",metavar="output_tag",default="")
    parser.add_option("-v", "--verbose",action="store_true",dest="verbose",help="increase the verbosity level of output",default=False)
    parser.add_option("-z", "--zmin",dest="zmin",help="min z-layer to be processed", metavar="zmin", default=0)
    parser.add_option("-Z", "--zmax",dest="zmax",help="max z-layer to be processed", metavar="zmax", default=sys.maxint)
    return parser

# ----------------------------------------------------------------------

def RhoanaTrainClassifier_command_line_parser(parser):
    parser.add_option("-A", "--project",dest="project_code",help="code to be used with qsub",metavar="project_code", default="flyTEM")
    parser.add_option("-c", "--classifier_type",dest="classifier_type",help="GB(=gradient boosting) or RF(=random forest)",metavar="classifier_type",default="GB")
    parser.add_option("-D", "--debug",dest="debug",help="don't delete intermediate outputs", action="store_true", default=False)
    parser.add_option("-f", "--features",dest="features",help="name of a folder containing computed features",metavar="features",
                                           default=os.path.join(ms_data,"features_150324_pedunculus"))
    parser.add_option("-M", "--membrane_labels",dest="membrane_labels",help="folder containing membrane training labels",\
                                                metavar="membrane_labels",\
                                                default=os.path.join(ms_data,"membrane_labels_150324_pedunculus"))
    parser.add_option("-m", "--mitochondria_labels",dest="mitochondria_labels",help="folder containing mitochondria training labels",\
                                                metavar="mitochondria_labels",\
                                                default=os.path.join(ms_data,"mitochonria_labels_150324_pedunculus"))
    parser.add_option("-o", "--output_file",dest="output_file",help="output file",metavar="output_file",
                                           default=os.path.join(ms_home, "RhoanaSegmentation", "GB_classifier.txt"))
    parser.add_option("-p", "--processing_start",dest="processing_start",help="start processing from probs(=1) or segm(=2)",\
                      metavar="processing_start",default=1)
    parser.add_option("-P", "--processing_end",dest="processing_end",help="complete processing at step probs(=1) or segm(=2) ",\
                      metavar="processing_end",default=3)
    parser.add_option("-r", "--raw_images",dest="raw_images",help="folder containing the raw images for training",\
                                           metavar="raw_images", default=os.path.join(ms_data,"raw_150324_pedunculus"))
    parser.add_option("-v", "--verbose",action="store_true",dest="verbose",help="increase the verbosity level of output",default=False)
    parser.add_option("-z", "--zmin",dest="zmin",help="min z-layer to be processed", metavar="zmin", default=0)
    parser.add_option("-Z", "--zmax",dest="zmax",help="max z-layer to be processed", metavar="zmax", default=sys.maxint)
    return parser

