#!/usr/bin/env python
# coding: utf-8

#################################################################################################################
#Filename: create_deepcell_dir_format_from_single_channel_fovs.py

# GOAL: Take input fom tile_czi.py to format tiffs for deepcell segmentation. It is done in place, so the same
#directory is used for deepcell

# INPUT: Path to directory containing single channel images. This is the output from tile_czi (output of tile_czi.py above)

# OUTPUT: Directory structure that is deepcell segmentation friendly. It can be used as data directly

#Output looks like:
# cell_to_segment/
#   deepcell_output/
#   input_data/
#       deepcell_input/
#       mibtiff_inputs/
#       single_channel_inputs/
#           fov0/
#               TIFs/
#                   channel0.tiff
#                   channel1.tiff
#           fov1/
#               TIFs/
#                   channel0.tiff
#                   channel1.tiff

# usage: tile_czi.py [-h] [--debug] [--channel] files [files ...]
#
# positional arguments:
#   files          Input 1 or more CZI File Paths to Tile
#

# optional arguments:
#   -h, --help     show this help message and exit
#   --debug, -d    Prints out information when tiling the czi. (Default is False)

#################################################################################################################
from pathlib import Path

import os
import pathlib
import re
import shutil

from segmentation_utils import print_colored, cell_segment_parser

DEBUG = False

cell_segment_parser.add_argument("--ark_target", "-t",
                    dest='ark_target',
                    action="store",
                    nargs='?', #Show it is optional
                    type=pathlib.Path,
                    help="Specify where the output of script should go to. Should be ark_analysis/data");

cell_segment_parser_args = cell_segment_parser.parse_args()

if cell_segment_parser_args.debug:
    DEBUG = True


def isempty(dir_path):
    # Checks a directory to see if it contains any files

    if len(os.listdir(dir_path)) == 0: #if is empty
        return True
    return False

def create_dir(dir_path):
    # Create an empty directory at dir_path if it doesn't yet exist
    try:
        Path.mkdir(dir_path)
    except FileExistsError:
        if not isempty(dir_path):
            print_colored("red", f'Directory {dir_path.name} exists and is not empty. Not creating')

def format_directory(directory_of_formatted_fovs):
    for formatted_dir in directory_of_formatted_fovs:
        fov_dirs = []
        for fov in Path(formatted_dir).iterdir(): #Loop thru each fov dir
            if fov.is_dir() and re.match("fov\d+", fov.name): #if it is a dir and it is the fov# format, then add it to
                fov_dirs.append(fov)

        # Create required directory structure for DeepCell
        deepcell_output_dir = Path(formatted_dir, "deepcell_output")
        create_dir(deepcell_output_dir)
        if DEBUG: print(f"DEBUG: Created {deepcell_output_dir}")

        input_dir = Path(formatted_dir, "input_data")
        create_dir(input_dir)
        if DEBUG: print(f"DEBUG: Created {input_dir}")

        single_tiff_dir = Path(input_dir, "single_channel_inputs")
        create_dir(single_tiff_dir)
        if DEBUG: print(f"DEBUG: Created {single_tiff_dir}")

        mibitiff_dir = Path(input_dir, "mibitiff_inputs")
        create_dir(mibitiff_dir)
        if DEBUG: print(f"DEBUG: Created {mibitiff_dir}")

        deepcell_input_dir = Path(input_dir, "deepcell_input")
        create_dir(deepcell_input_dir)
        if DEBUG: print(f"DEBUG: Created {deepcell_input_dir}")

        for fov_dir in fov_dirs:
            print_colored("cyan", f"Rearranging directory {fov_dir}")
            fov_tiff_dir = Path(single_tiff_dir, fov_dir.name)
            create_dir(fov_tiff_dir)
            if DEBUG: print(f"DEBUG: Created {fov_tiff_dir}")

            tiff_dir = Path(single_tiff_dir, fov_dir.name, "TIFs")
            create_dir(tiff_dir)
            if DEBUG: print(f"DEBUG: Created {tiff_dir}")

            tiffs = []
            for im_file in Path(fov_dir).iterdir():
                if (im_file.name.endswith(".tif") or im_file.name.endswith(".tiff")):
                    tiffs.append(im_file)

            for tf in tiffs:
                if DEBUG: print(f"DEBUG: Moving {str(tf)} to {str(Path(tiff_dir, tf.name))}")
                shutil.move(str(tf), str(Path(tiff_dir, tf.name)))

            if DEBUG: print(f"DEBUG: Removing directory {fov_dir}")
            os.rmdir(fov_dir)

        #Handling the -t target file path if specified with argparse
        if(cell_segment_parser_args.ark_target != None):
            if(os.path.exists(cell_segment_parser_args.ark_target)):
                shutil.move(str(formatted_dir), str(cell_segment_parser_args.ark_target))
                print_colored("green", f"Moved {formatted_dir} to {str(cell_segment_parser_args.ark_target)}")
            else:
                print_colored("red", f"{cell_segment_parser_args.ark_target} is not a valid path! Output will go to current directory")
                print_colored("green", f"Run:\n mv -v {formatted_dir}\n to wherever ark-analysis/data file location is")
        else:
            print_colored("green", f"Run:\n mv -v {formatted_dir}\n to wherever ark-analysis/data file location is")

input_files = cell_segment_parser_args.files

directory_of_formatted_fovs_to_rearrange = [Path(input_path) for input_path in input_files]

if __name__ == "__main__":
    format_directory(directory_of_formatted_fovs_to_rearrange)

#Can run: mv  -v cell_to_segment_dir LOCATION_OF_ARK_ANALYSIS/DATA
# to move folder and all contents to deepcell location