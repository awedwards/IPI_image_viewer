#!/usr/bin/env python
# coding: utf-8

#################################################################################################################
#Filename: tile_czi.py

# GOAL: Goal of this script is to take a czi file, and to create a directory of subdirectories containing
#single channel tiffs  to prepare their organization for deep_cell segmntation

# INPUT: .czi file (Example: cell_to_segment.czi)

# OUTPUT: Directory containing subdirectories for each FOV, themselves containing single-channel .tiff images
#Output loks like:
    #Cell_to_segment/
    #   Fov0/
    #       channel0.tiff
    #       channel1.tiff ...
    #   Fov1/
    #       channel0.tiff
    #       channel1.tiff ...

#Arguments that you can use!
# usage: tile_czi.py [-h] [--debug] [--channel] files [files ...]
#
# positional arguments:
#   files          Input 1 or more CZI File Paths to Tile
#
# optional arguments:
#   -h, --help     show this help message and exit
#   --debug, -d    Prints out information when tiling the czi. (Default is False)
#   --channel, -c  View, Select, and Add Channels!

#################################################################################################################

from aicspylibczi import CziFile
from pathlib import Path
from segmentation_utils import print_colored, cell_segment_parser
import json

import numpy as np
import tifffile as tf

import os

DEBUG = False #Used to output additional infomation when the czi is being filed
TILE_SIZE = 2048 # Using this tile size because it is the largest that can be done
DEFAULT_CHANNELs_TO_USE = 1
written = False

#All current channels. Can be added to for easier command line argument parsing.
all_channels = [
    ['DAPI','FoxP3','CD4','CD45','CD8'],
    ['DAPI', 'HLADR', 'CD8', 'CD163', 'CD4', 'XCR1', 'CD3', 'PDL1', 'PanCK'],
    ['DAPI', 'HLADR', 'CD8', 'CD163', 'CD4', 'XCR1', 'CD3', 'PDL1', 'EPCAM']  # only for CRC
]
channels_to_use = all_channels[1]

#Using argument parser to organize the input

cell_segment_parser.add_argument("--channel", "-c", dest='channel', action="store_true",
                    help="View, Select, and Add Channels!")

cell_segment_parser_args = cell_segment_parser.parse_args()

if cell_segment_parser_args.debug:
    DEBUG = True

#Gets proper channels from above for tiling
def get_channel_choice(len_of_channels):
    while True:
        try:
            number = int(input('Chose an option from menu: '))
            if 0 < number <= len_of_channels:
                return number
            else:
                raise ValueError("Not in valid number range!")
        except Exception as  e:
            print_colored("red", f"Invalid channel choice! {e}")

if cell_segment_parser_args.channel:
    print_colored("cyan", f"Select which channels you would like to use. Optionally, use your own. (Note: Current Default is option {DEFAULT_CHANNELs_TO_USE})")

    print_colored("green", f"CHANNEL CHOICES")
    for index, channel in enumerate(all_channels):
        print_colored("green", f"[{index}] {channel}")
    print_colored("green", f"[{len(all_channels)}] Input Custom Channels")
    choices = list(range(0, len(all_channels) + 1))

    num_to_use = get_channel_choice(len(all_channels))

    if(num_to_use < len(all_channels)):
        channels_to_use = all_channels[num_to_use]
    else:
        print_colored("cyan", f"Type each channel out separated by a comma, and custom channels will be created. (No trailing comma!)")
        channels_to_use = input("Give Channels: ").replace(" ", "").split(',')

input_czi_files = cell_segment_parser_args.files

#This is useful for restitching fovs back together when done.
def write_tile_breakdown(rows, cols, czi_filename):
    global written
    written = True
    data = {
        "dims": {'rows': rows, 'cols': cols},
        "filename": czi_filename.stem
    }

    print_colored("yellow", f"Writing {str(data)} to .{czi_filename.stem} in final_data/")
    with open(os.path.dirname(os.path.realpath(__file__)) + f'/final_data/.{czi_filename.stem}_metadata.json', 'w') as f:
        json.dump(data, f, ensure_ascii=False)

#Steps
# 1. loop through all czi_files in input
# 2. Loop thru each channel
# 3. Create FOV Directories if not already made
# 4. For the czi, create fov tiles with proper size
# 5. Create the tiff files

def tile_czi_file():
    for input_czi_file in input_czi_files:


        czi_file_path = Path(input_czi_file)
        czi = CziFile(czi_file_path)
        im_shape = czi.get_dims_shape()

        if DEBUG: print("DEBUG: CZI Image Shape:", im_shape)#  Ex: [{'X': (0, 2252), 'Y': (0, 1208), 'C': (0, 9), 'T': (0, 1), 'M': (0, 17), 'S': (0, 1), 'H': (0, 1)}]

        nchannels = im_shape[0]['C'][1]

        assert len(channels_to_use) == nchannels, "Number of channels for CZI and in file must match"

        dir_to_create = Path(os.path.curdir, czi_file_path.stem + '_dir')

        try:
            if DEBUG: print("DEBUG: Creating", dir_to_create)
            Path.mkdir(dir_to_create)
        except FileExistsError:
            print_colored("yellow", f"NOTE: Tried to create {dir_to_create}, but directory {dir_to_create} is already made!")

        with open(Path(dir_to_create, f"tile_metadata.txt"), "w") as f:
            f.write('fov,x1,x2,y1,y2\n')

            if DEBUG: print("DEBUG: Reading", czi_file_path.name)
            for channel in np.arange(nchannels):
                print_colored("cyan", "Reading " + channels_to_use[channel] + " channel")
                fov = 0
                im = czi.read_mosaic(C=channel)
                _, w, h = im.shape #Ex: (1, 7290, 4131)

                #Here we get breakdown with tile_sizes
                rows = list(np.arange(0, w, TILE_SIZE)) #Ex: [0, 2048, 4096, 6144]
                cols = list(np.arange(0, h, TILE_SIZE)) #Ex: [0, 2048, 4096]

                if(not written):
                    write_tile_breakdown(len(rows), len(cols), czi_file_path)

                for x in rows:
                    if w - x < TILE_SIZE: #Get cap on width
                        x_end = w
                    else:
                        x_end = x + TILE_SIZE

                    for y in cols:
                        # create empty tile, useful for padding out incomplete tiles at the edges with zeros
                        tile = np.zeros((TILE_SIZE, TILE_SIZE))

                        if h - y < TILE_SIZE: #Get cap on height
                            y_end = h
                        else:
                            y_end = y + TILE_SIZE

                        if DEBUG: print(f"DEBUG: x:{x}, y:{y}, x_end:{x_end}, y_end:{y_end}")
                        tile[0:x_end - x, 0:y_end - y] = im[0, x:x_end, y:y_end]

                        savedir = Path(dir_to_create, 'fov' + str(fov))

                        if(os.path.isdir(savedir)):
                            pass
                            # if DEBUG: print_colored("yellow", f"NOTE: Tried to create {savedir}, but directory {savedir} is already made!")
                        else:
                            if DEBUG: print(f"Created directory {savedir}")
                            Path.mkdir(savedir)

                        #use tifffile to write the tiff file
                        if DEBUG: print("DEBUG: Trying to write  Tiff file to " + str(savedir) + "/" + channels_to_use[channel] + ".tiff"+ " channel")

                        tf.imwrite(str(Path(savedir, channels_to_use[channel] + ".tiff")), tile)

                        if channel == 0:
                            f.write(",".join(["fov" + str(fov), str(x), str(x_end), str(y), str(y_end)]) + "\n")
                        fov += 1
        print_colored("green", f"Created {dir_to_create} with channel tiffs in fov directories!")
if __name__ == "__main__":
    tile_czi_file()

