# coding: utf-8

#################################################################################################################
#Filename: final_viewer_main.py

# GOAL: View CZI image with threshold points and segmentation boundaries

# INPUT: OPTIONAL input of czi image file (-i), cell_data points (-p), cell_boundaries file (-b)

# OUTPUT: None though OPTIONAL output of cell threshold and data information

# Parse arguments for napari viewer.
#
# optional arguments:
#   -h, --help            show this help message and exit
#   --debug, -d           Print additional information to the terminal when running script
#   --points POINTS, -p POINTS
#                         Specify a single file with points label to preload in
#   --image IMAGE, -i IMAGE
#                         Specify a single czi image file to preload in
#   --bounds BOUNDS, -b BOUNDS
#                         Specify a single file with segementation boundaries to load in

#################################################################################################################
#TODO: Fix update_celL_types
#TODO: Multithread the channels #NOTE: Version 1 of Multithreading is added (but there may be more to save)
#TODO: Test and Make Usable on windows OS (support)
#TODO: Threshold file is overwritten when opening or refreshing?
#TODO: Add DEBUG information
#TODO: Fix Same File Load Bug
#TODO: Color in Boundaries
#   Take bound and labels --> and change color based on if __
#TODO: Make large tiff files viewable, and integrate with imagej/fiji

#COMPLETED: Add channels from command line
#COMPLETED: Make slider threshold value be able to change range
#COMPLETED: Make boundaries get from tile_metadata.txt.
#COMPLETED: Make clear files optional
#COMPLETED: V1 of Multithreading image input

from PIL import Image
from aicsimageio import AICSImage
from contextlib import suppress
from datetime import datetime as dt
from magicgui import magicgui
from napari_properties_plotter import PropertyPlotter as propplot
from random import randrange

import aicspylibczi
import csv
import napari
import numpy as np
import os
import pandas as pd
import pathlib
import time
import traceback
import warnings
from segmentation_utils import print_colored, napari_viewer_parser, restricted_float
from napari.qt.threading import thread_worker

warnings.simplefilter(action='ignore', category=FutureWarning)

LUTs = ['blue', 'cyan', 'gray', 'green', 'magenta', 'red', 'yellow']
CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))

all_channels = [
    ['DAPI','HLADR','CD8','CD163','CD4','XCR1','CD3','PDL1','EPCAM'],
    ['DAPI','HLADR','CD8','CD163','CD4','XCR1','CD3','PDL1','PanCK']
]
DEFAULT_CHANNELS_TO_USE = 1
DEFAULT_SPINBOX_STEP = .5
DEBUG = False

channel_names = all_channels[DEFAULT_CHANNELS_TO_USE]
current_spinbox_step = DEFAULT_SPINBOX_STEP
cell_type_col = False

#Using argument parser to organize the input
napari_viewer_parser.add_argument("--channel", "-c", dest='channel', action="store_true",
                    help="View, Select, and Add Channels!")

napari_viewer_parser.add_argument("--step", "-s", dest='step', type=restricted_float,
                    help=f"This sets the step (or how much the value can be decreased or increased)"
                         f"for the threshold_value number. (Default is {DEFAULT_SPINBOX_STEP}")

napari_viewer_parser_args = napari_viewer_parser.parse_args()
if napari_viewer_parser_args.debug:
    DEBUG = True

if napari_viewer_parser_args.step:
    current_spinbox_step = napari_viewer_parser_args.step

#Gets proper channels from above for tiling
def get_channel_choice(len_of_channels):
    while True:
        try:
            number = int(input('Chose an option from menu: '))
            if 0 <= number <= len_of_channels:
                return number
            else:
                raise ValueError("Not in valid number range!")
        except Exception as  e:
            print_colored("red", f"Invalid channel choice! {e}")

if napari_viewer_parser_args.channel:
    print_colored("cyan", f"Select which channels you would like to use. Optionally, use your own. (Note: Current Default is option {DEFAULT_CHANNELS_TO_USE})")

    print_colored("green", f"CHANNEL CHOICES")
    for index, channel in enumerate(all_channels):
        print_colored("green", f"[{index}] {channel}")
    print_colored("green", f"[{len(all_channels)}] Input Custom Channels")
    choices = list(range(0, len(all_channels) + 1))

    num_to_use = get_channel_choice(len(all_channels))

    if(num_to_use < len(all_channels)):
        channel_names = all_channels[num_to_use]
    else:
        print_colored("cyan", f"Type each channel out separated by a comma, and custom channels will be created. (No trailing comma!)")
        channel_names = input("Give Channels: ").replace(" ", "").split(',')

def calculate_final_dimensions_from_metadata(tile_metadata_file_path):
    with open(tile_metadata_file_path, newline='') as metadata_file_path:
        tile_position_data = csv.DictReader(metadata_file_path)
        final_dimensions = metadata_file_path.readlines()[-1].split(',')
        rows = len(np.arange(0, int(final_dimensions[2]), 2048))  #final x dimension
        cols = len(np.arange(0, int(final_dimensions[4]), 2048))  #final y dimensions
        return rows, cols

# initialize dictionary to store thresholds for each channel
threshold_dict = {}

for c in channel_names:
    threshold_dict[c] = 0.0

@magicgui(
    # call_button connects to save method
    threshold_value={"widget_type": "FloatSpinBox", 'value': 0.0, 'step': current_spinbox_step, 'max': 20.0},
    call_button="Save thresholds",
    # name, value pairs for the marker drop down menu. The Tumor marker value is changed based on which type of sample is loaded in
    marker={"choices": [('DAPI', 'DAPI'), ('CD3', 'CD3'), ('CD4', 'CD4'), ('CD8', 'CD8'), ('CD163', 'CD163'),
                        ('XCR1', 'XCR1'), ('HLADR', 'HLADR'), ('PDL1', 'PDL1'), ('Tumor', 'PanCK')]},
    # cell type dropdown list
    cell_type={"choices": ['double_neg_t_cell',
                           'cd4_t_cell',
                           'cd8_t_cell',
                           'mac',
                           'cdc1',
                           'other_myeloid_and_b_cells',
                           'double_pos_t_cell']},
    threshold_slider={"widget_type": "FloatSlider", 'max': 20},
    clear_layers_button={"label": "Clear Layers with New Image",
                     "widget_type": "RadioButtons", 'choices': ["Yes", "No"]}

)
def threshold_widget(
        threshold_value: float,
        threshold_slider=0.0,
        marker='DAPI',
        cell_type='cd4_t_cell',
        czi_image_filename=pathlib.Path("<Select File>"),
        cell_data_filename=pathlib.Path("<Select File>"),
        cell_boundaries_filename=pathlib.Path("<Select File>"),
        clear_layers_button='No'
): pass

def update_layer(inputs):
    new_image = inputs[0]
    channel_name = inputs[1]
    viewer.add_image(new_image, name=channel_name, visible=False, contrast_limits=[0, 2 ** 16],
                     colormap=LUTs[randrange(len(LUTs))], opacity=1.0, blending='additive',
                     interpolation='gaussian')

@thread_worker(connect={"returned": update_layer})
def multiload_image(czi_file_path, ind, chan):
    czi_file = aicspylibczi.CziFile(czi_file_path) #This is the part that is not memory efficient wise, tho it is faster
    return czi_file.read_mosaic(C=ind, scale_factor=1), chan

@threshold_widget.czi_image_filename.changed.connect
def load_new_image(value: str):
    #Clears out old channel values when a new image is loaded using widget gui
    if(threshold_widget.clear_layers_button.value == 'Yes'):
        print_colored("yellow", f"Clearing Out Old Layers")
        while(len(viewer.layers) != 0):
            viewer.layers.pop()
    start_time = time.perf_counter()
    #Add each channel img to layers with associated name
    for index, channel_name in enumerate(channel_names):
        print_colored("cyan", f"Loading channel {index} - {channel_name}")
        #This line starts a new worker thread which reads in an image, and safely
        #adds it to the viewer
        multiload_image(value, index, channel_name) #
    finish_time = time.perf_counter()
    if DEBUG: print_colored("green", f"Final Time: {finish_time - start_time}")
    threshold_widget.marker.set_choice('Tumor', channel_names[-1])

@threshold_widget.cell_data_filename.changed.connect
def load_cell_data(cell_data_file_path: str):

    # clear old points data
    if ('points' in viewer.layers):
        del viewer.layers['points']

    if ('cell type results' in viewer.layers):
        del viewer.layers['cell type results']

    data = pd.read_csv(cell_data_file_path) #(15960, 56)
    x = np.array(data['centroid-0']) # (15960,)
    y = np.array(data['centroid-1']) # (15960,)

    # Initalized column for threshold results (data --> (15960, 64))
    for c in channel_names:
        data[c + "_expressed"] = np.array(data[c] > threshold_dict[c], dtype=np.int8)
    # format data for adding as layer
    points = np.stack((x, y)).transpose() #(15960, 2))

    points_layer = viewer.add_points(points,
        size=25,
        properties=data,
        face_color=threshold_widget.marker.value,
        name='points',
        visible=True,
        shown=np.array([True] * len(data)),
    )

    if("cell_type" in data.columns):
        cell_type_col = True
        shown_data = data['cell_type'] == threshold_widget.cell_type.value
    else:
        shown_data = np.array([True]) * len(data)

    cell_type_layer = viewer.add_points(points,
        size=25,
        properties=data,
        name='cell type results',
        visible=True,
        shown=shown_data)

@threshold_widget.threshold_slider.changed.connect
def threshold_slider_change(value: float):
    threshold_widget.threshold_value.value = value
    channel = threshold_widget.marker.value
    threshold_dict[channel] = value

    data = pd.DataFrame.from_dict(viewer.layers['points'].properties)
    thresholded_idx = data[channel] > value

    data[channel + "_expressed"] = np.array(thresholded_idx, dtype=np.uint8)
    viewer.layers['points'].shown = thresholded_idx
    viewer.layers['points'].properties = data

    update_cell_types()
    cell_type_changed(threshold_widget.cell_type.value)

@threshold_widget.threshold_value.changed.connect
def threshold_value_changed(value: float):
    threshold_widget.threshold_slider.value = value

@threshold_widget.marker.changed.connect
def marker_changed(value: str):
    threshold_widget.threshold_value.value = threshold_dict[value]
    threshold_widget.threshold_slider.value = threshold_dict[value]

@threshold_widget.call_button.clicked.connect
def save():
    filepath = threshold_widget.cell_data_filename.value.parent
    name = threshold_widget.czi_image_filename.value.name

    with open(pathlib.Path(filepath, name[:-4] + "_thresholds.txt"), "a") as f:
        for c in channel_names:
            line = dt.now().strftime('%Y-%m-%d %H:%M') + ": " + c + " gate set to " + str(threshold_dict[c]) + "\n"
            f.write(line)

    print("Thresholds saved at " + str(pathlib.Path(filepath, name[:-4] + "_thresholds.txt")))

    data = pd.DataFrame.from_dict(viewer.layers['points'].properties)
    for channel in threshold_dict.keys():
        data[channel + "_expressed"] = np.array(data[channel] > threshold_dict[channel], dtype=np.int8)

    data.to_csv(pathlib.Path(filepath, name[:-4] + "_single_cell_data_gated" + dt.now().strftime('%Y%m%d') + ".csv"))

@threshold_widget.cell_type.changed.connect
def cell_type_changed(value: str):
    if(cell_type_col):
        data = pd.DataFrame.from_dict(viewer.layers['points'].properties)
        ct_idx = data['cell_type'] == value
        viewer.layers['cell type results'].shown = ct_idx

def update_cell_types():
    data = pd.DataFrame.from_dict(viewer.layers['points'].properties)

    cell_types = ['double_neg_t_cell',
                  'cd4_t_cell',
                  'cd8_t_cell',
                  'mac',
                  'cdc1',
                  'other_myeloid_and_b_cells',
                  'double_pos_t_cell',
                  'double_pos_xrc1_cd163']

    data['cell_type'] = ['other'] * len(data)

    ct_idx = np.zeros((len(data), len(cell_types)), dtype=bool)

    ct_idx[:, 0] = (data['DAPI_expressed'] == 1) & \
                   (data['CD3_expressed'] == 1) & \
                   (data['CD4_expressed'] == 0) & \
                   (data['CD8_expressed'] == 0) & \
                   (data['XCR1_expressed'] == 0)

    ct_idx[:, 1] = (data['DAPI_expressed'] == 1) & \
                   (data['CD4_expressed'] == 1) & \
                   (data['CD3_expressed'] == 1) & \
                   (data['CD8_expressed'] == 0) & \
                   (data['XCR1_expressed'] == 0)

    ct_idx[:, 2] = (data['DAPI_expressed'] == 1) & \
                   (data['CD8_expressed'] == 1) & \
                   (data['CD3_expressed'] == 1) & \
                   (data['CD4_expressed'] == 0) & \
                   (data['XCR1_expressed'] == 0)

    ct_idx[:, 3] = (data['DAPI_expressed'] == 1) & \
                   (data['CD163_expressed'] == 1) & \
                   (data['HLADR_expressed'] == 1) & \
                   (data['XCR1_expressed'] == 0) & \
                   (data['CD3_expressed'] == 0)

    ct_idx[:, 4] = (data['DAPI_expressed'] == 1) & \
                   (data['XCR1_expressed'] == 1) & \
                   (data['HLADR_expressed'] == 1) & \
                   (data['CD3_expressed'] == 0) & \
                   (data['CD163_expressed'] == 0)

    ct_idx[:, 5] = (data['DAPI_expressed'] == 1) & \
                   (data['HLADR_expressed'] == 1) & \
                   (data['CD163_expressed'] == 0) & \
                   (data['CD3_expressed'] == 0) & \
                   (data['XCR1_expressed'] == 0)

    ct_idx[:, 6] = (data['DAPI_expressed'] == 1) & \
                   (data['CD3_expressed'] == 1) & \
                   (data['CD4_expressed'] == 1) & \
                   (data['CD8_expressed'] == 1) & \
                   (data['XCR1_expressed'] == 0)
    
    ct_idx[:, 7] = (data['DAPI_expressed'] == 1) & \
                   (data['CD3_expressed'] == 0) & \
                   (data['XCR1_expressed'] == 1) & \
                   (data['CD163_expressed'] == 1)
    
    assigned_twice = np.sum(ct_idx, axis=1) > 1

    for i, ct in enumerate(cell_types):
        data.loc[ct_idx[:, i], 'cell_type'] = ct

    data.loc[assigned_twice, 'cell_type'] = 'assigned_twice'

    viewer.layers['points'].properties = data

@threshold_widget.cell_boundaries_filename.changed.connect
def get_boundaries(boundaries_file_path: str):
    segmented_cell_borders_filename = boundaries_file_path.stem.split('-')[-1]
    tile_metadata_path = pathlib.Path(CURRENT_DIR, 'final_data',
                                      f"{segmented_cell_borders_filename}_dir", 'tile_metadata.txt')

    try:
        rows, cols = calculate_final_dimensions_from_metadata(tile_metadata_path)
    except Exception as e:
        print_colored("red", f"Could not open or load tile_metadata_file_path."
                             f"Ensure associated meta data is in the final_data directory!")
        print(f"{e} \n{traceback.format_exc()}")

    try:
        with open(boundaries_file_path, 'rb') as boundaries_file:
            a = np.load(boundaries_file).reshape(rows*cols, 2048, 2048)
    except Exception as e:
        print_colored("red", f"Could not open or load {boundaries_file_path}")
        print(f"{e} \n{traceback.format_exc()}")
        return

    #Need to get czi dimensions to properly stitch boundaries together

    final_concat = None
    for i in range(rows):
        row_concat = a[cols * i].reshape(2048, 2048)
        row_concat = np.where(row_concat > 0, 1, 0)  # fov + (3135 * i*3+j), 0)

        for j in range(1, cols):
            fov = a[i*cols+j].reshape(2048, 2048)
            fov = np.where(fov > 0, 1, 0) #ensure it is all 0 (Black) and 1(White)
            row_concat = np.concatenate([row_concat, fov], axis=1)

        if(i == 0):
            final_concat = row_concat
        else:
            final_concat = np.concatenate([final_concat, row_concat], axis=0)

    #Here we are taking the image, converting it to RBA so we can make deadspace transparent
    #Then adding it a an image to the layers

    temp_img = Image.fromarray((final_concat * 255).astype(np.uint8)) #Need to alter to 255 scale as it is grayscale
    temp_img_with_transparency = temp_img.convert("RGBA") #A is Alpha for transparency
    array_with_transparency = np.asarray(temp_img_with_transparency)
    array_with_transparency[:, :, 3] = final_concat * 255 #Adjust transparency #White
    # array_with_transparency[:, :, 2] = final_concat * 0 #Adjust transparency #Yellow
    # array_with_transparency[:, :, 1] = final_concat * 0 #Adjust transparency #Magenta
    # array_with_transparency[:, :, 0] = final_concat * 0 #Adjust transparency #Cyan
    viewer.add_image(array_with_transparency, name="NPY Bounds")


viewer = napari.Viewer()
viewer.window.add_dock_widget(threshold_widget)

pp = propplot(viewer)
viewer.window.add_dock_widget(pp, area='bottom')

if(napari_viewer_parser_args.image):
    threshold_widget.czi_image_filename.value = napari_viewer_parser_args.image[0]

if(napari_viewer_parser_args.points):
    threshold_widget.cell_data_filename.value = napari_viewer_parser_args.points[0]

if(napari_viewer_parser_args.bounds):
    threshold_widget.cell_boundaries_filename.value = napari_viewer_parser_args.bounds[0]

napari.run()

