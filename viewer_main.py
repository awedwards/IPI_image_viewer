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
#TODO: Add DEBUG information
#TODO: Fix Same File Load Bug

from PIL import Image
from contextlib import suppress
from datetime import datetime as dt
from magicgui import magicgui
from napari_properties_plotter import PropertyPlotter as propplot
from random import randrange

import aicspylibczi
import json
import napari
import numpy as np
import os
import pandas as pd
import pathlib
import traceback
import warnings
from segmentation_utils import print_colored, napari_viewer_parser

warnings.simplefilter(action='ignore', category=FutureWarning)

LUTs = ['blue', 'cyan', 'gray', 'green', 'magenta', 'red', 'yellow']
CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))

all_channels = [
    ['DAPI','HLADR','CD8','CD163','CD4','XCR1','CD3','PDL1','EPCAM'],
    ['DAPI','HLADR','CD8','CD163','CD4','XCR1','CD3','PDL1','PanCK']
]
DEFAULT_CHANNELS_TO_USE = 1

channel_names = all_channels[DEFAULT_CHANNELS_TO_USE]

cell_type_col = False

#Using argument parser to organize the input
napari_viewer_parser.add_argument("--channel", "-c", dest='channel', action="store_true",
                    help="View, Select, and Add Channels!")

napari_viewer_parser_args = napari_viewer_parser.parse_args()

if napari_viewer_parser_args.debug:
    DEBUG = True

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

# initialize dictionary to store thresholds for each channel
threshold_dict = {}

for c in channel_names:
    threshold_dict[c] = 0.0

@magicgui(
    # call_button connects to save method
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
    threshold_slider={"widget_type": "FloatSlider", 'max': 20}
)
def threshold_widget(
        threshold_value: float,
        threshold_slider=0.0,
        marker='DAPI',
        cell_type='cd4_t_cell',
        czi_image_filename = pathlib.Path("<Select File>"),
        cell_data_filename = pathlib.Path("<Select File>"),
        cell_boundaries_filename=pathlib.Path("<Select File>")
): pass
from aicsimageio import AICSImage

@threshold_widget.czi_image_filename.changed.connect
def load_new_image(value: str):
    #TODO: Make large tiff files viewable, and integrate with imagej/fiji
    contrast_limits = []

    czi_file = aicspylibczi.CziFile(value)

    #TODO: Make optional
    #Clears out old channel values when a new image is loaded using widget gui
    while(len(viewer.layers) != 0):
        viewer.layers.pop()
    #Add each channel img to layers with associated name
    for index, channel_name in enumerate(channel_names):
        print_colored("cyan", f"Loading channel {index} - {channel_name}")
        image = czi_file.read_mosaic(C=index, scale_factor=1)
        contrast_limits.append([0, 2**16])
        # add each channel
        viewer.add_image(image, name=channel_name, visible=False, contrast_limits=contrast_limits[index])
        viewer.layers[channel_name].colormap = LUTs[randrange(len(LUTs))]
        viewer.layers[channel_name].opacity = 1.0
        viewer.layers[channel_name].blending = 'additive'
        viewer.layers[channel_name].interpolation = 'gaussian'

    threshold_widget.marker.set_choice('Tumor',channel_names[-1])

@threshold_widget.cell_data_filename.changed.connect
def load_cell_data(cell_data_file_path: str):
    #TODO: Keep track of number points and output to napari

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
                  'double_pos_t_cell']

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

    assigned_twice = np.sum(ct_idx, axis=1) > 1

    for i, ct in enumerate(cell_types):
        data.loc[ct_idx[:, i], 'cell_type'] = ct

    data.loc[assigned_twice, 'cell_type'] = 'assigned_twice'

    viewer.layers['points'].properties = data

@threshold_widget.cell_boundaries_filename.changed.connect
def get_boundaries(boundaries_file_path: str):
    segmented_cell_borders_filename = boundaries_file_path.stem.split('-')[-1]

    try:
        with open(f'{CURRENT_DIR}/final_data/.{segmented_cell_borders_filename}_metadata.json', 'rb') as stitching_dims_file:
            data = json.load(stitching_dims_file)
    except Exception as e:
        print_colored("red", f"Could not open or load {CURRENT_DIR}/final_data/.{segmented_cell_borders_filename}.czi'."
                             f"Ensure that the czi file was tiled the associated meta data is in the final_data directory!")
        print(e)
        print(traceback.format_exc())
        return

    rows = data['dims']['rows']
    cols = data['dims']['cols']

    try:
        with open(boundaries_file_path, 'rb') as boundaries_file:
            a = np.load(boundaries_file).reshape(rows*cols, 2048, 2048)
    except Exception as e:
        print_colored("red", f"Could not open or load {boundaries_file_path}")
        print(e)
        print(traceback.format_exc())
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
    array_with_transparency[:, :, 3] = final_concat * 255 #Adjust transparency
    viewer.add_image(array_with_transparency, name="NPY Bounds")


viewer = None
viewer = napari.Viewer()
viewer.window.add_dock_widget(threshold_widget)

pp = propplot(viewer)
viewer.window.add_dock_widget(pp, area='bottom')

if(napari_viewer_parser_args.image):
    load_new_image(napari_viewer_parser_args.image[0])

if(napari_viewer_parser_args.points):
    load_cell_data(napari_viewer_parser_args.points[0])

if(napari_viewer_parser_args.bounds):
    get_boundaries(napari_viewer_parser_args.bounds[0])

print_colored("cyan", "Running napari...")

napari.run()
