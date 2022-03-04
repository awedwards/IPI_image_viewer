import napari
from napari_properties_plotter import PropertyPlotter as propplot
from magicgui import magicgui
import pathlib
import numpy as np
import pandas as pd
from random import randrange
import aicspylibczi
from datetime import datetime as dt
from contextlib import suppress

import time
    
LUTs = ['blue','cyan','gray','green','magenta','red','yellow']
# harcoded channel names because the CZIs only have the marke names
channel_names = ['DAPI','HLADR','CD8','CD163','CD4','XCR1','CD3','PDL1','PanCK']

# initialize dictionary to store thresholds for each channel
threshold_dict = {}

for c in channel_names:
    threshold_dict[c] = 0.0

@magicgui(
        # call_button connects to save method
        call_button = "Save thresholds",
        # name, value pairs for the marker drop down menu. The Tumor marker value is changed based on which type of sample is loaded in
        marker={"choices": [('DAPI', 'DAPI'),('CD3','CD3'),('CD4','CD4'),('CD8','CD8'),('CD163','CD163'),('XCR1','XCR1'),('HLADR','HLADR'),('PDL1','PDL1'),('Tumor','PanCK')]},
        # cell type dropdown list
        cell_type={"choices":['double_neg_t_cell',
              'cd4_t_cell',
              'cd8_t_cell',
              'mac',
              'cdc1',
              'other_myeloid_and_b_cells',
              'double_pos_t_cell']},
        threshold_slider={"widget_type": "FloatSlider",'max':20}
        )

def threshold_widget(
        # Initalize widget values

        threshold_value:float,
        threshold_slider=0.0,
        marker='DAPI',
        cell_type='cd4_t_cell',
        image_filename = pathlib.Path('/some/path.czi'),
        cell_data_filename = pathlib.Path('/some/path.csv')
        ):

    pass


@threshold_widget.image_filename.changed.connect
def load_new_image(value: str):

    im_data = []

    czifile = aicspylibczi.CziFile(value)
    dims = czifile.get_dims_shape()
    nchannels = dims[0]['C'][1]
    contrast_limits = []
    
    for c in np.arange(nchannels):
        
        #clear old channel image
        try:
            viewer.layers.pop(channel_names[c])
        except KeyError:
            pass
    
    # Check the type of sample and adjust the tumor marker value. Will need to add MUC-1 and MelanA soon.
    if "CRC" in threshold_widget.image_filename.value.name:
        channel_names[-1] = 'EPCAM'
        threshold_widget.marker.set_choice('Tumor','EPCAM')
        threshold_dict.pop('PanCK', None)
        threshold_dict['EPCAM'] = 0.0
    else:
        channel_names[-1] = 'PanCK'
        threshold_widget.marker.set_choice('Tumor','PanCK')

    for c in range(len(channel_names)):
        #load in new channel
        print("Loading channel " + str(c))
        image = czifile.read_mosaic(C=c,scale_factor=1)
        contrast_limits.append([0, 2**16])

        # add each channel
        viewer.add_image(image, name=channel_names[c], visible=False, contrast_limits=contrast_limits[c])
        viewer.layers[c].colormap = LUTs[randrange(len(LUTs))]
        viewer.layers[c].opacity = 1.0
        viewer.layers[c].blending = 'additive'
        viewer.layers[c].interpolation = 'gaussian'

@threshold_widget.cell_data_filename.changed.connect
def load_cell_data(value: str):
    
    #clear old points data
    try:
        viewer.layers['points'].pop()
        viewer.layers['cell type results'].pop()
    except KeyError:
        pass

    data = pd.read_csv(value)
    n_points = len(data)
    x = np.array(data['centroid-0'])
    y = np.array(data['centroid-1'])

    # Initalized column for threshold results
    for c in channel_names:
        data[c + "_expressed"] = np.array(data[c] > threshold_dict[c], dtype=np.int8)
    
    # format data for adding as layer
    points = np.stack((x, y)).transpose()
    points_layer = viewer.add_points(points,
            size=25,
            properties=data,
            face_color=threshold_widget.marker.value,
            name='points',
            visible=True,
            shown=[True]*len(data))
    if "cell_type" in data.columns:
        points_layer = viewer.add_points(points,
            size=25,
            properties=data,
            name='cell type results',
            visible=True,
            shown=data['cell_type']==threshold_widget.cell_type.value)
    else:
        points_layer = viewer.add_points(points,
            size=25,
            properties=data,
            name='cell type results',
            visible=True,
            shown=[True]*len(data))

@threshold_widget.threshold_slider.changed.connect
def threshold_slider_change(value: float):
    
    threshold_widget.threshold_value.value = value
    
    channel = threshold_widget.marker.value
    data = pd.DataFrame.from_dict(viewer.layers['points'].properties)
    
    x = np.array(data['centroid-0'])
    y = np.array(data['centroid-1'])
    all_points = np.stack((x,y)).transpose()

    thresholded_idx = data[channel] > value
    
    threshold_dict[channel] = value
    
    data[channel + "_expressed"] = np.array(thresholded_idx, dtype=np.uint8)
    viewer.layers['points'].shown = thresholded_idx
    viewer.layers['points'].properties = data

    update_cell_types()
    cell_type_changed(threshold_widget.cell_type.value)

@threshold_widget.cell_type.changed.connect
def cell_type_changed(value: str):
    
    data = pd.DataFrame.from_dict(viewer.layers['points'].properties)
    
    if 'cell_type' not in data.columns:
        return
    
    ct_idx = data['cell_type']==value

    viewer.layers['cell type results'].shown = ct_idx

@threshold_widget.threshold_value.changed.connect
def threshold_value_changed(value: float):
    threshold_widget.threshold_slider.value = value

@threshold_widget.call_button.clicked.connect
def save():

    filepath = threshold_widget.cell_data_filename.value.parent
    name = threshold_widget.image_filename.value.name
    
    with open(pathlib.Path(filepath, name[:-4] + "_thresholds.txt"),"a") as f:
        for c in channel_names:
            line = dt.now().strftime('%Y-%m-%d %H:%M') + ": "+c+" gate set to "+str(threshold_dict[c])+"\n"
            f.write(line)

    print("Thresholds saved at " + str(pathlib.Path(filepath, name[:-4] + "_thresholds.txt")))

    data = pd.DataFrame.from_dict(viewer.layers['points'].properties)
    for channel in threshold_dict.keys():
        data[channel + "_expressed"] = np.array(data[channel] > threshold_dict[channel], dtype=np.int8)
    
    data.to_csv(pathlib.Path(filepath, name[:-4] + "_single_cell_data_gated" + dt.now().strftime('%Y%m%d') + ".csv"))

def update_cell_types():
   
    data = pd.DataFrame.from_dict(viewer.layers['points'].properties)
    
    cell_types = ['double_neg_t_cell',
              'cd4_t_cell',
              'cd8_t_cell',
              'mac',
              'cdc1',
              'other_myeloid_and_b_cells',
              'double_pos_t_cell']
    
    data['cell_type'] = ['other']*len(data)

    ct_idx = np.zeros((len(data), len(cell_types)),dtype=bool)
    
    ct_idx[:,0]= (data['DAPI_expressed']==1) & \
                (data['CD3_expressed']==1) & \
                (data['CD4_expressed']==0) & \
                (data['CD8_expressed']==0) & \
                (data['XCR1_expressed']==0)

    ct_idx[:,1]=(data['DAPI_expressed']==1) & \
                (data['CD4_expressed']==1) & \
                (data['CD3_expressed']==1) & \
                (data['CD8_expressed']==0) & \
                (data['XCR1_expressed']==0)

    ct_idx[:,2]=(data['DAPI_expressed']==1) & \
                (data['CD8_expressed']==1) & \
                (data['CD3_expressed']==1) & \
                (data['CD4_expressed']==0) & \
                (data['XCR1_expressed']==0)

    ct_idx[:,3]=(data['DAPI_expressed']==1) & \
                (data['CD163_expressed']==1) & \
                (data['HLADR_expressed']==1) & \
                (data['XCR1_expressed']==0) & \
                (data['CD3_expressed']==0)

    ct_idx[:,4]=(data['DAPI_expressed']==1) & \
                (data['XCR1_expressed']==1) & \
                (data['HLADR_expressed']==1) & \
                (data['CD3_expressed']==0) & \
                (data['CD163_expressed']==0)

    ct_idx[:,5]=(data['DAPI_expressed']==1) & \
                (data['HLADR_expressed']==1) & \
                (data['CD163_expressed']==0) & \
                (data['CD3_expressed']==0) & \
                (data['XCR1_expressed']==0)

    ct_idx[:,6]=(data['DAPI_expressed']==1) & \
                (data['CD3_expressed']==1) & \
                (data['CD4_expressed']==1) & \
                (data['CD8_expressed']==1)& \
                (data['XCR1_expressed']==0)

    assigned_twice = np.sum(ct_idx,axis=1)>1
    
    for i,ct in enumerate(cell_types):
        
        data.loc[ct_idx[:,i],'cell_type']=ct
        
    data.loc[assigned_twice,'cell_type'] = 'assigned_twice'

    viewer.layers['points'].properties = data  

viewer = napari.Viewer()

viewer.window.add_dock_widget(threshold_widget)
with suppress(TypeError):
	pp = propplot(viewer)
	viewer.window.add_dock_widget(pp, area='bottom')
napari.run()
