import napari
from napari_properties_plotter import PropertyPlotter as propplot
from magicgui import magicgui
import pathlib
import numpy as np
import pandas as pd
from random import randrange
import aicspylibczi
from datetime import datetime as dt
    
LUTs = ['blue','cyan','gray','green','magenta','red','yellow']
channel_names = ['DAPI','HLADR','CD8','CD163','CD4','XCR1','CD3','PDL1','PanCK']

threshold_dict = {}

for c in channel_names:
    threshold_dict[c] = 0.0

@magicgui(
        call_button = "Save thresholds",
        marker={"choices": [('DAPI', 'DAPI'),('CD3','CD3'),('CD4','CD4'),('CD8','CD8'),('CD163','CD163'),('XCR1','XCR1'),('HLADR','HLADR'),('PDL1','PDL1'),('Tumor','PanCK')]},
        threshold_slider={"widget_type": "FloatSlider",'max':20}
        )

def threshold_widget(
        threshold_value:float,
        threshold_slider=0.0,
        marker='DAPI',
        image_filename = pathlib.Path('/some/path.czi'),
        cell_data_filename = pathlib.Path('/some/path.csv')
        ):

    pass


@threshold_widget.image_filename.changed.connect
def load_new_image(value: str):

    # clear old image

    im_data = []

    czifile = aicspylibczi.CziFile(value)
    dims = czifile.get_dims_shape()
    nchannels = dims[0]['C'][1]

    for c in np.arange(nchannels):
        print("Loading channel " + str(c))
        im_data.append(czifile.read_mosaic(C=c,scale_factor=1))
    
    image = np.stack(im_data)
    if "CRC" in threshold_widget.image_filename.value.name:
        channel_names[-1] = 'EPCAM'
        threshold_widget.marker.set_choice('Tumor','EPCAM')
        threshold_dict.pop('PanCK', None)
        threshold_dict['EPCAM'] = 0.0
    else:
        channel_names[-1] = 'PanCK'
        threshold_widget.marker.set_choice('Tumor','PanCK')

    for c in range(len(channel_names)):
        
        #clear old channel image
        try:
            viewer.layers.pop(channel_names[c])
        except KeyError:
            pass

        viewer.add_image(image[c], name=channel_names[c], visible=False)
        viewer.layers[c].colormap = LUTs[randrange(len(LUTs))]
        viewer.layers[c].opacity = 1.0
        viewer.layers[c].blending = 'additive'
        viewer.layers[c].interpolation = 'gaussian'

@threshold_widget.cell_data_filename.changed.connect
def load_cell_data(value: str):
    
    #clear old points data
    try:
        viewer.layers['points'].pop()
    except KeyError:
        pass

    data = pd.read_csv(value)
    n_points = len(data)
    x = np.array(data['centroid-0'])
    y = np.array(data['centroid-1'])

    points = np.stack((x, y)).transpose()
    
    points_layer = viewer.add_points(points,
            size=25,
            properties=data,
            face_color=threshold_widget.marker.value,
            name='points')

@threshold_widget.threshold_slider.changed.connect
def threshold_slider_change(value: float):
    
    try:
        viewer.layers.pop('threshold result')
    except KeyError:
        pass

    data = pd.DataFrame.from_dict(viewer.layers['points'].properties)
    data = data[data[threshold_widget.marker.value] > value]
    x = np.array(data['centroid-0'])
    y = np.array(data['centroid-1'])

    points = np.stack((x, y)).transpose()
    if len(points) > 0: 
        points_layer = viewer.add_points(points,
                size=25,
                properties=data,
                face_color=threshold_widget.marker.value,
                name='threshold result')
    threshold_widget.threshold_value.value = value
    threshold_dict[threshold_widget.marker.value] = value

@threshold_widget.marker.changed.connect
def update_points_facecolor(value: str):
    viewer.layers['points'].face_color = value
    threshold_widget.threshold_slider.value = threshold_dict[value]
    threshold_widget.threshold_value.value = threshold_dict[value]

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
    
    for c in channel_names:

        data[c + "_expressed"] = np.array(data[c] > threshold_dict[c], dtype=np.int8)
    
    data.to_csv(pathlib.Path(filepath, name[:-4] + "_single_cell_data_gated" + dt.now().strftime('%Y%m%d') + ".csv"))

viewer = napari.Viewer()

viewer.window.add_dock_widget(threshold_widget)
pp = propplot(viewer)
viewer.window.add_dock_widget(pp, area='bottom')
napari.run()
