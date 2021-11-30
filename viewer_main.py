import napari
from magicgui import magicgui
import pathlib
import numpy as np
import pandas as pd
from random import randrange
import aicspylibczi


LUTs = ['blue','cyan','gray','green','magenta','red','yellow']
channel_names = ['DAPI','HLADR','CD8','CD163','CD4','XCR1','CD3','PDL1','EPCAM']

@magicgui(
        call_button = "Save thresholds",
        marker={"choices": ['DAPI','CD8','CD4','CD163','XCR1','HLADR','CD3','PDL1','Tumor']},
        Threshold={"widget_type": "FloatSlider",'max':20}
        )

def threshold_widget(
        Threshold_value:float,
        Threshold=0.0,
        marker='DAPI',
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

    for c in np.arange(nchannels):
        print("Loading channel " + str(c))
        im_data.append(czifile.read_mosaic(C=c,scale_factor=1))
    
    image = np.stack(im_data)

    for c in range(len(channel_names)):

        viewer.add_image(image[c], name=channel_names[c], visible=False)
        viewer.layers[c].colormap = LUTs[randrange(len(LUTs))]
        viewer.layers[c].opacity = 1.0
        viewer.layers[c].blending = 'additive'
        viewer.layers[c].interpolation = 'gaussian'
   
    threshold_dict = {}
    for c in channel_names:
        threshold_dict[c] = 0.0

@threshold_widget.cell_data_filename.changed.connect
def load_cell_data(value: str):

    data = pd.read_csv(value)

    x = np.array(data['centroid-0'])
    y = np.array(data['centroid-1'])

    points = np.stack((x, y)).transpose()

    points_layer = viewer.add_points(points, size=25)

@threshold_widget.Threshold.changed.connect
def threshold_slider_change(value: float):

    threshold_dict[threshold_widget.marker] = value
    
    result = data[data[threshold_widget.marker] > value]

    print(len(result))



viewer = napari.Viewer()
viewer.window.add_dock_widget(threshold_widget)
   
napari.run()

