from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QWidget, QFileDialog, QSlider
from magicgui import magicgui
from magicgui.widgets import ComboBox

from random import randrange
import sys
import imageio
import napari
import exifread
import pathlib
import numpy as np
import pandas as pd
import aicspylibczi

class App(QWidget):

    #@magicgui(
    #        auto_call=True, threshold={"widget_type": QSlider, "minimum":1, "maximum":10}
    #        )

    def __init__(self):
        super().__init__()
        
        # Limits to the linear monochrome LUTs within napari
        self.available_LUTs = ['blue', 'cyan', 'gray', 'green', 'magenta', 'red', 'yellow']

        self.run_all()

    #def adjust_threshold(channel: str, threshold: int=5):
    #    print(threshold)

    # Prompt to select MIBItiff. No error handling.
    def get_working_directory(self, message):
        path_str = (QFileDialog.getOpenFileName(self, message, "/media/austin/DrosophilaMelanogaster/IPI/8plex/data"))
        return pathlib.Path(str(path_str[0]))

    def get_saving_directory(self, message):
        path_str = (QFileDialog.getSaveFileName(self, message))
        return pathlib.Path(str(path_str[0]))

    # Opens the tiff file
    def open_tiff(self, file_name):
        image = imageio.mimread(file_name, multifile=True)
        return image

    # Gets the channels from the tags
    def get_tags(self, file_name, num_images):
        f = open(file_name, 'rb')
        tags = exifread.process_file(f)
        channel_list = [str(tags['Image PageName']), str(tags['Thumbnail PageName'])]
        for n in range(2, num_images):
            channel_list.append(str(tags['IFD ' + str(n) + ' PageName']))

        return channel_list

    # Determine how many channels there should be
    def get_image_length(self, image):
        num_images = len(image)
        return num_images

    # launches napari with the fle and channel names
    def launch_napari(self, image, channel_names, LUTs):

        self.viewer = napari.Viewer()

        self.file_menu = self.viewer.window.file_menu
        self.add_load_csv_menu()
        self.threshold_value = 0 # initialize threshold
        for c in range(len(channel_names)):
            self.viewer.add_image(image[c], name=channel_names[c], visible=False)
            self.viewer.layers[c].name = channel_names[c]
            self.viewer.layers[c].colormap = LUTs[randrange(len(LUTs))]
            self.viewer.layers[c].opacity = 1.0
            self.viewer.layers[c].blending = 'additive'
            self.viewer.layers[c].interpolation = 'gaussian'
        
        self.create_channel_menu(channel_names)
        #self.slider = scrollThreshold.Gui()
        #self.slider.threshold_widget.channel_changed.connect(update_slider)

    #def update_slider(self,channel):
    #    
    #    self.slider.threshold_widget.setMaximum(np.max(self.cell_data['channel']))
    #    self.slider.threshold_widget.setMinimum(np.min(self.cell_data['channel']))

    def load_csv(self):
        csvfilter = "Comma-Separated Value File (*.csv)"
        path_str = QFileDialog.getOpenFileName(self, "Select data table for the current image", "/media/austin/DrosophilaMelanogaster/IPI/8plex/data", filter=csvfilter)
        self.cell_data = pd.read_csv(pathlib.Path(str(path_str[0])))

        #self.x = np.array(self.cell_data['centroid-0'])
        #self.y = np.array(self.cell_data['centroid-1'])
        #self.points = np.stack((self.x,self.y)).transpose()
        #self.points_layer = self.viewer.add_points(self.points, size = 20)

    def add_load_csv_menu(self):
        self.csv_menu = self.file_menu.addAction('Load cell data table...')
        self.csv_menu.triggered.connect(self.load_csv)
    
    def create_channel_menu(self, labels):
        self.channel_menu = ComboBox(label='channel_threshold', choices=labels)
        self.viewer.window.add_dock_widget(self.channel_menu)
        
        def label_changed(event):
            
            self.selected_channel = event.value
            #self.update_slider()
            self.update_scatter()
        self.channel_menu.changed.connect(label_changed)
    
    def update_scatter(self):
       
       result = self.cell_data[self.cell_data[self.selected_channel] > self.threshold_value]
       self.x = np.array(result['centroid-0'])
       self.y = np.array(result['centroid-1'])
       self.points = np.stack((self.x, self.y)).transpose()
       self.points_layer = self.viewer.add_points(self.points, size=20)
        
    # Runs the program
    def run_all(self):
        self.image_path = self.get_working_directory("Select .czi image")
        
        im_data = []

        if str(self.image_path).endswith(".czi"):
            img = aicspylibczi.CziFile(str(self.image_path))
            dims = img.get_dims_shape()
            nchannels = dims[0]['C'][1]
            
            for c in np.arange(nchannels):
                print("Loading channel " + str(c))
                im_data.append(img.read_mosaic(C=c,scale_factor=1))
        
            self.im = np.stack(im_data)
        if self.image_path.name.startswith("IPICRC"):
            self.channels = ['DAPI','HLADR','CD8','CD163','CD4','XCR1','CD3','PDL1','EPCAM']
        else:
            self.channels = ['DAPI','HLADR','CD8','CD163','CD4','XCR1','CD3','PDL1','PanCK']
        self.launch_napari(self.im, self.channels, self.available_LUTs)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
