# Final Cell Segmentation

## Info 
This small project contains file to help segment czi files and display them with napari

# Getting Started

### Craete python venv and install requirements.txt
1. Create a python Virtual Env
   1. https://docs.python.org/3/library/venv.html
2. Install pip requirements
```
pip install -r requirements.txt
```
3. Download Ark-Analysis (https://github.com/angelolab/ark-analysis)
   1. Look at the README for instructions. This is necessary for the project

## Running Project

1. Get a .czi file and tile it using tile_czi.py
```
python3 tile_czi.py path/to/.czi/file
```
This should leave you with czi_filename_dir in this directory, which will be used for the next step

You can optionally add flags to argument parsing. Read comments at the top of the file for this information

2. Run create_deepcell_dir_format_from_single_channel_fovs.py on the output directory from tile_czi. 
```
python3 create_deepcell_dir_format_from_single_channel_fovs.py czi_filename_dir
```
This will set up the directory so it can be run for deepcell. 
Note: If you specify the target path to the ark-analysis folder it will automatically move it to ark-analysis/data folder. Optionally you can run: 
```
mv  -v cell_to_segment_dir LOCATION_OF_ARK_ANALYSIS/DATA
```
4. Copy the Segment_Image_Data-Final.ipynb in this directory to the scripts file of the ark-analysis project

5. Run Segment_Image_Data-Final.ipynb in the ark-analysis jupyter_notebook. This should give you a bunch of output in the cell_to_segment_dir. Move the completed directory over to the final_data directory here. 
```
mv  -v LOCATION_OF_ARK_ANALYSIS/DATA/cell_to_segent_dir  ./final_data/
```
6. Run view_main.py to use napari. 
   1. For loading images: Use your initial czi files
   2. For loading points: look for /single_cell_output/cell_table_arcsinh_transformed_stitched-{cell_file_name}.csv
   3. For loading boundaries: Look for segmentation_borders_current-{cell_file_name}.npy
   4. Also note: you can load a completely overlayed and stitched version of the cell called final_overlay-cell_name.tiff

