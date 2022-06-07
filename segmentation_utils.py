import argparse
import pathlib

#Useful utilities for segmentation folder
# Includes
#   1. Argparser for taking different command line arguments
#   2. A print_colored function to print in color to help readability and debug

cell_segment_parser = argparse.ArgumentParser(description='Process some integers.')
cell_segment_parser.add_argument(
    '--debug', '-d', dest='debug',
    action='store_true',
    help='Print additional information to the terminal when running script'
)

cell_segment_parser.add_argument(
    'files', nargs='+', type=pathlib.Path,
    help='Print additional information to the terminal when running script'
)
class colors:
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    UNDERLINE = '\033[4m'
    RESET = '\033[0m'
    BOLD = '\033[1m'

colors_converter_dict = {
    'black': colors.BLACK,
    'red': colors.RED,
    'r': colors.RED,
    'green': colors.GREEN,
    'g': colors.GREEN,
    'yellow': colors.YELLOW,
    'y': colors.YELLOW,
    'blue': colors.BLUE,
    'magenta': colors.MAGENTA,
    'm': colors.MAGENTA,
    'cyan': colors.CYAN,
    'c': colors.CYAN,
    'white': colors.WHITE,
    'w': colors.WHITE,
    'underline': colors.UNDERLINE,
    'u': colors.UNDERLINE,
    'bold': colors.UNDERLINE,
}
def print_colored(output_color, *args):
    output_color=output_color.lower()
    if(output_color not in colors_converter_dict):
        print(colors.RED + f"Error: {output_color} is not an acceptable color option. Will print out in cyan!" + colors.RESET)
        output_color = 'cyan'
    print(colors_converter_dict[output_color], end='')
    print(*args)
    print(colors.RESET, end='')

napari_viewer_parser = argparse.ArgumentParser(description='Parse arguments for napari viewer.')

napari_viewer_parser.add_argument(
    '--debug', '-d', dest='debug',
    action='store_true',
    help='Print additional information to the terminal when running script'
)

napari_viewer_parser.add_argument("--points", "-p",
                    dest='points',
                    action="store",
                    nargs=1,
                    type=pathlib.Path,
                    help="Specify a single file with points label to preload in");

napari_viewer_parser.add_argument("--image", "-i",
                    dest='image',
                    action="store",
                    nargs=1,
                    type=pathlib.Path,
                    help="Specify a single czi image file to preload in");
#
#
napari_viewer_parser.add_argument("--bounds", "-b",
                    dest='bounds',
                    action="store",
                    nargs=1, #Show it is optional
                    type=pathlib.Path,
                    help="Specify a single file with segementation boundaries to load in");