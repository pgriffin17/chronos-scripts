#!/usr/bin/env python

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import sys, glob, os
import astropy.units as u
import astropy.wcs as wcs
from astropy.coordinates import SkyCoord
import click
from tabulate import tabulate

# suppress warnings
import warnings

warnings.filterwarnings("ignore")

__author__ = "Philip Griffin"

def cutout(data, x, y, size):
    """
    Create a cutout around a point in an image, with a specified size
    as the x and y dimensions of the cutout.

    Parameters
    ----------
    data : numpy.ndarray
        The image data

    x : int
        The x coordinate of the center of the cutout

    y : int
        The y coordinate of the center of the cutout

    size : int
        The size of the cutout in pixels

    Returns
    -------
    numpy.ndarray
        The cutout image data
    """
    size = int(size / 2) # Cut size in half to be the half-size of the cutout
    return data[y - size : y + size, x - size : x + size]


def check_saturation(image_path, cutout_size=40, saturation_threshold=50000):
    """
    Check if a star is saturated in an image

    Parameters
    ----------
    image_path : str
        Path to the image to check

    cutout_size : int
        Size of the cutout to take around the star in pixels

    saturation_threshold : int
        Threshold to consider a star as saturated

    Returns
    -------
    is_saturated : bool
        True if the star is saturated, False otherwise
    
    max_val : int
        Maximum value in the cutout around the star

    exptime : float
        Exposure time of the image

    filter_name : str
        Name of the filter used in the image
    """
    # print(f"Checking saturation for {image_path}")
    # Open the image and make wcs object
    hdul = fits.open(image_path)
    wcs_obj = wcs.WCS(hdul[0].header)

    # Get BLKRA from the header
    blkra = hdul[0].header["BLKRA"]
    # Get BLKDEC from the header
    blkdec = hdul[0].header["BLKDEC"]  

    # Get filter name
    filter_name = hdul[0].header["FILTER"]
    # Remove whitespace from filter name
    filter_name = filter_name.replace(" ", "")
    # Get exposure time
    exptime = hdul[0].header["EXPTIME"]

    # Store in skycoord object
    c = SkyCoord(blkra, blkdec, unit=(u.hourangle, u.degree))

    # Get pixel coordinates of the skycoord object
    pix = wcs_obj.world_to_pixel(c)

    # Get the data from the image
    data = hdul[0].data

    # Get the cutout around the star
    cutout_data = cutout(data, int(pix[0]), int(pix[1]), cutout_size)

    # Close the image
    hdul.close()

    # Find the maximum value in the cutout
    max_val = np.max(cutout_data)

    # Check if the cutout is saturated using the threshold
    is_saturated = max_val > saturation_threshold

    return is_saturated, max_val, exptime, filter_name


@click.command()
@click.option(
    "--image_path", "-i", help="Path to the image to check", required=True, type=str
)
@click.option(
    "--cutout_size",
    "-s",
    help="Size of the cutout to take around the star in pixels",
    default=40,
)
@click.option(
    "--saturation_threshold",
    "-t",
    help="Threshold to consider a star as saturated",
    default=50000,
)
def main(image_path, cutout_size, saturation_threshold):
    print(f"Checking saturation for images in {image_path} with cutout size {cutout_size}x{cutout_size} px and saturation threshold {saturation_threshold} ADU.")
    file_paths = glob.glob(image_path)
    # Sort the file paths
    file_paths.sort()
    non_grism_file_paths = []
    is_sat_list = []
    max_val_list = []
    exptime_list = []
    filter_name_list = []
    for image_path in file_paths:
        if not os.path.exists(image_path):
            print(f"File {image_path} does not exist")
            continue
        elif "_lrg_" in image_path or "_hrg_" in image_path:
            # Name contains _lrg_ or _hrg_, skip as it is a grism image
            # TODO: Add a check for grism images - extract the strip of the grism image and check for saturation
            continue
        is_saturated, max_val, exptime, filter_name = check_saturation(image_path, cutout_size, saturation_threshold)
        image_name = os.path.basename(image_path)
        non_grism_file_paths.append(image_path)
        is_sat_list.append(is_saturated)
        max_val_list.append(max_val)
        exptime_list.append(exptime)
        filter_name_list.append(filter_name)

        # if is_saturated:
        #     print(f"{image_name} is saturated with max value {max_val}")
        # else:
        #     print(f"{image_name} is not saturated with max value {max_val}")

    # Create better table that can be output to a terminal using tabulate
    
    table = []
    print("\n")
    for i in range(len(non_grism_file_paths)):
        sat_label = "Sat" if is_sat_list[i] else "Unsat"
        table.append([os.path.basename(non_grism_file_paths[i]), filter_name_list[i], exptime_list[i], sat_label, max_val_list[i]])
    print(tabulate(table, headers=["Image Name", "Fil", "Exp", "Sat?", "Peak"]))
    print("\n")

if __name__ == "__main__":
    main()
