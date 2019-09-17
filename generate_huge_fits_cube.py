#!/usr/bin/python3
"""
This script generates an empty dummy 4 dimensional data cube in fits format.
After the initialisation this cube gets filled with fits image data. The
cube header gets updated from the first image in PATHLIST_STOKESI.
This script can be used to generate fits data cubes of sizes that exceeds the
machine's RAM (tested with 234 GB RAM and 335 GB cube data).

Please adjust the INPUT section in this script to your needs.

The data in directory `images` is test data and consists of Gaussian noise only.  

Developed at: IDIA (Institure for Data Intensive Astronomy), Cape Town, ZA
Inspired by: https://github.com/idia-astro/image-generator
Source: https://github.com/lh-astro/IDIA-generate-huge-fits-cube
E-Mail: astro[et]lennartheino[.]de
"""

import itertools
import logging
from logging import info, error
import os
import csv
import datetime
from glob import glob

import numpy as np
from astropy.io import fits


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# INPUT

# directory with fits images per channel for Stokes IQUV.
DIR_IMAGES = "images/"

PATHLIST_STOKESI = sorted(glob(DIR_IMAGES + "*.I.im-image.fits"))
PATHLIST_STOKESQ = sorted(glob(DIR_IMAGES + "*.Q.im-image.fits"))
PATHLIST_STOKESU = sorted(glob(DIR_IMAGES + "*.U.im-image.fits"))
PATHLIST_STOKESV = sorted(glob(DIR_IMAGES + "*.V.im-image.fits"))

OBJECT_NAME = os.path.basename(PATHLIST_STOKESI[0].split(".")[0])
CUBE_NAME = "cube." + OBJECT_NAME + ".fits"

# flag data that is above the RMS noise estimate from Stokes V. Set is to a high
# value for no flagging
RMS_THRESHOLD = 20  # in [uJy/beam]
RMS_THRESHOLD = RMS_THRESHOLD * 1e-6  # to [Jy/beam], don't change this!

# Outputs a statistics file with estimates for RMS noise in Stokes I and V
WRITE_STATISTICS_FILE = True
FILEPATH_STATISTICS = "rms-noise-statistics.tab"

# INPUT
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# SETTINGS

logging.basicConfig(
    format="%(asctime)s\t[ %(levelname)s ]\t%(message)s", level=logging.INFO
)
SEPERATOR = "-----------------------------------------------------------------"

# SETTINGS
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def get_and_add_custom_header(header):
    """
    Gets header from fits file and updates the cube header.


    Parameters
    ----------
    header: astroy.io.fits header
       The header class that gets updated

    Returns
    -------
    header: astroy.io.fits header
       The header class that was updated

    """
    info(SEPERATOR)
    info("Getting header for data cube from: %s", PATHLIST_STOKESI[0])
    with fits.open(PATHLIST_STOKESI[0], memmap=True) as hud:
        header = hud[0].header
        # Optional: Update the header.
        header["OBJECT"] = OBJECT_NAME
        header["NAXIS3"] = len(PATHLIST_STOKESI)
        header["CTYPE3"] = ("FREQ", "")
    return header


def make_empty_image():
    """
    Generate an empty dummy fits data cube.

    The data cube dimensions are derived from the channel fits images. The
    resulting data cube can exceed the machine's RAM.

    """
    info(SEPERATOR)
    info("Getting image dimension for data cube from: %s", PATHLIST_STOKESI[0])
    with fits.open(PATHLIST_STOKESI[0], memmap=True) as hud:
        xdim, ydim = np.squeeze(hud[0].data).shape
    info("X-dimension: %s", xdim)
    info("Y-dimension: %s", ydim)

    info(
        "Getting channel dimension Z for data cube from number of entries in PATHLIST_STOKESI."
    )
    zdim = len(PATHLIST_STOKESI)
    info("Z-dimension: %s", zdim)

    info("Assuming full Stokes for dimension W.")
    wdim = 4
    info("W-dimension: %s", wdim)

    dims = tuple([xdim, ydim, zdim, wdim])

    # create header

    dummy_dims = tuple(1 for d in dims)
    dummy_data = np.zeros(dummy_dims, dtype=np.float32)
    hdu = fits.PrimaryHDU(data=dummy_data)

    header = hdu.header
    header = get_and_add_custom_header(header)
    for i, dim in enumerate(dims, 1):
        header["NAXIS%d" % i] = dim

    header.tofile(CUBE_NAME, overwrite=True)

    # create full-sized zero image

    header_size = len(
        header.tostring()
    )  # Probably 2880. We don't pad the header any more; it's just the bare minimum
    data_size = np.product(dims) * np.dtype(np.float32).itemsize
    # This is not documented in the example, but appears to be Astropy's default behaviour
    # Pad the total file size to a multiple of the header block size
    block_size = 2880
    data_size = block_size * ((data_size // block_size) + 1)

    with open(CUBE_NAME, "rb+") as f:
        f.seek(header_size + data_size - 1)
        f.write(b"\0")


def get_mad(a, axis=None):
    """
    Compute *Median Absolute Deviation* of an array along given axis.

    from: https://informatique-python.readthedocs.io/fr/latest/Exercices/mad.html

    Parameters
    ----------
    a: numpy.array
       The numpy array of which MAD gets calculated from

    Returns
    -------
    mad: float
       MAD from a

    """
    # Median along given axis, but *keeping* the reduced axis so that
    # result can still broadcast against a.
    med = np.median(a, axis=axis, keepdims=True)
    mad = np.median(np.absolute(a - med), axis=axis)  # MAD along given axis
    return mad


def get_std_via_mad(npArray):
    """
    Estimate standard deviation via Median Absolute Deviation.


    Parameters
    ----------
    npArray: numpy.array
       The numpy array of which the Standard Deviation gets calculated from

    Returns
    -------
    std: float
       Standard Deviation from MAD

    """
    mad = get_mad(npArray)
    std = 1.4826 * mad
    # std = round(std, 3)
    info("Got std via mad [uJy/beam]: %s ", round(std * 1e6, 2))
    return std


def check_rms(npArray):
    """
    Check if the Numpy Array is below RMS_THRESHOLD and above 1e-6 uJy/beam.

    If the Numpy Array is not within the range it gets assigned to not a number
    (np.nan).

    Parameters
    ----------
    npArray: numpy.array
       The numpy array to check

    Returns
    -------
    [npArray, std]: list with numpy.array and float
       List of length 2 with  the Numpy Array and the Standard Deviation

    """
    std = get_std_via_mad(npArray)
    if (std > RMS_THRESHOLD) or (std < 1e-6):
        npArray = np.nan
    return [npArray, std]


def write_statistics_file(rmsDict):
    """
    Takes the dictionary with Stokes I and V RMS noise and writes it to a file.

    Parameters
    ----------
    rmdDict: dict of lists with floats
       Dictionary with lists for Stokes I and V rms noise

    """
    legendList = ["RMS Stokes I [uJy/beam]", "RMS Sokes V [uJy/beam]"]
    info("Writing statistics file: %s", FILEPATH_STATISTICS)
    with open(FILEPATH_STATISTICS, "w") as csvFile:
        writer = csv.writer(csvFile, delimiter="\t")
        csvData = [legendList]
        for i, entry in enumerate(rmsDict["rmsI"]):
            rmsI = round(rmsDict["rmsI"][i] * 1e6, 4)
            rmsV = round(rmsDict["rmsV"][i] * 1e6, 4)
            csvData.append([rmsI, rmsV])
        writer.writerows(csvData)


def fill_cube_with_images():
    """
    Fills the empty data cube with fits data.


    """
    info(SEPERATOR)
    info("Opening data cube: %s", CUBE_NAME)
    hudCube = fits.open(CUBE_NAME, memmap=True, mode="update")
    dataCube = hudCube[0].data

    rmsDict = {}
    rmsDict["rmsI"] = []
    rmsDict["rmsV"] = []
    for i, filePathFits in enumerate(PATHLIST_STOKESI):
        # Switch
        stokesVflag = False

        info(SEPERATOR)
        info("Opening fits file: %s", PATHLIST_STOKESV[i])
        with fits.open(PATHLIST_STOKESV[i], memmap=True) as hud:
            checkedArray, std = check_rms(hud[0].data[:, :])
            rmsDict["rmsV"].append(std)
            dataCube[3, i, :, :] = checkedArray
            if np.isnan(np.sum(checkedArray)):
                stokesVflag = True

        if not stokesVflag:
            info("Opening fits file: %s", PATHLIST_STOKESI[i])
            with fits.open(PATHLIST_STOKESI[i], memmap=True) as hud:
                std = get_std_via_mad(hud[0].data[:, :])
                rmsDict["rmsI"].append(std)
                dataCube[0, i, :, :] = hud[0].data[:, :]

            info("Opening fits file: %s", PATHLIST_STOKESQ[i])
            with fits.open(PATHLIST_STOKESQ[i], memmap=True) as hud:
                dataCube[1, i, :, :] = hud[0].data[:, :]

            info("Opening fits file: %s", PATHLIST_STOKESU[i])
            with fits.open(PATHLIST_STOKESV[i], memmap=True) as hud:
                dataCube[2, i, :, :] = hud[0].data[:, :]

        if stokesVflag:
            info(
                "Stokes V RMS noise of %s [uJy/beam] above RMS_THRESHOLD of %s [uJy/beam]. Flagging Stokes IQUV",
                str(round(rmsDict["rmsV"][-1] * 1e6, 2)),
                str(round(RMS_THRESHOLD * 1e6, 3)),
            )
            dataCube[0, i, :, :] = np.nan
            dataCube[1, i, :, :] = np.nan
            dataCube[2, i, :, :] = np.nan
            dataCube[3, i, :, :] = np.nan

    hudCube.close()
    if WRITE_STATISTICS_FILE:
        write_statistics_file(rmsDict)


if __name__ == "__main__":
    # start timestamp
    info(SEPERATOR)
    TIMESTAMP_START = datetime.datetime.now()
    info("START script at: %s", TIMESTAMP_START)
    info(SEPERATOR)

    # call methods
    make_empty_image()
    fill_cube_with_images()

    # end timestamp
    TIMESTAMP_END = datetime.datetime.now()
    TIMESTAMP_DELTA = TIMESTAMP_END - TIMESTAMP_START
    info(SEPERATOR)
    info("END script at %s in %s", str(TIMESTAMP_END), str(TIMESTAMP_DELTA))
    info(SEPERATOR)
