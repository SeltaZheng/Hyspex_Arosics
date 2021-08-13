"""
Auxiliary file generation for hyspex data.
All auxiliary files are stored in one obs_ort file.
Modified based on hytools_prep.py
"""

import os
# import argparse

import json

from osgeo import gdal, osr
import glob

import numpy as np

from ENVI import empty_envi_header, read_envi_header, write_envi_header

'''
    Adapted from code by Zhiwei Ye .

    TODO:
     - Write ancillary products to disk

'''


### ----------------------------------------


def write_to_envi(out_fn, arr, srs, gt, *, descr='', nodata=-9999.0):

    hdr = empty_envi_header()
    hdr['file type'] = 'ENVI Standard'
    hdr['description'] = descr

    hdr['samples'] = arr.shape[1]
    hdr['lines'] = arr.shape[0]
    try:
        hdr['bands'] = arr.shape[2]
        # shift arr:
        arr = np.moveaxis(arr, -1, 0)
    except IndexError:
        hdr['bands'] = 1

    hdr['byte order'] = 0
    hdr['header offset'] = 0
    hdr['interleave'] = 'bsq'
    hdr['data type'] = 4

    hdr['data ignore value'] = nodata

    hdr['coordinate system string'] = srs.ExportToWkt()

    hdr['map info'] = [
        srs.GetAttrValue('projcs').replace(',', ''),
        1, 1, gt[0], gt[3], gt[1], gt[1], ' ', ' ',
        srs.GetAttrValue('datum').replace(',', ''),
        srs.GetAttrValue('unit')
    ]

    if srs.GetAttrValue('PROJECTION').lower() == 'transverse_mercator':
        hdr['map info'][7] = srs.GetUTMZone()
        if gt[3] > 0.0:
            hdr['map info'][8] = 'North'
        else:
            hdr['map info'][8] = 'South'

    write_envi_header(out_fn + '.hdr', hdr)

    # Write data to binary file

    with open(out_fn, mode='wb+') as f:
        f.write(arr.tobytes())


def main(nodata=-9999.0):
    data_dir = r'Z:\townsenduser-rw\HyspexPro\Output\Cheesehead_V3\CHEESEHEAD_20190806\Merge'
    anc_dir = data_dir
    # glob the files:
    imgs = glob.glob(data_dir + '/*_Refl')
    for img in imgs:

        # data_dir = Path(os.getcwd()) / basename
        print(img)
        basename = os.path.basename(img)
        basename = basename.replace('_Refl', '')
        msk_fn = f'{data_dir}/{basename}_BackgroundMask'
        dem_fn = f'{data_dir}/{basename}_DEM'
        sca_fn = f'{data_dir}/{basename}_SCA'

        # Read the background mask
        msk_ds = gdal.Open(str(msk_fn))
        msk = msk_ds.GetRasterBand(1).ReadAsArray().astype(bool)

        # While the dataset is open, grab the geotransform array
        gt = msk_ds.GetGeoTransform()

        msk_ds = None

        # Calculate DEM aspect
        aspect_option = gdal.DEMProcessingOptions(trigonometric=False, computeEdges=True)
        aspect_ds = gdal.DEMProcessing(f'/vsimem/_tmp_{basename}_Aspect.tif',
                                       str(dem_fn), 'aspect', options=aspect_option)

        # Calculate DEM slope
        slope_option = gdal.DEMProcessingOptions(computeEdges=True)
        slope_ds = gdal.DEMProcessing(f'/vsimem/_tmp_{basename}_Slope.tif',
                                      str(dem_fn), 'slope', options=slope_option)

        # Convert slope & aspect to radians
        # dem_aspect = np.radians(aspect_ds.GetRasterBand(1).ReadAsArray())
        # dem_slope = np.radians(slope_ds.GetRasterBand(1).ReadAsArray())

        dem_aspect = aspect_ds.GetRasterBand(1).ReadAsArray()
        dem_slope = slope_ds.GetRasterBand(1).ReadAsArray()


        # Close temporary files
        aspect_ds = None
        slope_ds = None

        # Load SCA image header
        sca_hdr = read_envi_header(sca_fn + '.hdr')

        # Read SCA image
        sca_shape = (sca_hdr['bands'], sca_hdr['lines'], sca_hdr['samples'])
        sca_image = np.memmap(sca_fn, dtype='float32', mode='r', shape=sca_shape)

        # Get sensor viewing angles from SCA image
        # view_zenith = np.radians(sca_image[0])
        # view_azimuth = np.radians(sca_image[1])
        view_zenith = np.copy(sca_image[0])
        view_azimuth = np.copy(sca_image[1])

        # Calculate solar angles
        # solar_zenith = np.radians(float(sca_hdr['sun zenith']))
        # solar_azimuth = np.radians(float(sca_hdr['sun azimuth']))
        solar_zenith = float(sca_hdr['sun zenith'])
        solar_azimuth = float(sca_hdr['sun azimuth'])

        relative_azimuth = np.radians(dem_aspect) - np.radians(solar_azimuth)

        # Generate flat SZA image
        sza_image = np.full_like(view_zenith, nodata)
        sza_image[~msk] = solar_zenith

        # Generate flat SAA image
        saa_image = np.full_like(view_azimuth, nodata)
        saa_image[~msk] = solar_azimuth

        # Cosine of I
        cos_i = (np.cos(np.radians(solar_zenith)) * np.cos(np.radians(dem_slope)) +
                 np.sin(np.radians(solar_zenith)) * np.sin(np.radians(dem_slope)) * np.cos(relative_azimuth))

        # Phase
        # Wanner et al. JGRA 1995, eq. 51
        # Schlapfer et al. IEEE-TGARS 2015, eq. 2
        cos_phase = np.cos(np.radians(solar_zenith)) * np.cos(np.radians(view_zenith)) + np.sin(np.radians(solar_zenith)) * np.sin(np.radians(view_zenith)) * np.cos(
            relative_azimuth)
        phase = np.degrees(np.arccos(cos_phase))

        # Set no data values
        dem_aspect[msk] = nodata
        dem_slope[msk] = nodata
        view_zenith[msk] = nodata
        view_azimuth[msk] = nodata
        relative_azimuth[msk] = nodata
        cos_i[msk] = nodata

        # stack the obs:
        obs = np.stack([view_azimuth, view_zenith, saa_image, sza_image, phase, dem_slope, dem_aspect, cos_i], axis=2)

        # Get SRS
        srs = osr.SpatialReference()
        srs.ImportFromWkt(sca_hdr['coordinate system string'])

        ### TODO: Write ancillary products to disk, INCLUDING HEADERS

        anc_fn = f'{anc_dir}/{basename}_obs_ort'
        descr = ['sensor_az', 'sensor_zn', 'solar_az', 'solar_zn', 'phase', 'slope', 'aspect', 'cosine_i']
        # Write obs file:
        write_to_envi(anc_fn, obs, srs, gt,
                      descr=descr)


# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument('-b', '--basename')
#
#     args = parser.parse_args()
#
#     main(args.basename, nodata=-9999.0)
#

if __name__ == '__main__':
    main(nodata=-9999.0)
