import os
import argparse

import json

from pathlib import Path
from osgeo import gdal, osr

import numpy as np

from HyPro_CHTC.ENVI import empty_envi_header, read_envi_header, write_envi_header

'''
    Adapted from code by Zhiwei Ye .
    
    TODO:
     - Write ancillary products to disk
    
'''

### ----------------------------------------

def write_to_envi(out_fn, arr, srs, gt, *, descr='', nodata=-9999.0):
    
    # Write data to binary file
    with out_fn.open(mode='wb+') as f:
        f.write(arr.tobytes())
    
    hdr = empty_envi_header()
    hdr['file type'] = 'ENVI Standard'
    hdr['description'] = descr
    
    hdr['samples'] = arr.shape[1]
    hdr['lines'] = arr.shape[0]
    try:
        hdr['bands'] = arr.shape[2]
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
        1, 1, gt[0], gt[3], gt[1], gt[1], ' ',' ',
        srs.GetAttrValue('datum').replace(',', ''),
        srs.GetAttrValue('unit')
    ]
    
    
    if srs.GetAttrValue('PROJECTION').lower() == 'transverse_mercator':
        hdr['map info'][7] = srs.GetUTMZone()
        if gt[3]>0.0:
            hdr['map info'][8] = 'North'
        else:
            hdr['map info'][8] = 'South'
    
    write_envi_header(str(out_fn.with_suffix('.hdr')), hdr)
    

def main(basename, nodata=-9999.0):
    
    data_dir = Path(os.getcwd())/basename
    print(data_dir)
    
    msk_fn = data_dir/f'{basename}_BackgroundMask'
    dem_fn = data_dir/f'{basename}_DEM'
    sca_fn = data_dir/f'{basename}_SCA'
    
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
    dem_aspect = np.radians(aspect_ds.GetRasterBand(1).ReadAsArray())
    dem_slope = np.radians(slope_ds.GetRasterBand(1).ReadAsArray())

    # Close temporary files
    aspect_ds = None
    slope_ds = None

    # Load SCA image header
    sca_hdr = read_envi_header(sca_fn.with_suffix('.hdr'))

    # Read SCA image
    sca_shape = (sca_hdr['bands'], sca_hdr['lines'], sca_hdr['samples'])
    sca_image = np.memmap(sca_fn, dtype='float32', mode='r', shape=sca_shape)
    
    # Get sensor viewing angles from SCA image
    view_zenith =  np.radians(sca_image[0])
    view_azimuth = np.radians(sca_image[1])
    
    # Calculate solar angles
    solar_zenith = np.radians(float(sca_hdr['sun zenith']))
    solar_azimuth = np.radians(float(sca_hdr['sun azimuth']))
    
    relative_azimuth = dem_aspect-solar_azimuth
    
    # Generate flat SZA image
    sza_image = np.full_like(view_zenith, nodata)
    sza_image[~msk] = solar_zenith

    # Generate flat SAA image
    saa_image = np.full_like(view_azimuth, nodata)
    saa_image[~msk] = solar_azimuth
    
    # Cosine of I
    cos_i = (np.cos(solar_zenith)*np.cos(dem_slope) + 
             np.sin(solar_zenith)*np.sin(dem_slope) * np.cos(relative_azimuth))

    # Phase
    # Wanner et al. JGRA 1995, eq. 51
    # Schlapfer et al. IEEE-TGARS 2015, eq. 2
    cos_phase = np.cos(solar_zenith)*np.cos(view_zenith) + np.sin(solar_zenith)*np.sin(view_zenith)*np.cos(relative_azimuth)
    phase = np.arccos(cos_phase)
    
    # Set no data values
    dem_aspect[msk] = nodata
    dem_slope[msk] = nodata
    view_zenith[msk] = nodata
    view_azimuth[msk] = nodata
    relative_azimuth[msk] = nodata
    cos_i[msk] = nodata

    # Get SRS
    srs = osr.SpatialReference()
    srs.ImportFromWkt(sca_hdr['coordinate system string'])
    
    ### TODO: Write ancillary products to disk, INCLUDING HEADERS
    
    anc_dir = Path(f'./{basename}_ancillary')
    anc_dir.mkdir()
    
	# Write DEM slope
    write_to_envi(anc_dir/f'{basename}_DSMSlope', dem_slope, srs, gt,
                  descr='DEM slope, in [rad]')
    
    # Write DSM aspect
    write_to_envi(anc_dir/f'{basename}_DSMAspect', dem_aspect, srs, gt,
                  descr='DEM aspect, in [rad]')
    
    # Write sensor zenith angles
    write_to_envi(anc_dir/f'{basename}_VZA', view_zenith, srs, gt,
                  descr='Sensor viewing zenith angle, in [rad]')
    
    # Write sensor azimuth angles
    write_to_envi(anc_dir/f'{basename}_VAA', view_azimuth, srs, gt,\
                  descr='Sensor viewing azimuth angle, in [rad]')
    
    # Write solar zenith angle
    write_to_envi(anc_dir/f'{basename}_SZA', sza_image, srs, gt,
                  descr='Solar zenith angle, in [rad]')
    
    # Write solar azimuth angle
    write_to_envi(anc_dir/f'{basename}_SAA', saa_image, srs, gt,\
                  descr='Solar azimuth angle, in [rad]')
	
	# Write relative azimuth angles
    write_to_envi(anc_dir/f'{basename}_RAA', relative_azimuth, srs, gt,\
                  descr='Relative viewing azimuth angle, in [rad]')
    
	# Write cosine I
    write_to_envi(anc_dir/f'{basename}_CosI', cos_i, srs, gt,\
                  descr='Cosine of I')
    
	# Write phase
    write_to_envi(anc_dir/f'{basename}_Phase', phase, srs, gt,\
                  descr='Phase')
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--basename')
    
    args = parser.parse_args()
    
    main(args.basename, nodata=-9999.0)
