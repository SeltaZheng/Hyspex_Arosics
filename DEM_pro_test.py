"""
Test out gdal DEM processing with the 'edge' option
"""

import gdal

dir_in = r'Z:\townsenduser-rw\HyspexPro\Output\Cheesehead_V3\CHEESEHEAD_20190629\Merge'
dem_fn = '/CHEESEHEAD_20190629_03_DEM'

aspect_option = gdal.DEMProcessingOptions(trigonometric=False, computeEdges=True)
gdal.DEMProcessing(dir_in+'/aspect_03.tif', dir_in+dem_fn, 'aspect', options=aspect_option)