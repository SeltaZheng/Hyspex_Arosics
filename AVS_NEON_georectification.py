"""
Testing script for Arosics
"""
from timeit import default_timer
import numpy as np
from ENVI import read_envi_header, write_envi_header
from arosics.DeShifter import DESHIFTER
from arosics import COREG_LOCAL
from osgeo import gdal

def hyspex_extract(src_img, bands, out_img):
    """
    extract certain bands from hyspex images (note: images are in BSQ)
    :param src_img:
    :param out_img:
    :param bands: list; bands to be extracted, starting from 0
    :return:
    """

    src_header = read_envi_header(src_img + '.hdr')
    # src_image = np.memmap(src_img,
    #                       dtype='int16',
    #                       mode='r',
    #                       offset=src_header['header offset'],
    #                       shape=(src_header['lines'],
    #                              src_header['samples'],
    #                              src_header['bands'],))

    # Create empty array for the output image:
    # out = np.memmap(out_img, dtype='int16',
    #                 mode='w+',
    #                 shape=(src_header['lines'],
    #                        src_header['samples'],
    #                        len(bands))
    #                 )
    # Read the given bands:
    i = 0
    with open(out_img, 'wb') as outfile:
        for b in bands:
            print('Now reading band {}'.format(b))
            src_image = np.memmap(src_img,
                                  dtype='int16',
                                  # dtype='float32',
                                  mode='r',
                                  shape=(src_header['lines'],
                                         src_header['samples']),
                                  offset=int(src_header['lines'] * src_header['samples'] * b * 16 / 8)) # 16 for int16

            # # Create empty array for the output image:
            # out = np.memmap(out_img, dtype='int16',
            #                 mode='w+',
            #                 shape=(src_header['lines'],
            #                        src_header['samples']),
            #                 offset=int(src_header['lines'] * src_header['samples'] * i * 16 / 8))

            # out[:, :, i] = src_image[:]
            outfile.write(src_image[:])
            del src_image

        i = i + 1
    outfile.close()
    # flush memory:
    # out.flush()

    # src_image.flush()

    # Modify the header for out_image:
    tgt_header = src_header
    # number of bands:
    tgt_header['bands'] = len(bands)
    # wavelength:
    wvl = [src_header['wavelength'][b] for b in bands]
    tgt_header['wavelength'] = wvl

    # fwhm:
    fwhm = [src_header['fwhm'][b] for b in bands]
    tgt_header['fwhm'] = fwhm
    write_envi_header(out_img + '.hdr', tgt_header)

def Gdal_write_band(dst_raster,band_number,band_array,band_name,no_data=-9999):
    '''
    Write band to a gdal raster.
    band_number: int
    band number in the dst raster
    band_array: np array
    the array needs to be written
    band_name: str
    name of the band
    no_data: int or float
    no data value, default to -9999
    '''

    outband = dst_raster.GetRasterBand(band_number)
    outband.SetNoDataValue(no_data)
    outband.SetDescription(band_name)
    outband.WriteArray(band_array)
    outband.FlushCache()

def img_band_extraction(src_img, band, out_img, nodata=-9999):
    src = gdal.Open(src_img)
    src_arr = src.GetRasterBand(band).ReadAsArray()
    driver = gdal.GetDriverByName('ENVI')
    dst_out = driver.Create(out_img, src_arr.shape[1], src_arr.shape[0], 1, gdal.GDT_Float32)
    dst_out.SetGeoTransform(src.GetGeoTransform())
    dst_out.SetProjection(src.GetProjection())
    Gdal_write_band(dst_out, 1, src_arr, 'green', no_data=nodata)




def modify_hyspex_band(img, band, arr):
    src_header = read_envi_header(img + '.hdr')
    # if int16:
    # src_image = np.memmap(img,
    #                       dtype='int16',
    #                       mode='r+',
    #                       shape=(src_header['lines'],
    #                              src_header['samples']),
    #                       offset=int(src_header['lines'] * src_header['samples'] * band * 16 / 8))
    # if float 32
    src_image = np.memmap(img,
                          dtype='float32',
                          mode='r+',
                          shape=(src_header['lines'],
                                 src_header['samples']),
                          offset=int(src_header['lines'] * src_header['samples'] * band * 32 / 8))
    src_image[:] = arr
    src_image.flush()
    del src_image, arr


def band_select(header):
    """
    find two bands (one from VNIR, one from SWIR) that closest in wavelength
    :param header: dict.
    :return: band numbers (starting from 0)
    """
    wvl = np.array(header['wavelength'])
    fwhm = np.array(header['fwhm'])
    dif = wvl[1:] - wvl[0:-1]
    dif = abs(dif)
    idx = np.where(dif == np.min(dif))
    idx = idx[0][0]
    idx_vnir = np.where(fwhm < 4)[0]
    idx_swir = np.where(fwhm > 4)[0]
    if np.isin(idx, idx_vnir) and np.isin((idx + 1), idx_swir):
        band_swir = idx + 1
        band_vnir = idx
    elif np.isin((idx + 1), idx_vnir) and np.isin(idx, idx_swir):
        band_swir = idx
        band_vnir = idx + 1
    else:
        print('Both bands are in vnir or in swir')
        band_vnir = -1
        band_swir = -1

    return band_vnir, band_swir, idx_vnir, idx_swir
def Arosics_test_local(ref_img, tgt_img, kwargs):
    """
    testing function for local registration in Arosics
    :param ref_img: str; full path and filename of the reference image
    :param tgt_img: str; full path and filename of the target image
    :param kwargs: dict; dictionary containing the detailed arguments for local correction.
    :return:
    """

    from arosics import COREG_LOCAL
    start = default_timer()
    CRL = COREG_LOCAL(ref_img, tgt_img, **kwargs)
    CRL.correct_shifts()
    end = default_timer()
    print(end - start)
    return CRL

def Arosics_test_global(ref_img, tgt_img, kwargs):
    """
    testing function for local registration in Arosics
    :param ref_img: str; full path and filename of the reference image
    :param tgt_img: str; full path and filename of the target image
    :param kwargs: dict; dictionary containing the detailed arguments for local correction.
    :return:
    """
    from arosics import COREG
    start = default_timer()
    CR = COREG(ref_img, tgt_img, **kwargs)
    CR.calculate_spatial_shifts()
    CR.correct_shifts()
    end = default_timer()
    print(end - start)
    return CR

##-------------- main function ------------------------------
def main():
    # dir_avs = r'D:\CA_NEON\AVIRIS\Warped'
    # dir_NEON = r'D:\CA_NEON\NEON\Merged'
    # dest_dir = dir_avs + '/Deshift'
    dir_NEON = r'Z:\townsenduser-rw\CA_project\CA_NEON\NEON\Merged'
    dir_in = r'Z:\townsenduser-rw\CA_project\CA_NEON\CA_spectra_unmixing\testing'
    # img = '/CHEESEHEAD_20190806_03_Refl'
    # ref_img = dir_in + '/Merge' + img
    # dest_file = dest_dir + img
    # img = '/AVIRIS_TEAK_2017_solar4_warp'
    # header = read_envi_header(dir_avs + '/AVIRIS_TEAK_2017_solar4_warp.hdr')
    # band_vnir, band_swir, idx_vnir, idx_swir = band_select(header)
    # hyspex_extract(dir_in + '/Merge' + img, [186], dest_dir + img + '_refband')
    # hyspex_extract(dir_avs + '/AVIRIS_TEAK_2017_solar4_warp', [16], dir_avs + img + '_g')
    # img_band_extraction(dir_avs+img, 17, dir_avs + img + '_g', nodata=-9999)
    # ref_img = dir_NEON + '/TEAK_rgb_2017_v2'
    # tgt_img = dir_avs + img + '_g'
    ref_img = dir_NEON + '/TEAK_rgb_2017_v2'
    tgt_img = dir_in + '/f170607t01p00r16_rfl_v1g_img_topo_brdf_solar_g_TEAK_clipped'
    dest_img = dir_in + '/f170607t01p00r16_rfl_v1g_img_topo_brdf_solar_SVDI_v0_TEAK_clipped'
    dest_img_deshift = dest_img + '_deshift'
    # out = dir_in + '/shifted_global'
    grid_ls = [150]#[150, 200, 300, 350]
    windows = [70]#[70, 120, 128, 150]
    band_idx = [0, 1, 2, 3, 4]
    bandnames = ['Substrate', 'Canopy', 'Dark', 'Snow', 'RMSE']

    for i in range(len(grid_ls)):

        # out = dir_avs + '/shifted_Teak_2017_{}_{}'.format(grid_ls[i], windows[i])
        out = dir_in +'/shifted_Teak_2017_{}_{}'.format(grid_ls[i], windows[i])
        g = grid_ls[i]
        w = windows[i]
        kwargs = kwargs = {
            'grid_res': g,
            'window_size': (w, w),
            'path_out': out,
            'nodata': (-9999, -9999),
            'r_b4match': 2,
            's_b4match': 1,
            'max_iter': 10,
            'max_shift': 5,
            'min_reliability': 25,
            # 'tieP_filter_level': 1,
            'resamp_alg_calc': 'nearest',
            'q': False
        }
        CRL = COREG_LOCAL(ref_img, tgt_img, **kwargs)
        CRL.calculate_spatial_shifts()
        # CRL.correct_shifts()


        for j, band in enumerate(band_idx):

            result = DESHIFTER(dest_img, CRL.coreg_info,
                               band2process=int(band+1),
                               align_grids=True,
                               match_gsd=True,
                               nodata=0,
                               resamp_alg='cubic',
                               q=True).correct_shifts()
            # modify_hyspex_band(dest_img_deshift, band, result['arr_shifted'])
            # Use gdal if the deshifted array has different dimensions compared to the original one.
            if j == 0:
                driver = gdal.GetDriverByName('ENVI')

                dst_raster = driver.Create(dest_img_deshift, result['arr_shifted'].shape[1], result['arr_shifted'].shape[0], len(band_idx), gdal.GDT_Float32)
                dst_raster.SetGeoTransform(result['updated geotransform'])
                dst_raster.SetProjection(result['updated projection'])
            Gdal_write_band(dst_raster, int(band+1), result['arr_shifted'], bandnames[j], no_data=0)
            del result

        # calculate the shift info
        # try:
        #     CRL = COREG_LOCAL(ref_img, tgt_img, **kwargs)
        #     CRL.calculate_spatial_shifts()
        #     CRL.correct_shifts()
        # except:
        #     continue

        # # apply the shift to a single band:
        # dest_img = dest_file + '_deshift'
        # band = 185
        # # shift_arg = {
        # #         'band2process': int(band + 1),
        # #         'resamp_alg': 'nearest',
        # #         'nodata': 0,
        # #         'align_grids': True,
        # #         'q': True
        # #     }
        # result = DESHIFTER(dir_in + '/Merge' + img, CRL.coreg_info,
        #                    band2process=int(band+1),
        #                    align_grids=True,
        #                    match_gsd=True,
        #                    nodata=0,
        #                    resamp_alg='cubic',
        #                    q=True).correct_shifts()
        # modify_hyspex_band(dest_img, band, result['arr_shifted'])

    # # apply the shift to all swir bands:
    # for band in idx_swir:
    #     print('Now shifting band {}'.format(band))
    #     shift_arg = {
    #         'band2process': int(band + 1),
    #         'resamp_alg': 'nearest',
    #         'nodata': 0,
    #         'q': True
    #     }
    #     try:
    #         result = DESHIFTER(ref_img, CRL.coreg_info, **shift_arg).correct_shifts()
    #         modify_hyspex_band(dest_img, band, result['arr_shifted'])
    #         del result
    #     except:
    #         print('{} is not processed properly'.format(img))
    #         break
    #
    # del CRL

 # # If using the global correction
 #    kwargs = {
 #        'path_out': out,
 #        'r_b4match': 2,
 #        's_b4match': 1,
 #        'max_shift': 5,
 #        'nodata': (-9999, -9999),
 #        'resamp_alg_calc': 'nearest',
 #        'q': False
 #    }
 #    Arosics_test_global(ref_img, tgt_img, kwargs)
    # return CRL.coreg_info


if __name__ == '__main__':
    main()
