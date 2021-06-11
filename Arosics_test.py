"""
Testing script for Arosics
"""
from timeit import default_timer
import numpy as np
from ENVI import read_envi_header, write_envi_header
from arosics.DeShifter import DESHIFTER
from arosics import COREG_LOCAL

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
                                  mode='r',
                                  shape=(src_header['lines'],
                                         src_header['samples']),
                                  offset=int(src_header['lines'] * src_header['samples'] * b * 16 / 8))

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


def modify_hyspex_band(img, band, arr):
    src_header = read_envi_header(img + '.hdr')
    src_image = np.memmap(img,
                          dtype='int16',
                          mode='r+',
                          shape=(src_header['lines'],
                                 src_header['samples']),
                          offset=int(src_header['lines'] * src_header['samples'] * band * 16 / 8))
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
    # dir_in = r'Z:\townsenduser-rw\HyspexPro\Output\Cheesehead_V3\CHEESEHEAD_20190626\CHEESEHEAD_20190626_02_quicklooks'
    dir_in = r'Z:\townsenduser-rw\HyspexPro\Output\Cheesehead_V3\CHEESEHEAD_20190806'
    dest_dir = dir_in + '/Deshift'
    img = '/CHEESEHEAD_20190806_03_Refl'
    # ref_img = dir_in + '/Merge' + img
    dest_file = dest_dir + img
    header = read_envi_header(dir_in + '/Merge' + img + '.hdr')
    band_vnir, band_swir, idx_vnir, idx_swir = band_select(header)
    # hyspex_extract(dir_in + '/Merge' + img, [186], dest_dir + img + '_refband')
    # hyspex_extract(dir_in + '/Merge' + img, [185], dest_dir + img + '_swir_test')

    ref_img = dest_dir + img + '_refband'
    tgt_img = dest_dir + img + '_swir_test'

    grid_ls = [350]#[300, 350, 400, 450]
    windows = [150]#[128, 150, 264, 280]

    for i in range(len(grid_ls)):

        out = dest_dir + '/shifted_local_03_{}_{}'.format(grid_ls[i], windows[i])
        g = grid_ls[i]
        w = windows[i]
        kwargs = kwargs = {
            'grid_res': g,
            'window_size': (w, w),
            'path_out': out,
            'nodata': (0, 0),
            'r_b4match': 1,
            's_b4match': 1,
            'max_iter': 10,
            'max_shift': 10,
            # 'min_reliability': 25,
            # 'tieP_filter_level': 1,
            'resamp_alg_calc': 'nearest',
            'q': False
        }


        dest_img = dest_file + '_deshift'

        # calculate the shift info
        try:
            CRL = COREG_LOCAL(ref_img, tgt_img, **kwargs)
            CRL.calculate_spatial_shifts()
            # CRL.correct_shifts()
        except:
            continue

        # apply the shift to a single band:
        band = 185
        # shift_arg = {
        #         'band2process': int(band + 1),
        #         'resamp_alg': 'nearest',
        #         'nodata': 0,
        #         'align_grids': True,
        #         'q': True
        #     }
        result = DESHIFTER(dir_in + '/Merge' + img, CRL.coreg_info,
                           band2process=int(band+1),
                           align_grids=True,
                           match_gsd=True,
                           nodata=0,
                           resamp_alg='cubic',
                           q=True).correct_shifts()
        modify_hyspex_band(dest_img, band, result['arr_shifted'])

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

 # If using the global correction
    # kwargs = kwargs = {
    #     'path_out': out,
    #     'r_b4match': 1,
    #     's_b4match': 1,
    #     'max_shift': 5,
    #     'nodata': (0, 0),
    #     'resamp_alg_calc': 'nearest',
    #     'q': False
    # }
    # Arosics_test_global(ref_img, tgt_img, kwargs)
    # return CRL.coreg_info


if __name__ == '__main__':
    main()
