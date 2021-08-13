"""
This script is used to georectify the mismatched SWIR bands usually due to poor boresight to VNIR bands
"""
import os, glob, shutil
import numpy as np
from ENVI import read_envi_header, write_envi_header
from arosics.DeShifter import DESHIFTER
from arosics import COREG_LOCAL


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


# # ---------- main -------------------
# dir_in = r'Z:\townsenduser-rw\HyspexPro\Output\Cheesehead_V3\CHEESEHEAD_20190629\Merge'
# img = '/CHEESEHEAD_20190629_03_Refl'
# header = read_envi_header(dir_in + img + '.hdr')
# band_vnir, band_swir, idx_vnir, idx_swir = band_select(header)
# # hyspex_extract(dir_in + img, [band_vnir], dir_in + img + '_refband')
# hyspex_extract(dir_in + img, [186], dir_in + img + '_refband')
#
# # src_img = dir_in + img
# # out_img = dir_in + img + '_swir'
# hyspex_extract(dir_in + img, [185], dir_in + img + '_swir_test')

##-------------- main function ------------------------------
def main():
    # parameters for the geo rectification:
    kwargs = {
        'grid_res': 400, #450, #400, # 350
        'window_size': (256, 256), #280, #256, #150
        # 'path_out': out,
        'r_b4match': 187,
        's_b4match': 186,
        'match_gsd': True,
        'max_iter': 8,
        'max_shift': 10,
        'nodata': (0, 0),
        'tieP_filter_level': 3, # avoid filter 2 and 3, fix the SSIM filter error
        # 'min_reliability': 25,
        'resamp_alg_calc': 'nearest',
        'q': False
    }
    # dir_in = r'Z:\townsenduser-rw\HyspexPro\Output\Cheesehead_V3\CHEESEHEAD_20190626\CHEESEHEAD_20190626_02_quicklooks'
    dir_in = r'Z:\townsenduser-rw\HyspexPro\Output\Cheesehead_V3\CHEESEHEAD_20190830'
    # create deshift folder for deshifted Refl data:
    dest_dir = dir_in + '/Deshift'
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
    imgs = glob.glob(dir_in + '/Merge/*_Refl')
    imgs.sort()
    # get swir bands:
    header = read_envi_header(imgs[0] + '.hdr')
    band_vnir, band_swir, idx_vnir, idx_swir = band_select(header)
    error_ls = []

    # for left over imgs:
    imgs = [imgs[x] for x in [12]]
    for img in imgs:
        # # copy file in /Merge to Deshift:
        shutil.copy(img, dest_dir)
        # file header:
        shutil.copy(img + '.hdr', dest_dir)

        # rename files:
        img_base = os.path.basename(img)
        dest_file = os.path.join(dest_dir, img_base)
        dest_img = dest_file + '_deshift'

        #
        os.rename(dest_file, dest_img)
        dest_hdr = os.path.join(dest_dir, img_base + '.hdr')
        os.rename(dest_hdr, dest_file + '_deshift.hdr')
        del dest_file, dest_hdr

        # calculate the shift info
        CRL = COREG_LOCAL(img, img, **kwargs)
        try:
            CRL.calculate_spatial_shifts()
        except:
            print('{} is not processed properly'.format(img))
            error_ls.append(img)
            continue

        # apply the shift to all swir bands:
        for band in idx_swir:
            print('Now shifting band {}'.format(band))
            shift_arg = {
                'band2process': int(band + 1),
                'resamp_alg': 'cubic',
                'nodata': 0,
                'align_grids': True,
                'q': True
            }
            result = DESHIFTER(img, CRL.coreg_info, **shift_arg).correct_shifts()
            modify_hyspex_band(dest_img, band, result['arr_shifted'])
            del result


        del CRL

    # If files are not processed properly:
    if len(error_ls) > 0:
        with open(dest_dir + '/error_files.txt', 'w') as f:
            for item in error_ls:
                f.write("%s\n" % item)


if __name__ == '__main__':
    main()
