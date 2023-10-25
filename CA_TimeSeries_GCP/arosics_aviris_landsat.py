# Usage:
# python arosics_test.py -h
# Example:
# python arosics_test.py  -r ref_f170607t01p00r16_rfl_v1g_img_topo_brdf_b29.tif -t f170607t01p00r16_rfl_v1g_img_topo_brdf_b29.tif -o ./  --GCP_Only

import os, sys, argparse, json, glob
import numpy as np
from numpy.linalg import inv
from osgeo import gdal, gdalconst
from arosics.DeShifter import DESHIFTER
from arosics import COREG_LOCAL

import warnings

warnings.filterwarnings("ignore")


def own_copy(coreg_info):
    new_dict = {}
    new_dict['GCPList'] = [
        gdal.GCP(gcp_swig_obj.GCPX, gcp_swig_obj.GCPY, gcp_swig_obj.GCPZ, gcp_swig_obj.GCPPixel, gcp_swig_obj.GCPLine)
        for gcp_swig_obj in coreg_info['GCPList']]
    new_dict['mean_shifts_px'] = {'x': coreg_info['mean_shifts_px']['x'], 'y': coreg_info['mean_shifts_px']['y']}
    new_dict['mean_shifts_map'] = {'x': coreg_info['mean_shifts_map']['x'], 'y': coreg_info['mean_shifts_map']['y']}
    new_dict['updated map info means'] = [x for x in coreg_info['updated map info means']]
    new_dict['original map info'] = [x for x in coreg_info['original map info']]
    new_dict['reference projection'] = coreg_info['reference projection']
    new_dict['success'] = coreg_info['success']
    new_dict['reference extent'] = {'cols': coreg_info['reference extent']['cols'],
                                    'rows': coreg_info['reference extent']['rows']}
    new_dict['reference grid'] = [
        [coreg_info['reference grid'][0][0], coreg_info['reference grid'][0][1]],
        [coreg_info['reference grid'][1][0], coreg_info['reference grid'][1][1]]
    ]
    new_dict['reference geotransform'] = [x for x in coreg_info['reference geotransform']]
    return new_dict


def update_gcp_rotation(coreg_dict, inv_rotation_matrix, geotransform, rotation_matrix_warp, save_gcp_flag):
    coreg_dict['updated map info means'][3] = geotransform[0] + coreg_dict['mean_shifts_map']['x']
    coreg_dict['updated map info means'][4] = geotransform[3] + coreg_dict['mean_shifts_map']['y']

    for item, gcp_swig_obj in enumerate(coreg_dict['GCPList']):
        x, y, z, pixel, line = gcp_swig_obj.GCPX, gcp_swig_obj.GCPY, gcp_swig_obj.GCPZ, gcp_swig_obj.GCPPixel, gcp_swig_obj.GCPLine

        new_pixel, new_line, _ = (inv_rotation_matrix @ rotation_matrix_warp) @ np.array([pixel, line, 1])
        if save_gcp_flag:
            # in order to save the dictionary in json, gdal.GCP object is converted to list
            coreg_dict['GCPList'][item] = [gcp_swig_obj.GCPX, gcp_swig_obj.GCPY, gcp_swig_obj.GCPZ, new_pixel, new_line]
        else:
            # keep gdal.GCP object
            coreg_dict['GCPList'][item] = gdal.GCP(gcp_swig_obj.GCPX, gcp_swig_obj.GCPY, gcp_swig_obj.GCPZ, new_pixel,
                                                   new_line)


def main():  # argv

    # set up directories
    yrs = ['2013']#['2013', '2014', '2015', '2016', '2017', '2018']
    boxes = ['f130626']#['f130612', 'f140603', 'f150601', 'f150602', 'f160621', 'f170607', 'f180622']
    box_name = 'Yosemite_NEON'
    dir_landsat = r'Y:\CA_timeseries\LandSAT'
    dir_avs = r'Z:\townsenduser-rw\CA_project\Raw'
    dir_out = r'G:\My Drive\Projects_ongoing\HyspIRI_Validation\TraitValidation\CA_TimeSeries\GCPs'
    gcp_save_flag = True
    nodata_tgt = -9999

    # two for loops: 1. each box; 2. each flight line

    for i, box in enumerate(boxes):
        yr = yrs[i]
        ref_img = glob.glob(f'{dir_landsat}/*{yr}*_sr_brdf.tif')[0]
        # flightline list for the box:
        dir_l = f'{dir_avs}/{yr}/{box_name}/{box}/L2'
        flights = sorted(glob.glob(f'{dir_l}/{box}*_img'))

        # only use the green band for ref:
        ref_ds = gdal.Open(ref_img)
        translate_option = gdal.TranslateOptions(bandList=[2], outputType=gdalconst.GDT_Float32)

        # ref_filename = '/vsimem/ref.tif'
        ref_filename = './ref.tif'
        gdal.Translate(ref_filename, ref_ds, options=translate_option)


        # loop 2: each flightline in the box:
        for j, flight in enumerate(flights):

            tgt_img = flight
            # nodata_tgt = args.nodata

            print(ref_img, tgt_img, dir_out, gcp_save_flag)

            # warp_img = '/vsimem/warp_tmp.tif'  # out_dir +'/'+os.path.basename(tgt_img).split('.tif')[0]+'_auto_warp.tif'
            warp_img = './warp_tmp.tif'
            # extract one band from target img, for AVIRIS-Classic, using band 21, wvl: 560nm
            tgt_ds = gdal.Open(tgt_img)
            gt = tgt_ds.GetGeoTransform()

            translate_option = gdal.TranslateOptions(bandList=[21], outputType=gdalconst.GDT_Float32)
            # translate_option = gdal.TranslateOptions(bandList=[3],outputType=gdalconst.GDT_Float32) #[3] gdalconst.GDT_Int16

            dst_filename = '/vsimem/translate.tif'

            gdal.Translate(dst_filename, tgt_ds, options=translate_option)

            gdal.Warp(warp_img, dst_filename)
            gdal.GetDriverByName('GTiff').Delete(dst_filename)



            kwargs = {
                'grid_res': 80,  # 400, # 350
                'window_size': (100, 100),  # (64, 64), #256, #150
                'path_out': None,  # out_img,
                # 'r_b4match': 187,
                # 's_b4match': 186,
                'match_gsd': False,
                'max_iter': 16,
                'max_shift': 205,
                'nodata': (0, nodata_tgt),  # (0, -0.9998),
                'CPUs': 8,
                # 'tieP_filter_level': 1, # avoid filter 2 and 3, fix the SSIM filter error
                # 'min_reliability': 25,
                'resamp_alg_calc': 'nearest',
                'resamp_alg_deshift': 'nearest',
                'align_grids': False,
                'q': True
            }

            # CRL = COREG_LOCAL(ref_img, warp_img, **kwargs)
            CRL = COREG_LOCAL(ref_filename, warp_img, **kwargs)

            try:
                CRL.calculate_spatial_shifts()
            except:
                print('{} is not processed properly'.format(tgt_img))
                # error_ls.append(tgt_img)
                # quit()
                continue
            new_coreg_info = own_copy(CRL.coreg_info)

            rot_mat = np.array([[gt[1], gt[2], gt[0]], [gt[4], gt[5], gt[3]], [0, 0, 1]])
            gt_warp = gdal.Open(warp_img).GetGeoTransform()
            warp_rot_mat = np.array([[gt_warp[1], gt_warp[2], gt_warp[0]], [gt_warp[4], gt_warp[5], gt_warp[3]], [0, 0, 1]])

            inv_rot_mat = inv(rot_mat)
            update_gcp_rotation(new_coreg_info, inv_rot_mat, gt, warp_rot_mat, gcp_save_flag)

            gdal.GetDriverByName('GTiff').Delete(warp_img)

            print("Totally {} GCPs are found.".format(len(new_coreg_info['GCPList'])))

            if not gcp_save_flag:

                out_img = dir_out + '/' + os.path.basename(tgt_img).split('.tif')[0] + '_geocorr_rot.bin'
                out_pixel_size = np.sqrt(gt[1] ** 2 + gt[2] ** 2)
                shift_arg = {
                    'band2process': None,  # int(band + 1),
                    'resamp_alg': 'nearest',
                    'CPUs': 8,
                    'out_gsd': [out_pixel_size, out_pixel_size],
                    'nodata': nodata_tgt,
                    'path_out': out_img,
                    'align_grids': False,
                    'q': True
                }

                DESHIFTER(tgt_img, new_coreg_info, **shift_arg).correct_shifts()
            else:
                # save GCPs and new_coreg_info
                out_json = dir_out + '/' + os.path.basename(tgt_img).split('.tif')[0] + '_GCPs.json'
                with open(out_json, "w") as outfile:
                    json.dump(new_coreg_info, outfile)
                    print("JSON file {} for GCP is saved.".format(out_json))

            del CRL

        gdal.GetDriverByName('GTiff').Delete(ref_filename)


if __name__ == '__main__':
    main()