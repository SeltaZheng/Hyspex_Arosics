# python apply_geocorr_json2.py -t f140613t01p00r07_rfl_v1b_img_subset.tif -o ./ -g f140613t01p00r07_rfl_v1b_img_subset_GCPs.json --bound 540000 3743360 550000 3748000

import os, sys, argparse, json
import numpy as np

import gdal, gdalconst, osr

import warnings

warnings.filterwarnings("ignore")


def load_json_gcp_dict(json_file):
    with open(json_file, "r") as outfile:
        coreg_info = json.load(outfile)

        for item, gcp_info in enumerate(coreg_info['GCPList']):
            x, y, z, pixel, line = gcp_info
            coreg_info['GCPList'][item] = gdal.GCP(x, y, z, pixel, line)

        print("{} GCP points loaded.".format(len(coreg_info['GCPList'])))
        return coreg_info


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', type=int, required=False, help="Band number in the image to warp (index start at 1)",
                        nargs='+')
    parser.add_argument('-t', type=str, required=True, help="Target image")
    parser.add_argument('-o', type=str, required=True, help="Output directory")
    parser.add_argument('-g', type=str, required=True, help="GCP JSON file")

    parser.add_argument('--nodata', type=float, default=-9999, required=False, help="Nodata value in ouput image")
    parser.add_argument('-p', type=float, required=False, help="Output pixel size (same unit with input image)")
    parser.add_argument('--bound', type=float, nargs='+', required=False, help="Order: MinX, MinY, MaxX, MaxY")
    parser.add_argument('--ext', type=str, required=False, default="_geocorr_rot",
                        help="Output suffix (default '_geocorr_rot')")

    args = parser.parse_args()

    tgt_img = args.t
    out_dir = args.o
    gcp_file = args.g
    out_nodata = args.nodata
    output_suffix = args.ext
    # print(args.b)
    # outputBounds = None #[543900,3743360,  550000, 3748000]

    if args.b is not None:
        band2process = args.b  # int(args.b[0])+1
    else:
        band2process = None

    if args.bound is not None:
        if len(args.bound) >= 4:
            outputBounds = args.bound
        else:
            outputBounds = None
    else:
        outputBounds = None

    print(args.bound)
    # quit()

    ds = gdal.Open(tgt_img, gdal.GA_Update)
    gt = ds.GetGeoTransform()

    if args.p is not None:
        out_pixel_size = args.p
    else:
        out_pixel_size = np.sqrt(gt[1] ** 2 + gt[2] ** 2)

    proj = osr.SpatialReference(wkt=ds.GetProjection())

    coreg_info = load_json_gcp_dict(gcp_file)
    GCPs = coreg_info['GCPList']

    destname = "/vsimem/translate.tif"

    translate_option = gdal.TranslateOptions(GCPs=GCPs, bandList=band2process)
    # translate_ds = gdal.Translate(destname,ds)
    translate_ds = gdal.Translate(destname, ds, options=translate_option)
    warp_option = gdal.WarpOptions(outputBounds=outputBounds,
                                   format="ENVI",
                                   dstNodata=out_nodata,
                                   dstSRS=proj,
                                   xRes=out_pixel_size,
                                   yRes=out_pixel_size,
                                   multithread=True
                                   )
    out_img = out_dir + '/' + os.path.basename(tgt_img).split('.tif')[0] + output_suffix
    warp_ds = gdal.Warp(out_img, translate_ds, options=warp_option)

    for i_band in range(1, warp_ds.RasterCount + 1):
        warp_ds.GetRasterBand(i_band).SetNoDataValue(out_nodata)

    gdal.GetDriverByName('GTiff').Delete(destname)
    warp_ds = None
    ds = None


if __name__ == '__main__':
    main()
