#!/usr/bin/python

#  OCO2toGeoTIFF.py
#  Hannah Wong
#  Converts OCO-2 data from (https://co2.jpl.nasa.gov/) to GeoTIFF raster.

# pylint: disable=no-name-in-module

import getopt
import multiprocessing as mp
import os
import sys
import time
from itertools import chain
from traceback import print_exc

import geopandas as gpd
import numpy as np
import shapely
from shapely.geometry import Polygon, MultiPolygon
import rasterio
from rasterio import features
from geopy.distance import geodesic
from netCDF4 import Dataset
from osgeo import gdal, osr


# this is default, but just to make sure
shapely.speedups.enable()


def timing(fn):
    def wrap(*args, **kwargs):
        start_time = time.time()
        result = fn(*args, **kwargs)
        end_time = time.time()
        duration = (end_time - start_time) * 1000.0
        f_name = fn.__name__
        print("{} took {:.3f} ms".format(f_name, duration))

        return result
    return wrap


def peek(iterable):
    try:
        first = next(iterable)
    except StopIteration:
        return None
    return chain([first], iterable)


def filename(file_in, directory, resolution, var_code):
    # split up file path to get file name, then split by underscore
    split = os.path.split(file_in)[-1].split('_')

    # split[2] is truncated date in YYMMDD
    return os.path.join(
        directory,
        f'oco2_{var_code}_{resolution}_20{split[2]}'
    )


def row_geometry(row):
    coordinates = np.column_stack((row[1], row[0]))
    polygon = Polygon(
        lon_lat
        for lon_lat
        in coordinates
    )

    if coordinates[:, 0].max() - coordinates[:, 0].min() > 180:
        dateline = MultiPolygon(
            ((180.0, 89), (179.95, 89), (179.95, -89), (180.0, -89), (180.0, 89)),
            ((-180.0, 89), (-179.95, 89), (-179.95, -89), (-180.0, -89), (-180.0, 89))
        )
        polygon = dateline.difference(polygon)

    return polygon


# @timing
def oco2_to_shapes(file_in, var_code, epsg, bounds):
    print(f'Reading \'{file_in}\'.')
    nc = Dataset(file_in)

    # get nparray
    vertex_lat = nc.variables['vertex_latitude'][:].filled()
    vertex_lon = nc.variables['vertex_longitude'][:].filled()
    var = nc.variables[var_code][:].filled()

    min_lon, min_lat, max_lon, max_lat = bounds

    return (
        (
            row_geometry(row),
            row[2]
        )
        for row
        in zip(vertex_lat, vertex_lon, var)
        if not any((
            row[0].min() > max_lat,
            row[0].max() < min_lat,
            (
                row[1].min() > max_lon or
                row[1].max() < min_lon
            ) if row[1].max() - row[1].min() <= 180 else (
                row[1].min() % 360 > max_lon % 360 or
                row[1].max() % 360 < min_lon % 360
            )
        ))
    )


# @timing
def create_template_raster(file_name, resolution, projection, bounds):
    min_lon, min_lat, max_lon, max_lat = bounds

    equatorial = min_lat < 0 and max_lat > 0

    width = round(
        max(
            geodesic((min_lat, min_lon), (min_lat, max_lon)).meters,
            geodesic((max_lat, min_lon), (max_lat, max_lon)).meters,
            geodesic((0, min_lon), (0, max_lon)).meters if equatorial else 0
        ) / resolution
    )
    height = round(
        geodesic((min_lat, min_lon), (max_lat, min_lon)).meters / resolution
    )
    xres = (max_lon - min_lon) / float(width)
    yres = (max_lat - min_lat) / float(height)
    geotransform = (min_lon, xres, 0, max_lat, 0, -yres)

    raster = gdal\
        .GetDriverByName('GTiff')\
        .Create(
            file_name,
            width,
            height,
            1,
            gdal.GDT_Float32
        )
    raster.SetProjection(projection)
    raster.SetGeoTransform(geotransform)
    raster.GetRasterBand(1).SetNoDataValue(0)

    # write to disk
    raster.FlushCache()
    # close file
    raster = None
    del raster


# @timing
def shapes_to_raster(shapes, file_in, meta, file_out):
    if shapes is None:
        print(f'\'{file_in}\' contains no data near bounds, skipping.')
        return

    print('Generating raster.')
    with rasterio.open(file_out, 'w+', **meta) as out:
        out_arr = out.read(1)

        rasterized = features.rasterize(
            shapes=shapes,
            fill=0,
            out=out_arr,
            transform=out.transform
        )

        # write rasterized shapes to raster
        out.write_band(1, rasterized)
        print(f'Raster \'{file_out}\' generated.')


@timing
def nc4_to_raster(file_in, directory, resolution, var_code, epsg, bounds, meta):
    shapes_to_raster(
        peek(oco2_to_shapes(file_in, var_code, epsg, bounds)),
        file_in,
        meta,
        filename(file_in, directory, f'{resolution}M2', var_code) + '.tif'
    )


def raster_wrapper(*args, **kwargs):
    try:
        nc4_to_raster(*args, **kwargs)
    except KeyboardInterrupt:
        return KeyboardInterruptError()

def average_rasters(inputs, meta, file_out):
    if len(inputs) < 1:
        return

    for i, file in enumerate(inputs):
        gd_obj = gdal.Open(file)
        array = gd_obj.ReadAsArray()
        array[array == 0] = np.nan
        array = np.expand_dims(array, 2)
        if i == 0:
            allarrays = array
        else:
            allarrays = np.concatenate((allarrays, array), axis=2)

    mean = np.nanmean(allarrays, axis=2)

    with rasterio.open(file_out, 'w+', **meta) as out:
        out.write_band(1, mean)
        print(f'Raster \'{file_out}\' generated.')


def average_wrapper(*args, **kwargs):
    try:
        average_rasters(*args, **kwargs)
    except KeyboardInterrupt:
        return KeyboardInterruptError()


def usage():
    print('There are 6 options:')
    print('    -r <number> or --resolution <number> to specify the size of each pixel, in meters. Default: 1000')
    print('    -v <string> or --variable <dir> to specify a variable from the dataset to extract. Default "xco2"')
    print('    -e <EPSG code> or --epsg <EPSG code> to specify an output projection. Default: 4326')
    print('    -b <min_lon, min_lat, max_lon, max_lat> or --bounds min_lon, min_lat, max_lon, max_lat> to specify an raster bounds. Default: "-45,45,-180,180"')
    print('    -o <dir> or --output <dir> to specify a path to an output directory. Default: ""')
    print('    -f <filename> or --file <filename> to specify an input file. Omit to process all files. Default: ""')
    print('    -a or --average to generate yearly and monthly averages. Default: False')
    print('    -h to display this message ')


def main(argv):
    directory = ''
    resolution = 1000
    var_code = 'xco2'
    epsg = 4326
    bounds = '-45, -180, 45, 180'
    file_input = None
    average = False

    try:
        opts, _args = getopt.getopt(argv, 'har:v:e:b:o:f:', [
            'help',
            'average',
            'output=',
            'resolution=',
            'variable=',
            'epsg=',
            'bounds=',
            'file=',
        ])
    except:
        print("Exception getting arguments.")
        usage()

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('-a', '--average'):
            average = True
        elif opt in ('-o', '--output'):
            directory = arg.strip(' "\'')
            print(f'Output directory is {directory}')
        elif opt in ('-r', '--resolution'):
            resolution = arg
            print(f'Output resolution is {resolution}')
        elif opt in ('-v', '--variable'):
            var_code = arg
            print(f'Variable to generate raster on is {var_code}')
        elif opt in ('-e', '--epsg'):
            try:
                epsg = int(arg)
                print(f'Output projection is EPSG:{epsg}')
            except ValueError:
                print(f'Option {opt} expects an integer value.')
                os._exit(0)
        elif opt in ('-b', '--bounds'):
            try:
                [float(x.strip()) for x in arg.strip(' "\'').split(',')]
                bounds = arg.strip(' "\'')
                print(f'Raster bounds are set to {bounds}')
            except:
                print(f'Option {opt} expects a comma separated list.')
                os._exit(0)
        elif opt in ('-f', '--file'):
            file_input = arg
            print(f'Input file is {file_input}')
        else:
            pass

    spatial_reference = osr.SpatialReference()
    spatial_reference.ImportFromEPSG(epsg)
    projection = spatial_reference.ExportToWkt()

    bounds = [float(x.strip()) for x in bounds.split(',')]

    nc4_files = [
        f
        for f
        in os.listdir(os.getcwd())
        if f.startswith('oco2') and f.endswith('.nc4')
    ] if not file_input else [file_input]

    file_out = filename(nc4_files[0], directory, resolution, var_code)
    template_path = f'{"_".join(file_out.split("_")[:3])}_template.tif'

    if not os.path.isfile(template_path):
        print('Template raster does not exist, creating now.')
        create_template_raster(
            template_path, resolution, projection, bounds)
        print(f'Template \'{template_path}\' created.')

    with rasterio.open(template_path) as template:
        meta = template.meta.copy()
        meta.update(compress='lzw')

    fn_args = (
        (f, directory, resolution, var_code, epsg, bounds, meta)
        for f
        in nc4_files
    )

    if file_input:
        raster_wrapper(*next(fn_args))
    else:
        with mp.Pool(4) as pool:
            try:
                pool.starmap(raster_wrapper, fn_args)
            except KeyboardInterrupt:
                pool.terminate()
            except Exception as e:
                print_exc()
                pool.terminate()

        if average:
            rasters = [
                os.path.join(directory, f)
                for f
                in os.listdir('1000m Rasters')
                if (
                    f.split('.')[-1] == 'tif' and
                    'template' not in f and
                    'average' not in f
                )
            ]
            oldest = min(rasters)
            latest = max(rasters)

            to_year = lambda x: int(x.split('_')[-1][:4])
            to_month = lambda x: int(x.split('_')[-1][4:6])
            average_filename_base = os.path.join(
                directory,
                f'oco2_average_{var_code}_{resolution}M2_'
            )

            year_args = [
                (
                    [
                        f
                        for f
                        in rasters
                        if to_year(f) == year
                    ],
                    meta,
                    average_filename_base + f'{year}.tif'
                )
                for year
                in range(to_year(oldest), to_year(latest) + 1)
            ]
            
            for args in year_args:
                average_rasters(*args)

            month_args = [
                [
                    (
                        [
                            f
                            for f
                            in rasters
                            if to_year(f) == year and to_month(f) == month
                        ],
                        meta,
                        average_filename_base + f'{year}{month:02}.tif'
                    )
                    for month
                    in range(1, 12)
                ]
                for year
                in range(to_year(oldest), to_year(latest) + 1)
            ]
            month_args = [item for subl in month_args for item in subl]

            with mp.Pool(4) as pool:
                try:
                    pool.starmap(average_wrapper, month_args)
                except KeyboardInterrupt:
                    pool.terminate()
                except Exception as e:
                    print_exc()
                    pool.terminate()


if __name__ == '__main__':
    main(sys.argv[1:])
