#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 11:10:04 2019

@author: Yongzhen Fan
"""

from pathlib import Path
import os
from os.path import isdir
from os import makedirs
import sys
import numpy as np
import h5py
from netCDF4 import Dataset
import time
from global_land_mask import globe

from src.L1B import L1B
from src.auxUtils import AUXData
from src.sensorinfo import sensorinfo
from src.ancillary import ANCILLARY
from src.rayleigh import Rayleigh
from src.cloudmask import Cloudmask
#from glint import get_glint_coeff
from src.mlnn import MLNN
from src.ocparam import CHL, TSM, CDOM
from src.config import Config
import warnings

DEBUG_ANCILLARY_ARRAYS = False  # Set to True if ancillary arrays (e.g. ozone concentration, no2_concentration, etc.) need to be populated (e.g. for easy viewing in the debugger)

def print_L2_write_warning():
    GAS_BITS = {'oz': 1, 'co2': 2, 'no2': 4, 'h2o': 8, 'co': 16, 'ch4': 32, 'n2o': 64, 'o2': 128}
    
    for prod in l2_prod:
        for prefix in ('tg_sol_', 'tg_sen_'):
            if prod.startswith(prefix):
                gas = prod[len(prefix):]
                if gas in GAS_BITS and not (sinfo.gasid & GAS_BITS[gas]):
                    print(f"\033[33mWarning: '{prod}' was requested in l2_prod but {sinfo.sensor} "
                        f"does not apply a '{gas}' correction (gasid={sinfo.gasid}). "
                        f"This product will not be written to the L2 file.\033[0m")

#######  Read Input parameters  #####################
input_param = {}
OCSMART_script_dir = str(Path(sys.argv[0]).resolve().parent)
with open(OCSMART_script_dir + "/OCSMART_Input.txt") as inputfile:
    for line in inputfile:
        if "=" in line:
            name, var = line.split("=")
            if "path" in name:
                var = var.replace(".", OCSMART_script_dir)
            input_param[name.strip()] = var.strip()
####################################################
        
#######  directory setup  ##########################
if 'l1b_path' in input_param.keys():
    L1_path = input_param['l1b_path']
else:
    print('\033[1;31;47mError: no L1B data directory, please set L1B_path in OCSMART_input.txt','\033[m')
    print('Quit OC-SMART ...')
    time.sleep(5)
    sys.exit()
if 'l2_path' in input_param.keys():
    L2_path = input_param['l2_path']
else:
    print('Warning: L2 data directory not specified, using default L2 data directory:' + OCSMART_script_dir +  '/L2/')
    L2_path = OCSMART_script_dir + '/L2/'
####################################################
    
#######  parameter setup  ##########################
if 'solz_limit' in input_param.keys():
    solz_limit = float(input_param['solz_limit'])
else:    
    solz_limit = 70.0
if 'senz_limit' in input_param.keys():
    senz_limit = float(input_param['senz_limit'])
else:    
    senz_limit = 70.0
if 'l2_prod' in input_param.keys():
    l2_prod=input_param['l2_prod'].split(',')
    l2_prod = [x.strip() for x in l2_prod]
else:
    print('List of level-2 products not provided, using default list: rrs, chl')
    l2_prod = ['rrs','chl']
#######  Floating Point Precision Options  ##########################
if 'floating_point_datatype' in input_param.keys():
    default_floating_point_type = 'float32'
    datatype = input_param.get('floating_point_datatype', default_floating_point_type)
    if datatype.lower() not in ['float32', 'float64']:
        print(f"Warning: Supported floating point datatypes are float32 and float64. The floating_point_datatype value in OCSMART_Input.txt was {datatype}, setting to default type {default_floating_point_type}.")
        datatype = default_floating_point_type
    Config.datatype = datatype

    if datatype in ['float64']:
        # If the user intentionally chose float64, don't clutter up the console output with unnecessary warnings.
        warnings.filterwarnings("ignore", module="gas_corrections_lib")
#######  L2 file write settings  ##########################
gzip_compression_opt = int(input_param.get('gzip_compression_opt', 4))
if gzip_compression_opt < 1 or gzip_compression_opt > 9:
    gzip_compression_opt = min(9, max(1, gzip_compression_opt))
    print(f"gzip_compression_opt must be between 1 and 9, setting gzip_compression_opt to {gzip_compression_opt}")

shuffle_opt = (str(input_param.get('shuffle_opt', True)).lower() == "true")
compression_algorithm = input_param.get('compression_algorithm', "gzip")

compression_kwargs = dict(
    compression=compression_algorithm,
    shuffle=shuffle_opt
    )
if compression_algorithm == "gzip":
    compression_kwargs["compression_opts"] = gzip_compression_opt

water_subpixl_limit = 0.95 # for land/water mask
#glint_max = 0.05 (not needed)
####################################################

#######  subimage setup  ##############################################
# 3 options to define the subimage in the input file
# Option A. (north, south, east, west)
# Option B. (latitude_center, longitude_center, box_width, box_height)
# Option C. (start_line, end_line, start_pixel, end_pixel)
#######################################################################

mission = \
{ 
    'MODISA'    : 'Aqua MODIS',
    'MODIST'    : 'Terra MODIS',
    'SeaWiFS'   : 'SeaWiFS',
    'SNPP'      : 'Suomi-NPP VIIRS',
    'JPSS1'     : 'NOAA-20 VIIRS',
    'JPSS2'     : 'NOAA-21 VIIRS',
    'GOCI'      : 'COMS GOCI',
    'SGLI'      : 'GCOM-C SGLI',
    'L08'       : 'Landsat-8 OLI',
    'L09'       : 'Landsat-9 OLI',
    'S3A'       : 'Sentinel-3A OLCI',
    'S3B'       : 'Sentinel-3B OLCI',
    'S2A'       : 'Sentinel-2A MSI',
    'S2B'       : 'Sentinel-2B MSI',
    'EPIC'      : 'DSCOVR EPIC',
    'MERSI2'    : 'FengYun-3D MERSI-II',
    'HICO'      : 'ISS HICO',
    'PACE'      : 'PACE OLI',
    'HYPSO'     : 'HYPSO HSI'
}

print('Start OC-SMART ... \n')
if isdir(L1_path):
    print('Input level-1 data directory : {}'.format(L1_path))
else:
    print('Input level-1 data directory : {} does not exist, please check input file ...'.format(L1_path)) 
    print('Quit OC-SMART ...')
    sys.exit()
if isdir(L2_path):
    print('Output level-2 data directory : {}'.format(L2_path))
else:
    print('Level-2 data directory: {} does not exist, creating level-2 data directory ... '.format(L2_path))
    makedirs(L2_path)
#print('Output level 2 data format : {} \n'.format(L2_format))

# list all files in the L1 directory
#L1files = [f for f in os.listdir(L1_path) if isfile(join(L1_path, f))]
L1files = sorted([f for f in os.listdir(L1_path)])
nfiles = len(L1files)
print('{} files found in the level-1 directory. \n'.format(nfiles))

print('Level-2 products: {} \n'.format(', '.join(l2_prod)))
print(f"Compression Algorithm: {compression_algorithm}")
if compression_algorithm == "gzip":
    print(f"GZIP Compression Option: {gzip_compression_opt}")
print(f"Compression Shuffle Option: {shuffle_opt}\n")

# read auxilary data (land/water mask)
aux=AUXData()

# start processing all files
for ifile in np.arange(nfiles):
    fname = L1files[ifile]
    t_start = time.time()    
    print('Processing file {}  {}'.format(ifile+1, fname))    
    
    # get sensor information
    sinfo = sensorinfo(L1_path + fname)    

    if 'geo_path' in input_param.keys():    
        GEO_path = input_param['geo_path']
    else:
        if sinfo.sat in ['MODISA','MODIST','MERSI2']:
            print('\033[1;31;47mError: Geo file path must be provided for {}, please set GEO_path in OCSMART_input.txt'.format(mission[sinfo.sat]),'\033[m')
            print('Quit OC-SMART ...')
            time.sleep(5)
            sys.exit()
        else:                
            GEO_path = ''
    if sinfo.sensor_status == 0:
        print_L2_write_warning()
                    
        print('Sensor :',mission[sinfo.sat]) 
        
        # read level1B data
        l1b=L1B(sensorinfo=sinfo,L1Bname=L1_path+fname,GEOpath=GEO_path)
        if l1b.sensor in ['OLI','OLI2']:
            OLI_timestamp = l1b.l8metadata['LANDSAT_METADATA_FILE']['IMAGE_ATTRIBUTES']['SCENE_CENTER_TIME']
            OLI_timestamp = OLI_timestamp.replace(':', '')
            OLI_timestamp = OLI_timestamp.split('.')[0]

            fname = fname.split('.')[0] + '_' + OLI_timestamp + '.' + fname.split('.')[1]
        l1b.readgeo()
        if l1b.geoloc_status == 0:
            if 'north' in input_param.keys() and 'south' in input_param.keys() and 'east' in input_param.keys() and 'west' in input_param.keys():
                north = float(input_param['north'])
                south = float(input_param['south'])
                east = float(input_param['east'])
                west = float(input_param['west'])                
                l1b.latlon2linepixl(north=north, south=south, east=east, west=west)
            elif 'latitude_center' in input_param.keys() and 'longitude_center' in input_param.keys() and 'box_width' in input_param.keys() and 'box_height' in input_param.keys():
                lat_center = float(input_param['latitude_center'])
                lon_center = float(input_param['longitude_center'])
                box_width = float(input_param['box_width'])
                box_height = float(input_param['box_height'])
                l1b.latlon2linepixl(lat_center=lat_center, lon_center=lon_center, box_width=box_width, box_height=box_height)
            elif 'start_line' in input_param.keys() and 'end_line' in input_param.keys() and 'start_pixel' in input_param.keys() and 'end_pixel' in input_param.keys():
                sline = int(input_param['start_line'])
                eline = int(input_param['end_line'])
                spixl = int(input_param['start_pixel'])
                epixl = int(input_param['end_pixel'])    
                l1b.latlon2linepixl(start_line=sline, end_line=eline, start_pixel=spixl, end_pixel=epixl) 
            else:
                l1b.sline = 0
                l1b.spixl = 0
                l1b.eline = l1b.imagedim[0]
                l1b.epixl = l1b.imagedim[1]

            if l1b.lp_status ==0:
                l1b.readl1b()

                #######  block processing setup  ###################################################
                # block processing can be used to save memory, but it takes longer processing time.
                # if you have enough memory, set block_size = -1 to turn off the block processing
                block_size=-1
                if 'block_size' in input_param.keys():
                    block_size = int(input_param['block_size'])
                
                # check ancillary files and download from NASA if needed
                anc=ANCILLARY(l1b)
                anc.download()
                anc.read_no2()
                anc.read_ozone()
                anc.read_met()
                anc.read_gas_transmittance_auxdata(sinfo)
                
                #read Rayleigh table
                ray = Rayleigh(info=sinfo)
                
                # get dimension
                l1b_dim = l1b.reflectance.shape
                
                #initialize cloud mask
                cm = Cloudmask(sensorinfo=sinfo)
                
                # initialize the Multilayer Neural Networks
                mlnn = MLNN(sensorinfo=sinfo)
                
                # initialize CHLa, CDOM and TSM algorithm
                chl = CHL(sensorinfo=sinfo)
                tsm = TSM(sensorinfo=sinfo)
                cdom = CDOM(sensorinfo=sinfo)
                
                # initialize all needed matrices 
                # niopbands = np.sum(sinfo.band < 700)       
                
                # mask invalid radiances data
                mask_valid = np.sum(l1b.reflectance <= 0.0, 2) == 0
                l2_mask = np.zeros([l1b_dim[0], l1b_dim[1]], dtype='int16')
                l2_mask[~mask_valid] = 1
                
                # mask high solar and sensor zenith angels
                if sinfo.sensor == 'EPIC':
                    senz_limit = 65.0
                mask_solz=l1b.solz < solz_limit
#                mask_nsolz=l1b.solz >= solz_limit
                mask_senz=l1b.senz < senz_limit
#                mask_nsenz=l1b.senz >= senz_limit
                mask_valid_geo = (mask_valid & (mask_solz & mask_senz)) 
#                l2_mask[mask_valid & (mask_nsolz | mask_nsenz)] = 4
                l2_mask[mask_valid & ~(mask_solz & mask_senz)] = 4
                
                #compute finest resolution for land/water mask
                midimg_x=int(l1b.latitude.shape[0]/2)
                midimg_y=int(l1b.latitude.shape[1]/2)
                lat_delta=np.array([np.abs(l1b.latitude[midimg_x,midimg_y]-l1b.latitude[midimg_x+1,midimg_y]),\
                                           np.abs(l1b.latitude[midimg_x,midimg_y]-l1b.latitude[midimg_x-1,midimg_y]),\
                                           np.abs(l1b.latitude[midimg_x,midimg_y]-l1b.latitude[midimg_x,midimg_y+1]),\
                                           np.abs(l1b.latitude[midimg_x,midimg_y]-l1b.latitude[midimg_x,midimg_y-1])])
                lon_delta=np.array([np.abs(l1b.longitude[midimg_x,midimg_y]-l1b.longitude[midimg_x+1,midimg_y]),\
                                           np.abs(l1b.longitude[midimg_x,midimg_y]-l1b.longitude[midimg_x-1,midimg_y]),\
                                           np.abs(l1b.longitude[midimg_x,midimg_y]-l1b.longitude[midimg_x,midimg_y+1]),\
                                           np.abs(l1b.longitude[midimg_x,midimg_y]-l1b.longitude[midimg_x,midimg_y-1])])
                
                l1b_dxdy=np.amax([np.mean(lat_delta[np.where(lat_delta<1)[0]]),\
                                  np.mean(lon_delta[np.where(lon_delta<1)[0]])])
                
                # compute land/water mask
                if sinfo.sensor in ['EPIC', 'VIIRS', 'MODIS-Aqua', 'MODIS-Terra', 'SeaWiFS','MERSI2', 'OCI']:
                    water_portion = np.zeros([l1b_dim[0], l1b_dim[1]], dtype=Config.datatype)
                    water_portion[mask_valid_geo]=aux.maskland(l1b.latitude[mask_valid_geo], l1b.longitude[mask_valid_geo], l1b_dxdy)
                    mask_water = water_portion > water_subpixl_limit
                    mask_nwater = water_portion <= water_subpixl_limit
                elif sinfo.sensor in ['SGLI', 'OLCI', 'OLI', 'OLI2', 'S2A', 'S2B', 'GOCI','HICO']:
                    mask_water = l1b.landmask == 0
                    mask_nwater = ~mask_water    
                elif sinfo.sensor in ['HYPSO_HSI']:
                    mask_water = globe.is_ocean(l1b.latitude, l1b.longitude)
                    mask_nwater = ~mask_water
                
                # some images are too large, like OLI, OLI2, GOCI, set block processing anyways                 
                if block_size < 0:
                    if(l1b_dim[0]>3250 or l1b_dim[1]>3250):                                               
                        block_size=3250                    
                    
                #mask all land pixels
                mask_valid_geo_water = mask_valid_geo & mask_water
                l2_mask[mask_valid_geo & mask_nwater] = 16

                training_bands = sinfo.training_bands

                tg_sol = np.ones(l1b_dim,dtype=Config.datatype)
                tg_sen = np.ones(l1b_dim,dtype=Config.datatype)

                # write L2 file in H5 format
                print('Open level-2 file {} ... '.format(os.path.splitext(fname)[0] + '_L2_OCSMART.h5 for writing'))
                print('Values will be written to the L2 file as they are computed and then removed from RAM to reduce peak memory usage.')

                with h5py.File(L2_path + os.path.splitext(fname)[0] + '_L2_OCSMART.h5', 'w') as hf:
                    # compute transmittance of gases
                    if sinfo.gasid & 1 > 0: 
                        solar_zenith_oz_masked, sensor_zenith_oz_masked = anc.compute_ozone_transmittance(sinfo.koz,\
                                                                                        l1b.latitude[mask_valid_geo_water],\
                                                                                        l1b.longitude[mask_valid_geo_water], \
                                                                                        l1b.solz[mask_valid_geo_water], \
                                                                                        l1b.senz[mask_valid_geo_water])

                        tg_sol[mask_valid_geo_water, :] *= solar_zenith_oz_masked
                        tg_sen[mask_valid_geo_water, :] *= sensor_zenith_oz_masked

                        # solar_zenith_oz_masked_py, sensor_zenith_oz_masked_py = anc.trans_ozone(sinfo.koz,\
                        #                                                                    l1b.latitude[mask_valid_geo_water],\
                        #                                                                    l1b.longitude[mask_valid_geo_water], \
                        #                                                                    l1b.solz[mask_valid_geo_water], \
                        #                                                                    l1b.senz[mask_valid_geo_water])

                        if DEBUG_ANCILLARY_ARRAYS:
                            ozone_concentration = np.zeros([l1b_dim[0], l1b_dim[1]], dtype=Config.datatype)
                            ozone_concentration[mask_valid_geo_water] = anc.ozone_concentration

                        if 'tg_sol_oz' in l2_prod:
                            solar_zenith_oz = np.ones(l1b_dim,dtype=Config.datatype)
                            # solar_zenith_oz_py = np.ones(l1b_dim,dtype=Config.datatype)

                            solar_zenith_oz[mask_valid_geo_water, :] = solar_zenith_oz_masked

                            gw = hf.create_group('tg_sol_oz')
                            for i in np.arange(mlnn.aodnn_layers[-1]):
                                gw.create_dataset('tg_sol_oz_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = solar_zenith_oz[:,:,i], **compression_kwargs)

                            del solar_zenith_oz

                        if 'tg_sen_oz' in l2_prod:
                            sensor_zenith_oz = np.ones(l1b_dim,dtype=Config.datatype)
                            # sensor_zenith_oz_py = np.ones(l1b_dim,dtype=Config.datatype)

                            sensor_zenith_oz[mask_valid_geo_water, :] = sensor_zenith_oz_masked
                            
                            gw = hf.create_group('tg_sen_oz')
                            for i in np.arange(mlnn.aodnn_layers[-1]):
                                gw.create_dataset('tg_sen_oz_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = sensor_zenith_oz[:,:,i], **compression_kwargs)

                            del sensor_zenith_oz

                        del solar_zenith_oz_masked, sensor_zenith_oz_masked
                    
                    if sinfo.gasid & 4 > 0:
                        solar_zenith_no2_masked, sensor_zenith_no2_masked = anc.compute_no2_transmittance(sinfo.kno2, \
                                                                                            l1b.month, \
                                                                                            l1b.latitude[mask_valid_geo_water],\
                                                                                            l1b.longitude[mask_valid_geo_water], \
                                                                                            l1b.solz[mask_valid_geo_water], \
                                                                                            l1b.senz[mask_valid_geo_water]) 

                        # solar_zenith_no2_masked, sensor_zenith_no2_masked = anc.compute_no2_transmittance(sinfo.kno2, \
                        #                                                     l1b.month, \
                        #                                                     l1b.latitude[np.full(mask_valid_geo_water.shape, True)],\
                        #                                                     l1b.longitude[np.full(mask_valid_geo_water.shape, True)], \
                        #                                                     l1b.solz[np.full(mask_valid_geo_water.shape, True)], \
                        #                                                     l1b.senz[np.full(mask_valid_geo_water.shape, True)])
                        
                        # solar_zenith_no2_masked_py, sensor_zenith_no2_masked_py = anc.trans_no2(sinfo.kno2, \
                        #                                                                      l1b.month, \
                        #                                                                      l1b.latitude[mask_valid_geo_water],\
                        #                                                                      l1b.longitude[mask_valid_geo_water], \
                        #                                                                      l1b.solz[mask_valid_geo_water], \
                        #                                                                      l1b.senz[mask_valid_geo_water])

                        tg_sol[mask_valid_geo_water, :] *= solar_zenith_no2_masked
                        tg_sen[mask_valid_geo_water, :] *= sensor_zenith_no2_masked

                        if DEBUG_ANCILLARY_ARRAYS:
                            fraction_tropospheric_no2_above_200m = np.zeros([l1b_dim[0], l1b_dim[1]], dtype=Config.datatype)
                            stratospheric_no2_concentration = np.zeros([l1b_dim[0], l1b_dim[1]], dtype=Config.datatype)
                            tropospheric_no2_concentration = np.zeros([l1b_dim[0], l1b_dim[1]], dtype=Config.datatype)

                            fraction_tropospheric_no2_above_200m[mask_valid_geo_water] = anc.fraction_tropospheric_no2_above_200m
                            stratospheric_no2_concentration[mask_valid_geo_water] = anc.stratospheric_no2_concentration
                            tropospheric_no2_concentration[mask_valid_geo_water] = anc.tropospheric_no2_concentration

                        if 'tg_sol_no2' in l2_prod:
                            solar_zenith_no2 = np.ones(l1b_dim,dtype=Config.datatype)
                            # solar_zenith_no2_py = np.ones(l1b_dim,dtype=Config.datatype)

                            gw = hf.create_group('tg_sol_no2')
                            for i in np.arange(mlnn.aodnn_layers[-1]):
                                gw.create_dataset('tg_sol_no2_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = solar_zenith_no2[:,:,i], **compression_kwargs)

                            del solar_zenith_no2

                        if 'tg_sen_no2' in l2_prod:
                            sensor_zenith_no2 = np.ones(l1b_dim,dtype=Config.datatype)
                            # sensor_zenith_no2_py = np.ones(l1b_dim,dtype=Config.datatype)

                            gw = hf.create_group('tg_sen_no2')
                            for i in np.arange(mlnn.aodnn_layers[-1]):
                                gw.create_dataset('tg_sen_no2_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = sensor_zenith_no2[:,:,i], **compression_kwargs)

                            del sensor_zenith_no2

                        del solar_zenith_no2_masked, sensor_zenith_no2_masked
                    

                    if sinfo.gasid & 2 > 0:
                        solar_zenith_co2_masked, sensor_zenith_co2_masked = \
                        anc.compute_co2_transmittance(l1b.solz[mask_valid_geo_water], l1b.senz[mask_valid_geo_water], sinfo.band)

                        tg_sol[mask_valid_geo_water, :] *= solar_zenith_co2_masked
                        tg_sen[mask_valid_geo_water, :] *= sensor_zenith_co2_masked

                        if 'tg_sol_co2' in l2_prod:
                            solar_zenith_co2 = np.ones(l1b_dim, dtype=Config.datatype)

                            gw = hf.create_group('tg_sol_co2')
                            for i in np.arange(mlnn.aodnn_layers[-1]):
                                gw.create_dataset('tg_sol_co2_'+str(training_bands[i])+'nm',dtype=Config.datatype,data = solar_zenith_co2[:,:,i],**compression_kwargs)

                            del solar_zenith_co2

                        if 'tg_sen_co2' in l2_prod:
                            sensor_zenith_co2 = np.ones(l1b_dim, dtype=Config.datatype)

                            gw = hf.create_group('tg_sen_co2')
                            for i in np.arange(mlnn.aodnn_layers[-1]):
                                gw.create_dataset('tg_sen_co2_'+str(training_bands[i])+'nm',dtype=Config.datatype,data = sensor_zenith_co2[:,:,i],**compression_kwargs)

                            del sensor_zenith_co2

                        del solar_zenith_co2_masked, sensor_zenith_co2_masked

                    if sinfo.gasid & 16 > 0:
                        solar_zenith_co_masked, sensor_zenith_co_masked = \
                        anc.compute_co_transmittance(l1b.solz[mask_valid_geo_water], l1b.senz[mask_valid_geo_water], sinfo.band)

                        tg_sol[mask_valid_geo_water, :] *= solar_zenith_co_masked
                        tg_sen[mask_valid_geo_water, :] *= sensor_zenith_co_masked

                        if 'tg_sol_co' in l2_prod:
                            solar_zenith_co = np.ones(l1b_dim, dtype=Config.datatype)

                            gw = hf.create_group('tg_sol_co')
                            for i in np.arange(mlnn.aodnn_layers[-1]):
                                gw.create_dataset('tg_sol_co_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = solar_zenith_co[:,:,i], **compression_kwargs)

                            del solar_zenith_co

                        if 'tg_sen_co' in l2_prod:
                            sensor_zenith_co = np.ones(l1b_dim, dtype=Config.datatype)

                            gw = hf.create_group('tg_sen_co')
                            for i in np.arange(mlnn.aodnn_layers[-1]):
                                gw.create_dataset('tg_sen_co_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = sensor_zenith_co[:,:,i], **compression_kwargs)

                            del sensor_zenith_co

                        del solar_zenith_co_masked, sensor_zenith_co_masked

                    if sinfo.gasid & 32 > 0:
                        solar_zenith_ch4_masked, sensor_zenith_ch4_masked = \
                        anc.compute_ch4_transmittance(l1b.solz[mask_valid_geo_water], l1b.senz[mask_valid_geo_water], sinfo.band)

                        tg_sol[mask_valid_geo_water, :] *= solar_zenith_ch4_masked
                        tg_sen[mask_valid_geo_water, :] *= sensor_zenith_ch4_masked

                        if 'tg_sol_ch4' in l2_prod:
                            solar_zenith_ch4 = np.ones(l1b_dim, dtype=Config.datatype)

                            gw = hf.create_group('tg_sol_ch4')
                            for i in np.arange(mlnn.aodnn_layers[-1]):
                                gw.create_dataset('tg_sol_ch4_' + str(training_bands[i]) + 'nm',dtype=Config.datatype, data = solar_zenith_ch4[:,:,i], **compression_kwargs)

                            del solar_zenith_ch4

                        if 'tg_sen_ch4' in l2_prod:
                            sensor_zenith_ch4 = np.ones(l1b_dim, dtype=Config.datatype)
                            
                            gw = hf.create_group('tg_sen_ch4')
                            for i in np.arange(mlnn.aodnn_layers[-1]):
                                gw.create_dataset('tg_sen_ch4_' + str(training_bands[i]) + 'nm',dtype=Config.datatype, data = sensor_zenith_ch4[:,:,i], **compression_kwargs)

                            del sensor_zenith_ch4

                        del solar_zenith_ch4_masked, sensor_zenith_ch4_masked

                    if sinfo.gasid & 128 > 0:
                        solar_zenith_o2_masked, sensor_zenith_o2_masked = \
                        anc.compute_o2_transmittance(l1b.solz[mask_valid_geo_water], l1b.senz[mask_valid_geo_water], l1b.reflectance[mask_valid_geo_water], sinfo.band)

                        tg_sol[mask_valid_geo_water, :] *= solar_zenith_o2_masked
                        tg_sen[mask_valid_geo_water, :] *= sensor_zenith_o2_masked

                        if 'tg_sol_o2' in l2_prod:
                            solar_zenith_o2 = np.ones(l1b_dim, dtype=Config.datatype)

                            gw = hf.create_group('tg_sol_o2')
                            for i in np.arange(mlnn.aodnn_layers[-1]):
                                gw.create_dataset('tg_sol_o2_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = solar_zenith_o2[:,:,i], **compression_kwargs)

                            del solar_zenith_o2

                        if 'tg_sen_o2' in l2_prod:
                            sensor_zenith_o2 = np.ones(l1b_dim, dtype=Config.datatype)

                            gw = hf.create_group('tg_sen_o2')
                            for i in np.arange(mlnn.aodnn_layers[-1]):
                                gw.create_dataset('tg_sen_o2_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = sensor_zenith_o2[:,:,i], **compression_kwargs)

                            del sensor_zenith_o2

                        del solar_zenith_o2_masked, sensor_zenith_o2_masked

                    if sinfo.gasid & 64 > 0:
                        solar_zenith_n2o_masked, sensor_zenith_n2o_masked = \
                        anc.compute_n2o_transmittance(l1b.solz[mask_valid_geo_water], l1b.senz[mask_valid_geo_water], sinfo.band)

                        tg_sol[mask_valid_geo_water, :] *= solar_zenith_n2o_masked
                        tg_sen[mask_valid_geo_water, :] *= sensor_zenith_n2o_masked

                        if 'tg_sol_n2o' in l2_prod:
                            solar_zenith_n2o = np.ones(l1b_dim, dtype=Config.datatype)

                            gw = hf.create_group('tg_sol_n2o')
                            for i in np.arange(mlnn.aodnn_layers[-1]):
                                gw.create_dataset('tg_sol_n2o_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = solar_zenith_n2o[:,:,i], **compression_kwargs)

                            del solar_zenith_n2o

                        if 'tg_sen_n2o' in l2_prod:
                            sensor_zenith_n2o = np.ones(l1b_dim, dtype=Config.datatype)

                            gw = hf.create_group('tg_sen_n2o')
                            for i in np.arange(mlnn.aodnn_layers[-1]):
                                gw.create_dataset('tg_sen_n2o_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = sensor_zenith_n2o[:,:,i], **compression_kwargs)

                            del sensor_zenith_n2o

                        del solar_zenith_n2o_masked, sensor_zenith_n2o_masked

                    if sinfo.gasid & 8 > 0:
                        solar_zenith_h2o_masked, sensor_zenith_h2o_masked = \
                        anc.compute_h2o_transmittance(l1b.solz[mask_valid_geo_water], l1b.senz[mask_valid_geo_water], l1b.reflectance[mask_valid_geo_water], sinfo.band)

                        tg_sol[mask_valid_geo_water, :] *= solar_zenith_h2o_masked
                        tg_sen[mask_valid_geo_water, :] *= sensor_zenith_h2o_masked

                        if 'tg_sol_h2o' in l2_prod:
                            solar_zenith_h2o = np.ones(l1b_dim, dtype=Config.datatype)

                            gw = hf.create_group('tg_sol_h2o')
                            for i in np.arange(mlnn.aodnn_layers[-1]):
                                gw.create_dataset('tg_sol_h2o_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = solar_zenith_h2o[:,:,i], **compression_kwargs)

                            del solar_zenith_h2o

                        if 'tg_sen_h2o' in l2_prod:
                            sensor_zenith_h2o = np.ones(l1b_dim, dtype=Config.datatype)

                            gw = hf.create_group('tg_sen_h2o')
                            for i in np.arange(mlnn.aodnn_layers[-1]):
                                gw.create_dataset('tg_sen_h2o_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = sensor_zenith_h2o[:,:,i], **compression_kwargs)

                            del sensor_zenith_h2o

                        del solar_zenith_h2o_masked, sensor_zenith_h2o_masked

                    # get met data: pressure, RH, windspeed
                    l1b_pressure = np.zeros([l1b_dim[0], l1b_dim[1]], dtype=Config.datatype)
                    l1b_rh = np.zeros([l1b_dim[0], l1b_dim[1]], dtype=Config.datatype)
                    l1b_ws = np.zeros([l1b_dim[0], l1b_dim[1]], dtype=Config.datatype)
                    l1b_pressure[mask_valid_geo_water], l1b_rh[mask_valid_geo_water], l1b_ws[mask_valid_geo_water] = anc.get_metdata(l1b.latitude[mask_valid_geo_water],l1b.longitude[mask_valid_geo_water])

                    # get Rayleigh reflectance and correct for real time pressure
                    l1b_ray = np.full(l1b_dim, np.nan, dtype=Config.datatype)
                    l1b_ray[mask_valid_geo_water,:] = ray.corr_ray(l1b.solz[mask_valid_geo_water],\
                                                            l1b.senz[mask_valid_geo_water],\
                                                            l1b.relaz[mask_valid_geo_water],\
                                                            l1b_ws[mask_valid_geo_water])

                    # press_fac = np.full(l1b_dim, np.nan, dtype=Config.datatype)  # press_fac is only used as in intermediate variable no need to pre-allocate it (which costs unnecessary memory for large images)
                    
                    press_fac_masked = ray.corr_ray_press(l1b.solz[mask_valid_geo_water],\
                                                                    l1b.senz[mask_valid_geo_water],\
                                                                    l1b_pressure[mask_valid_geo_water])
                    
                    l1b_ray[mask_valid_geo_water,:] = l1b_ray[mask_valid_geo_water,:] * press_fac_masked
                    del press_fac_masked
                    
                    # Oxygen absorption correction for SeaWiFS band 7
                    if sinfo.sensor == 'SeaWiFS':
                        tg_o2 = np.ones([l1b_dim[0],l1b_dim[1]], dtype=Config.datatype)
                        tg_o2[mask_valid_geo_water] = anc.trans_o2_ray(l1b.solz[mask_valid_geo_water], l1b.senz[mask_valid_geo_water])
                        l1b_ray[mask_valid_geo_water,6] = l1b_ray[mask_valid_geo_water,6] * tg_o2[mask_valid_geo_water]


                    # get whitecaps reflectance
                    l1b_wcaps = np.full(l1b_dim, np.nan, dtype=Config.datatype)
                    l1b_wcaps[mask_valid_geo_water,:] = anc.whitecaps(sinfo.band,l1b.solz[mask_valid_geo_water],\
                                                                l1b.senz[mask_valid_geo_water], \
                                                                l1b_ws[mask_valid_geo_water], \
                                                                l1b_pressure[mask_valid_geo_water], ray.taur)
                    
                    # compute Rayleigh corrected reflectance and mask negative value, if any
                    lrc = np.full(l1b_dim, -999., dtype=Config.datatype)

                    # lrc[mask_valid_geo_water,:] = l1b.reflectance[mask_valid_geo_water,:]/solar_zenith_oz[mask_valid_geo_water,:]/sensor_zenith_oz[mask_valid_geo_water,:]\
                    #                            /solar_zenith_no2[mask_valid_geo_water,:]/sensor_zenith_no2[mask_valid_geo_water,:]-l1b_wcaps[mask_valid_geo_water,:]-l1b_ray[mask_valid_geo_water,:]
    #                lrc[mask_valid_geo_water,:] = l1b.reflectance[mask_valid_geo_water,:]/solar_zenith_oz[mask_valid_geo_water,:]/sensor_zenith_oz[mask_valid_geo_water,:]\
    #                                           /solar_zenith_no2[mask_valid_geo_water,:]/sensor_zenith_no2[mask_valid_geo_water,:]-l1b_wcaps[mask_valid_geo_water,:]


                    print("Applying Corrections...")
                    lrc[mask_valid_geo_water,:] = l1b.reflectance[mask_valid_geo_water,:]
                    
                    lrc[mask_valid_geo_water,:] = lrc[mask_valid_geo_water,:]/tg_sol[mask_valid_geo_water,:]/tg_sen[mask_valid_geo_water,:]\
                                                - l1b_wcaps[mask_valid_geo_water,:] - l1b_ray[mask_valid_geo_water,:]

                    # Delete some arrays to save memory, but only after writing them to the file if they are requested
                    if 'tg_sol' in l2_prod:
                        gw = hf.create_group('tg_sol')
                        for i in np.arange(mlnn.aodnn_layers[-1]):
                            gw.create_dataset('tg_sol_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = tg_sol[:,:,i], **compression_kwargs)
                    del tg_sol

                    if 'tg_sen' in l2_prod:
                        gw = hf.create_group('tg_sen')
                        for i in np.arange(mlnn.aodnn_layers[-1]):
                            gw.create_dataset('tg_sen_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = tg_sen[:,:,i], **compression_kwargs)
                    del tg_sen

                    if 'Lr' in l2_prod:
                        gw = hf.create_group('Lr')
                        for i in np.arange(mlnn.aodnn_layers[-1]):
                            gw.create_dataset('Lr_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = l1b_ray[:,:,i], **compression_kwargs)
                    del l1b_ray

                    if 'Lwp' in l2_prod:
                        gw = hf.create_group('Lwp')
                        for i in np.arange(mlnn.aodnn_layers[-1]):
                            gw.create_dataset('Lwp_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = l1b_wcaps[:,:,i], **compression_kwargs)
                    del l1b_wcaps

                    if 'Lt' in l2_prod:
                        gw = hf.create_group('Lt')
                        for i in np.arange(mlnn.aodnn_layers[-1]):
                            gw.create_dataset('Lt_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = l1b.reflectance[:,:,i], **compression_kwargs)
                    del l1b.reflectance

                    neg = np.sum(lrc < 0.0, axis=2)
                    mask_valid_geo_water_lrcposi = (neg == 0)
                    lrc[lrc == -999] = 999
                    neg1 = np.sum(lrc < 0.0, axis=2)
                    mask_valid_geo_water_lrcneg = neg1 > 0
                    l2_mask[mask_valid_geo_water_lrcneg] = 1024
                    
                    # compute sunglint risk (not needed)
                    # glint_coeff = np.zeros([l1b_dim[0], l1b_dim[1]], dtype=Config.datatype)
                #    glint_coeff[mask_valid_geo_water] = get_glint_coeff(l1b.solz[mask_valid_geo_water], \
                #                                                        l1b.senz[mask_valid_geo_water], \
                #                                                        l1b.relaz[mask_valid_geo_water], \
                #                                                        l1b_ws[mask_valid_geo_water])
                #    mask_glint = glint_coeff > glint_max
                #    mask_nglint = glint_coeff <= glint_max   
                #    mask_valid_geo_water_lrcposi_nglint = mask_valid_geo_water_lrcposi & mask_nglint
                #    l2_mask[mask_valid_geo_water_lrcposi & mask_glint] = 128
                    
                    # cloud mask, use bands near 412,555,670,865 
                    cmask = np.zeros([l1b_dim[0], l1b_dim[1]], dtype='bool')
                    if sinfo.sensor in ['OLI', 'OLI2']:
                        cmask = l1b.cloud
                    else:
                        cmask[mask_valid_geo_water_lrcposi]=cm.run_cloudmask(lrc[mask_valid_geo_water_lrcposi,:])
                    #expand the cloud mask by 1 pixel
    #                if sinfo.sensor in ['EPIC']:
    #                    cmask_idx = np.concatenate((np.where(cmask==1)[0]+1,np.where(cmask==1)[0]-1))
    #                    cmask_idy = np.concatenate((np.where(cmask==1)[1],np.where(cmask==1)[1]))
    #                    cmask[cmask_idx,cmask_idy] = True
    #                    cmask_idx = np.concatenate((np.where(cmask==1)[0],np.where(cmask==1)[0]))
    #                    cmask_idy = np.concatenate((np.where(cmask==1)[1]+1,np.where(cmask==1)[1]-1))
    #                    cmask[cmask_idx,cmask_idy] = True
                    mask_valid_geo_water_lrcposi_cloud = mask_valid_geo_water_lrcposi & cmask
                    mask_valid_geo_water_lrcposi_nocloud = mask_valid_geo_water_lrcposi & ~cmask
                    l2_mask[mask_valid_geo_water_lrcposi_cloud] = 64

                    training_band_ids = [i for i, band in enumerate(sinfo.band) if band in sinfo.training_bands]
                    wavelength_mask = np.full(lrc.shape[2], False, dtype=bool)
                    wavelength_mask[training_band_ids] = True

                    # The conditional statements for the following are so that we can save memory if the user doesn't request that product in OCSMART_Input.txt
                    aods = np.full([l1b_dim[0], l1b_dim[1], int(mlnn.aodnn_layers[-1])], np.nan, dtype=Config.datatype) if 'aod' in l2_prod else None
                    aph = np.full([l1b_dim[0], l1b_dim[1], int(mlnn.aphnn_layers[-1])], np.nan, dtype=Config.datatype) if 'aph' in l2_prod else None
                    adg = np.full([l1b_dim[0], l1b_dim[1], int(mlnn.adgnn_layers[-1])], np.nan, dtype=Config.datatype) if 'adg' in l2_prod else None
                    bbp = np.full([l1b_dim[0], l1b_dim[1], int(mlnn.bbpnn_layers[-1])], np.nan, dtype=Config.datatype) if 'bbp' in l2_prod else None
                    ap = np.full([l1b_dim[0], l1b_dim[1], int(mlnn.apnn_layers[-1])], np.nan, dtype=Config.datatype) if 'at' in l2_prod else None
                    bp = np.full([l1b_dim[0], l1b_dim[1], int(mlnn.bpnn_layers[-1])], np.nan, dtype=Config.datatype) if 'bt' in l2_prod else None

                    # run Multilayer Neural Network (MLNN) retrieval on Lrc data
                    #if image is too large, separate into blocks to process

                    oos_flag = np.zeros([l1b_dim[0], l1b_dim[1]], dtype='bool')
                    rrs = np.full([l1b_dim[0],l1b_dim[1], int(mlnn.rrsnn_layers[-1])], np.nan, dtype=Config.datatype)
                    # lrc_aann = np.full([l1b_dim[0], l1b_dim[1], int(mlnn.aann_layers[-1])], np.nan, dtype=Config.datatype)  # lrc_aann is an output of compute_aann, but is only actually used to compute the other output, oos_flag. No need to allocate lrc_aann and we can save memory

                    if block_size < 0:
                        rows, cols = np.where(mask_valid_geo_water_lrcposi_nocloud)
                        band_idx = np.where(wavelength_mask)[0]
                        lrc_masked_bands = lrc[rows[:, None], cols[:, None], band_idx[None, :]]

                        oos_flag[mask_valid_geo_water_lrcposi_nocloud], _ = \
                                                                                        mlnn.compute_aann(l1b.solz[mask_valid_geo_water_lrcposi_nocloud],\
                                                                                        l1b.senz[mask_valid_geo_water_lrcposi_nocloud],\
                                                                                        l1b.relaz[mask_valid_geo_water_lrcposi_nocloud],\
                                                                                        lrc_masked_bands,\
                                                                                        l1b_rh[mask_valid_geo_water_lrcposi_nocloud])
                        if 'aod' in l2_prod:
                            aods[mask_valid_geo_water_lrcposi_nocloud,:] = mlnn.compute_aodnn(l1b.solz[mask_valid_geo_water_lrcposi_nocloud],\
                                                                                            l1b.senz[mask_valid_geo_water_lrcposi_nocloud],\
                                                                                            l1b.relaz[mask_valid_geo_water_lrcposi_nocloud],\
                                                                                            lrc_masked_bands,\
                                                                                            l1b_rh[mask_valid_geo_water_lrcposi_nocloud])
                        # rrs must be retrieved        
                        # if sinfo.sensor == 'HYPSO_HSI':
                        #     lrc = lrc/4
                        rrs[mask_valid_geo_water_lrcposi_nocloud,:] = mlnn.compute_rrsnn(l1b.solz[mask_valid_geo_water_lrcposi_nocloud],\
                                                                                        l1b.senz[mask_valid_geo_water_lrcposi_nocloud],\
                                                                                        l1b.relaz[mask_valid_geo_water_lrcposi_nocloud],\
                                                                                        lrc_masked_bands)
                        del lrc_masked_bands, rows, cols, band_idx

                        if 'aph' in l2_prod:
                            aph[mask_valid_geo_water_lrcposi_nocloud,:] = mlnn.compute_aphnn(rrs[mask_valid_geo_water_lrcposi_nocloud,:])
                        if 'adg' in l2_prod:
                            adg[mask_valid_geo_water_lrcposi_nocloud,:] = mlnn.compute_adgnn(rrs[mask_valid_geo_water_lrcposi_nocloud,:])
                        if 'bbp' in l2_prod:
                            bbp[mask_valid_geo_water_lrcposi_nocloud,:] = mlnn.compute_bbpnn(rrs[mask_valid_geo_water_lrcposi_nocloud,:])
                        if 'bt' in l2_prod:
                            bp[mask_valid_geo_water_lrcposi_nocloud,:] = mlnn.compute_bpnn(rrs[mask_valid_geo_water_lrcposi_nocloud,:])
                        if 'at' in l2_prod:
                            ap[mask_valid_geo_water_lrcposi_nocloud,:] = mlnn.compute_apnn(rrs[mask_valid_geo_water_lrcposi_nocloud,:])
                    else:
                        blockmask = np.zeros([l1b_dim[0], l1b_dim[1]], dtype='bool')
                        nblocks_x = int(np.ceil(l1b_dim[1] / block_size))
                        nblocks_y = int(np.ceil(l1b_dim[0] / block_size))
                        print('Processing image in {} blocks ... '.format(nblocks_x*nblocks_y))
                        block_boundary_x = np.zeros(nblocks_x + 1, dtype=int)
                        block_boundary_y = np.zeros(nblocks_y + 1, dtype=int)
                        if nblocks_x == 1:
                            block_boundary_x[0] = 0
                            block_boundary_x[1] = l1b_dim[1]
                        else:
                            block_boundary_x[0:nblocks_x] = np.arange(0,l1b_dim[1], block_size)
                            block_boundary_x[nblocks_x] = l1b_dim[1]
                        if nblocks_y == 1:
                            block_boundary_y[0] = 0
                            block_boundary_y[1] = l1b_dim[0]
                        else:
                            block_boundary_y[0:nblocks_y] = np.arange(0,l1b_dim[0], block_size)
                            block_boundary_y[nblocks_y] = l1b_dim[0]
                        for iby in np.arange(nblocks_y):
                            for ibx in np.arange(nblocks_x):                            
                                print('Porcessing block ',iby*nblocks_y + ibx + 1)
                                blockmask[:, :] = False
                                blockmask[block_boundary_y[iby]:block_boundary_y[iby + 1], block_boundary_x[ibx]:block_boundary_x[ibx + 1]] = True                        
                                process_mask = mask_valid_geo_water_lrcposi_nocloud & blockmask
                                if np.sum(process_mask) > 0:

                                    oos_flag[process_mask], _ = mlnn.compute_aann(l1b.solz[process_mask],\
                                                                            l1b.senz[process_mask],\
                                                                            l1b.relaz[process_mask],\
                                                                            lrc[process_mask,:],\
                                                                            l1b_rh[process_mask])
                                    if 'aod' in l2_prod: 
                                        aods[process_mask,:] = mlnn.compute_aodnn(l1b.solz[process_mask],\
                                                                                l1b.senz[process_mask],\
                                                                                l1b.relaz[process_mask],\
                                                                                lrc[process_mask,:],\
                                                                                l1b_rh[process_mask])
                                    # rrs must be retrieved
                                    rrs[process_mask,:] = mlnn.compute_rrsnn(l1b.solz[process_mask],\
                                                                            l1b.senz[process_mask],\
                                                                            l1b.relaz[process_mask],\
                                                                            lrc[process_mask,:])
                                    if 'aph' in l2_prod:
                                        aph[process_mask,:] = mlnn.compute_aphnn(rrs[process_mask,:])
                                    if 'adg' in l2_prod:
                                        adg[process_mask,:] = mlnn.compute_adgnn(rrs[process_mask,:])
                                    if 'bbp' in l2_prod:
                                        bbp[process_mask,:] = mlnn.compute_bbpnn(rrs[process_mask,:])
                                    if 'at' in l2_prod:
                                        ap[process_mask,:] = mlnn.compute_apnn(rrs[process_mask,:])
                                    if 'bt' in l2_prod:
                                        bp[process_mask,:] = mlnn.compute_bpnn(rrs[process_mask,:])
                
                    mask_valid_geo_water_lrcposi_nocloud_oos = mask_valid_geo_water_lrcposi_nocloud & oos_flag
                    l2_mask[mask_valid_geo_water_lrcposi_nocloud_oos] = 256

                    # Variables are explicitly deleted after writing to the L2 file so that they don't persist while processing subsequent files and use up a bunch of memory

                    if 'Lrc' in l2_prod:
                        gw = hf.create_group('Lrc')
                        for i in np.arange(mlnn.aodnn_layers[-1]):
                            gw.create_dataset('Lrc_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = lrc[:,:,i], **compression_kwargs)
                    del lrc

                    if 'aod' in l2_prod:
                        gw = hf.create_group('AOD')
                        for i in np.arange(mlnn.aodnn_layers[-1]):
                            gw.create_dataset('AOD_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data=aods[:,:,i], **compression_kwargs)
                    del aods

                    if 'aph' in l2_prod:
                        gw = hf.create_group('aph')    
                        for i in np.arange(mlnn.aphnn_layers[-1]):
                            gw.create_dataset('aph_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data=aph[:,:,i], **compression_kwargs)
                    del aph

                    if 'adg' in l2_prod:
                        gw = hf.create_group('adg')    
                        for i in np.arange(mlnn.adgnn_layers[-1]):
                            gw.create_dataset('adg_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data=adg[:,:,i], **compression_kwargs)
                    del adg

                    if 'bbp' in l2_prod:
                        gw = hf.create_group('bbp')    
                        for i in np.arange(mlnn.bbpnn_layers[-1]):
                            gw.create_dataset('bbp_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data=bbp[:,:,i], **compression_kwargs)
                    del bbp

                    if 'at' in l2_prod:
                        gw = hf.create_group('at')    
                        for i in np.arange(mlnn.apnn_layers[-1]):
                            gw.create_dataset('at_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data=ap[:,:,i], **compression_kwargs)
                    del ap

                    if 'bt' in l2_prod:
                        gw = hf.create_group('bt')    
                        for i in np.arange(mlnn.bpnn_layers[-1]):
                            gw.create_dataset('bt_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data = bp[:,:,i], **compression_kwargs)
                    del bp
                    
                    # retrieve CHL, CDOM and TSM
                    
                    # chl_ocx = np.full([l1b_dim[0], l1b_dim[1]], np.nan, dtype=Config.datatype)
                    

                    # chl_ocx[mask_valid_geo_water_lrcposi_nocloud]=chl.get_chl_ocx(rrs[mask_valid_geo_water_lrcposi_nocloud,:])
                    if 'chl' in l2_prod:
                        chl_oci = np.full([l1b_dim[0], l1b_dim[1]], np.nan, dtype=Config.datatype)
                        chl_oci[mask_valid_geo_water_lrcposi_nocloud] = chl.get_chl_oci(rrs[mask_valid_geo_water_lrcposi_nocloud, :])
                        hf.create_dataset('chlor_a(oci)', dtype=Config.datatype,data = chl_oci, **compression_kwargs)
                        del chl_oci

                        chl_yoc = np.full([l1b_dim[0], l1b_dim[1]], np.nan, dtype=Config.datatype)
                        chl_yoc[mask_valid_geo_water_lrcposi_nocloud] = chl.get_chl_yoc(rrs[mask_valid_geo_water_lrcposi_nocloud, :])
                        hf.create_dataset('chlor_a(yoc)', dtype=Config.datatype,data = chl_yoc, **compression_kwargs)
                        del chl_yoc

                    if 'tsm' in l2_prod:
                        tsm_yoc = np.full([l1b_dim[0],l1b_dim[1]], np.nan, dtype=Config.datatype)
                        tsm_yoc[mask_valid_geo_water_lrcposi_nocloud] = tsm.get_tsm_yoc(rrs[mask_valid_geo_water_lrcposi_nocloud, :])
                        hf.create_dataset('tsm(yoc)', dtype=Config.datatype, data = tsm_yoc, **compression_kwargs)
                        del tsm_yoc
                    
                    hf.create_dataset('Latitude', dtype='float32', data = l1b.latitude, **compression_kwargs)
                    del l1b.latitude
                    hf.create_dataset('Longitude', dtype='float32', data = l1b.longitude, **compression_kwargs)
                    del l1b.longitude
                    hf.create_dataset('Solar_zenith', dtype='float32', data = l1b.solz, **compression_kwargs)
                    del l1b.solz
                    hf.create_dataset('Sensor_zenith', dtype='float32', data = l1b.senz, **compression_kwargs)
                    del l1b.senz
                    hf.create_dataset('Relative_azimuth', dtype='float32', data = l1b.relaz, **compression_kwargs)
                    del l1b.relaz
                    hf.create_dataset('L2_flags', dtype='int16', data = l2_mask,**compression_kwargs)
                    del l2_mask

                    if 'pressure' in l2_prod:
                        hf.create_dataset('pressure', dtype=Config.datatype, data = l1b_pressure[:, :], **compression_kwargs)
                    del l1b_pressure

                    if 'relative_humidity' in l2_prod:
                        hf.create_dataset('relative_humidity', dtype=Config.datatype, data = l1b_rh[:, :], **compression_kwargs)
                    del l1b_rh

                    if 'wind_speed' in l2_prod:
                        hf.create_dataset('wind_speed', dtype=Config.datatype, data = l1b_ws[:, :], **compression_kwargs)
                    del l1b_ws 

                    if 'rrs' in l2_prod:
                        gw = hf.create_group('Rrs')    
                        for i in np.arange(mlnn.rrsnn_layers[-1]):
                            gw.create_dataset('Rrs_' + str(training_bands[i]) + 'nm', dtype=Config.datatype, data=rrs[:,:,i], **compression_kwargs)
                    del rrs

                    print('Finished writing level-2 file {} ... '.format(os.path.splitext(fname)[0] + '_L2_OCSMART.h5'))
                
                print('Processing finished in %.2f second.\n'%(time.time()-t_start))
            else:
                print('\033[1;31;47mWARNING:Unable to extract subimage, processing terminated... ', '\033[m')
                continue
        else:
            continue
    else:
        continue
