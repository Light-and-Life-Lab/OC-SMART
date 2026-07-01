#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Wed Jan 16 11:09:04 2019

@author: Yongzhen Fan
"""

import numpy as np
from os.path import basename
import sys

class sensorinfo(object):
    
    def __init__(self, L1Bname, sensor=None):
        self.sensor=sensor
        self.sat=None
        self.l1bname=L1Bname
        self.l1bbasename = basename(L1Bname)
        self.sensor_status = 0
        if sensor is None:
            self.autodetect()
        if self.sensor_status == 0:          
            info_fname='./auxdata/sensorinfo/'+self.sensor+'.txt'
            info=np.loadtxt(info_fname,dtype=np.float64)
            self.band=info[:,0].astype(int)
            if self.sensor == 'OCI':
                NN_wavelengths_fname = './auxdata/sensorinfo/' + self.sensor + '_NN_wavelengths.txt'
                NN_wavelengths_file = np.loadtxt(NN_wavelengths_fname, dtype=np.float64)
                self.training_bands = NN_wavelengths_file[:].astype(int)
            else:
                self.training_bands = self.band

            # self.training_bands = self.band

            self.koz=info[:,1]
            self.tauray=info[:,2]
            self.kno2=info[:,3]
            self.vgain=info[:,4]
            self.vgaino=info[:,5] 
            self.vgainc=info[:,6] 

    def autodetect(self):
        
        b = self.l1bbasename
        self.datalevel = None
        self.datasource = None

        if (b.startswith('MER_RR') or b.startswith('MER_FR')) and b.endswith('.N1'):
            self.sensor = 'MERIS'
            self.sat = 'MERIS'
            self.gasid = 111

        elif b.startswith('S3A_OL_1') and b.endswith('.SEN3'):
            self.sensor = 'OLCI'
            self.sat = 'S3A'
            self.gasid = 111
        
        elif b.startswith('S3B_OL_1') and b.endswith('.SEN3'):
            self.sensor = 'OLCI'
            self.sat = 'S3B'
            self.gasid = 111

        elif b.startswith('SNPP_VIIRS') and b.endswith('.nc'):
            self.sensor = 'VIIRS'
            self.datasource = 'OBPG'
            self.sat = 'SNPP'
            self.gasid = 127
        
        elif b.startswith('JPSS1_VIIRS') and b.endswith('.nc'):
            self.sensor = 'VIIRS'
            self.datasource = 'OBPG'
            self.sat = 'JPSS1'
            self.gasid = 127
        
        elif b.startswith('JPSS2_VIIRS') and b.endswith('.nc'):
            self.sensor = 'VIIRS'
            self.datasource = 'OBPG'
            self.sat = 'JPSS2'
            self.gasid = 127
            
        # elif b.startswith('NPP') and b.endswith('.hdf'):
            # self.sensor = 'VIIRS'
            # self.datasource = 'LAADS DAAC'
            # self.sat = 'SNPP'
        
        elif b.startswith('VNP02MOD') and b.endswith('.nc'):
            self.sensor = 'VIIRS'
            self.datasource = 'LAADS DAAC'
            self.sat = 'SNPP'
            self.gasid = 127
            
        elif b.endswith('noaa_ops.h5'):
            self.sensor = 'VIIRS'
            self.datasource = 'NOAA'
            self.sat = 'VIIRS'
            self.gasid = 127

        elif b.startswith('AQUA_MODIS') and b.endswith('.L1B.hdf'):
            self.sensor = 'MODIS-Aqua'
            self.datasource = 'OBPG'
            self.sat = 'MODISA'
        
        elif b.startswith('MYD021KM') and b.endswith('.hdf'):
            self.sensor = 'MODIS-Aqua'
            self.datasource = 'LAADS DAAC'
            self.sat = 'MODISA'
            self.gasid = 255
        
        elif b.startswith('TERRA_MODIS') and b.endswith('.L1B.hdf'):
            self.sensor = 'MODIS-Terra'
            self.datasource = 'OBPG'
            self.sat = 'MODIST'
            self.gasid = 127
        
        elif b.startswith('MOD021KM') and b.endswith('.hdf'):
            self.sensor = 'MODIS-Terra'
            self.datasource = 'LAADS DAAC'
            self.sat = 'MODIST'
            self.gasid = 127

        elif b.startswith('S') and '.L1B' in b:
            self.sensor = 'SeaWiFS'
            self.sat = 'SeaWiFS'
            self.gasid = 13

        elif b.startswith('COMS_GOCI_L1B') and b.endswith('.he5'):
            self.sensor = 'GOCI'
            self.sat = 'GOCI'
            self.gasid = 127
            
        elif b.startswith('GC1SG1') and b.endswith('.h5'):
            self.sensor = 'SGLI'
            self.sat = 'SGLI'
            self.gasid = 5
            
        elif b.startswith('epic_1b') and b.endswith('.h5'):
            self.sensor = 'EPIC'
            self.sat = 'EPIC'
            self.gasid = 5
            
        elif b.startswith('LC8') or b.startswith('LC08'):
            self.sensor = 'OLI'
            self.sat = 'L08'
            self.gasid = 5
            
        elif b.startswith('LC9') or b.startswith('LC09'):
            self.sensor = 'OLI2'
            self.sat = 'L09'
            self.gasid = 5
            
        elif b.startswith('S2A_MSIL1C'):
            self.sensor = 'S2A'
            self.datalevel = 'L1C'
            self.sat = 'S2A'
            self.gasid = 5
            
        elif b.startswith('S2A_MSIL2A'):
            self.sensor = 'S2A'
            self.datalevel = 'L2A'
            self.sat = 'S2A'
            self.gasid = 5
            
        elif b.startswith('S2B_MSIL1C'):
            self.sensor = 'S2B'
            self.datalevel = 'L1C'
            self.sat = 'S2B'
            self.gasid = 5
            
        elif b.startswith('S2B_MSIL2A'):
            self.sensor = 'S2B'
            self.datalevel = 'L2A'
            self.sat = 'S2B'
            self.gasid = 5
            
        elif b.startswith('FY3D') and b.endswith('1000M_MS.HDF'):
            self.sensor = 'MERSI2'
            self.sat = 'MERSI2'
            self.gasid = 5
        
        elif b.startswith('H') and b.endswith('ISS.nc'):
            self.sensor = 'HICO'
            self.datasource = 'OBPG'
            self.sat = 'HICO'
            self.gasid = 5
        
        elif b.startswith('PACE_OCI') and b.endswith('.nc'):
            self.sensor = 'OCI'
            self.datasource = 'NASA_DAAC'
            self.sat = 'PACE'
            self.gasid = 127

        elif b.startswith('HYPSO_HSI') and b.endswith('.nc'):
            self.sensor = 'HYPSO_HSI'
            self.datasource = 'NASA_DAAC' # TODO: Update this if necessary
            self.sat = 'HYPSO'
            self.gasid = 127
            
        else:
            print('\033[1;31;47mWARNING: Unable to detect sensor from file "{}", processing terminated ... \n'.format(b),'\033[m')
            self.sensor_status=1
        