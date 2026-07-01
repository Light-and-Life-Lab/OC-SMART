#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 11:05:04 2019

@author: Yongzhen Fan
"""
import os
from os.path import exists, dirname, basename, isdir
from os import makedirs
import numpy as np
from netCDF4 import Dataset
import h5py
from pyhdf.SD import SD, SDC
from datetime import datetime, timedelta
import bz2
import urllib.request
import urllib.error
from scipy import interpolate
import time
import calendar
from src.l8_utils import load_mtl
from glob import glob
from lxml import objectify
from src.obdaac_download import httpdl

from lib.gas_corrections_lib.src import gas_corrections

class ANCILLARY(object):

    def __init__(self, l1b_data):
        self.datestr = l1b_data.datestr
        self.path='./anc/'
        self.metsource='MERRA2' # NASA OBPG ancillary data source, options: NCEP, MERRA2
        self.gas_transmittance_manager = gas_corrections.Gas_Correction_Manager()

        print('Locating ancillary files ...')
        if not isdir(self.path):
            print('Ancillary directory does not exist, creating ancillary directory {} ... '.format(self.path))
            makedirs(self.path)
            

    def read_gas_transmittance_auxdata(self, sensor_info):
        auxpath = './auxdata/common'
        if sensor_info.sensor == 'GOCI':
            self.gas_transmittance_filepath = auxpath + '/' + "goci_gas_transmittance.nc"
        elif sensor_info.sensor == 'HICO':
            self.gas_transmittance_filepath = auxpath + '/' + "hico_gas_transmittance.nc"
        elif sensor_info.sensor == 'HYPSO_HSI':
            self.gas_transmittance_filepath = auxpath + '/' + "oci_gas_transmittance_cia_amf_v3.3.nc"  # Using OCI transmittance data for HYPSO
        elif sensor_info.sensor == 'MERIS':
            self.gas_transmittance_filepath = auxpath + '/' + "meris_gas_transmittance.nc"
        elif sensor_info.sensor == 'MODIS-Aqua':
            self.gas_transmittance_filepath = auxpath + '/' + "modisa_gas_transmittance.nc"
        elif sensor_info.sensor == 'MODIS-Terra':
            self.gas_transmittance_filepath = auxpath + '/' + "modist_gas_transmittance.nc"
        elif sensor_info.sensor == 'OCI':
            self.gas_transmittance_filepath = auxpath + '/' + "oci_gas_transmittance_cia_amf_v3.3.nc"
        elif sensor_info.sensor == 'OLCI' and sensor_info.sat == 'S3A':
            self.gas_transmittance_filepath = auxpath + '/' + "olcis3a_gas_transmittance.nc"
        elif sensor_info.sensor == 'OLCI' and sensor_info.sat == 'S3B':
            self.gas_transmittance_filepath = auxpath + '/' + "olcis3b_gas_transmittance.nc"
        elif sensor_info.sensor == 'OLI':
            self.gas_transmittance_filepath = auxpath + '/' + "olil8_gas_transmittance.nc"
        elif sensor_info.sensor == 'SeaWiFS':
            self.gas_transmittance_filepath = auxpath + '/' + "seawifs_gas_transmittance.nc"
        elif sensor_info.sensor == 'VIIRS' and sensor_info.sat == 'JPSS1':
            self.gas_transmittance_filepath = auxpath + '/' + "viirsj1_gas_transmittance.nc"
        elif sensor_info.sensor == 'VIIRS' and sensor_info.sat == 'JPSS2':
            self.gas_transmittance_filepath = auxpath + '/' + "viirsj2_gas_transmittance.nc"
        elif sensor_info.sensor == 'VIIRS' and sensor_info.sat == 'SNPP':
            self.gas_transmittance_filepath = auxpath + '/' + "viirsn_gas_transmittance.nc"
        # else:
        #     raise FileExistsError(f"Transmittance data is not available for {sensor_info.sensor} on satellite {sensor_info.sat}")
        
        # self.gas_transmittance_manager.read_gas_transmittance_table(self.gas_transmittance_filepath)

            
    def read_no2(self): 
        auxpath='./auxdata/common/'
        print('Reading NO2 data ...')        
        #read NO2 data       
        months=range(1,13)
        nmonths=12
        #set latitude and longitude grid
        self.no2_frac_lat=np.arange(91,-93,-2)
        self.no2_frac_lon=np.arange(-181,183,2)
        self.no2_lat=np.arange(90.125,-90.375,-0.25)
        self.no2_lon=np.arange(-180.125,180.375,0.25)
        no2_frac_nline=len(self.no2_frac_lat)
        no2_frac_npixl=len(self.no2_frac_lon)
        no2_nline=len(self.no2_lat)
        no2_npixl=len(self.no2_lon)
        
        self.no2_total = np.zeros((nmonths,no2_nline,no2_npixl), dtype='float64')
        self.no2_tropo = np.zeros((nmonths,no2_nline,no2_npixl), dtype='float64')
        self.no2_strat = np.zeros((nmonths,no2_nline,no2_npixl), dtype='float64')
        self.no2_frac = np.zeros((no2_frac_nline,no2_frac_npixl), dtype='float64')
        
        no2_fname=auxpath+'no2_climatology_v2013.hdf'
        no2_frac_fname=auxpath+'trop_f_no2_200m.hdf'
        
        #read no2 fraction data
        f=SD(no2_frac_fname, SDC.READ)
        self.no2_frac[1:no2_frac_nline-1,1:no2_frac_npixl-1] = f.select('f_no2_200m')[:,:]
        self.no2_frac[:,0]=self.no2_frac[:,no2_frac_npixl-2]
        self.no2_frac[:,no2_frac_npixl-1]=self.no2_frac[:,1]
        self.no2_frac[0,:]=self.no2_frac[1,:]
        self.no2_frac[no2_frac_nline-1,:]=self.no2_frac[no2_frac_nline-2,:]
        
        # read total and tropospheric no2 data
        f=SD(no2_fname, SDC.READ)
        for i, m in enumerate(months):               
            self.no2_tropo[i,1:no2_nline-1,1:no2_npixl-1] = f.select('trop_no2_{:02d}'.format(m))[:,:]
            self.no2_tropo[i,:,0]=self.no2_tropo[i,:,no2_npixl-2]
            self.no2_tropo[i,:,no2_npixl-1]=self.no2_tropo[i,:,1]
            self.no2_tropo[i,0,:]=self.no2_tropo[i,1,:]
            self.no2_tropo[i,no2_nline-1,:]=self.no2_tropo[i,no2_nline-2,:]
            self.no2_total[i,1:no2_nline-1,1:no2_npixl-1] = f.select('tot_no2_{:02d}'.format(m))[:,:]
            self.no2_total[i,:,0]=self.no2_total[i,:,no2_npixl-2]
            self.no2_total[i,:,no2_npixl-1]=self.no2_total[i,:,1]
            self.no2_total[i,0,:]=self.no2_total[i,1,:]
            self.no2_total[i,no2_nline-1,:]=self.no2_total[i,no2_nline-2,:]
            self.no2_strat[i,:,:]=self.no2_total[i,:,:]-self.no2_tropo[i,:,:]
        self.no2_strat[self.no2_strat<0.0]=0.0
        self.no2_total=self.no2_total * 1.0e15 
        self.no2_tropo=self.no2_tropo * 1.0e15
        self.no2_strat=self.no2_strat * 1.0e15        

    def download(self):
#        t=time.time()
#        ancurl_prefix='https://oceandata.sci.gsfc.nasa.gov/ob/getfile/' #NASA ancillary data
        server = 'oceandata.sci.gsfc.nasa.gov'
        self.l1btime=datetime.strptime(self.datestr,'%Y%m%d%H%M%S').timestamp()
        self.doy=datetime.strptime(self.datestr,'%Y%m%d%H%M%S').timetuple().tm_yday
        
        if self.metsource == 'NCEP':
            oz_postfix='_O3_AURAOMI_24h.hdf'
            oz_postfix2='_O3_EPTOMS_24h.hdf'
            met_postfix='_MET_NCEP_6h.hdf'            
            
            if(int(self.datestr[8:14])<120000): #before mid day
                dt=(datetime.strptime(self.datestr[0:8]+'000000','%Y%m%d%H%M%S')-timedelta(days=1)).timetuple()
                self.oz1_name='N'+str(dt.tm_year)+'{:03d}'.format(dt.tm_yday)+'00'+oz_postfix
                self.oz2_name='N'+self.datestr[0:4]+'{:03d}'.format(self.doy)+'00'+oz_postfix
                self.ozdt=(self.l1btime-(datetime.strptime(self.datestr[0:8]+'000000','%Y%m%d%H%M%S').timestamp()-12*3600))/3600/24
            else:
                dt=(datetime.strptime(self.datestr[0:8]+'000000','%Y%m%d%H%M%S')+timedelta(days=1)).timetuple()
                self.oz1_name='N'+self.datestr[0:4]+'{:03d}'.format(self.doy)+'00'+oz_postfix
                self.oz2_name='N'+str(dt.tm_year)+'{:03d}'.format(dt.tm_yday)+'00'+oz_postfix
                self.ozdt=(self.l1btime-datetime.strptime(self.datestr[0:8]+'120000','%Y%m%d%H%M%S').timestamp())/3600/24
          
            if(int(self.datestr[8:14])<=60000): #before 06:00
                self.met1_name='N'+self.datestr[0:4]+'{:03d}'.format(self.doy)+'00'+met_postfix
                self.met2_name='N'+self.datestr[0:4]+'{:03d}'.format(self.doy)+'06'+met_postfix
                self.metdt=(self.l1btime-datetime.strptime(self.datestr[0:8]+'000000','%Y%m%d%H%M%S').timestamp())/3600/24/0.25
            elif(int(self.datestr[8:14])<=120000): #before 12:00
                self.met1_name='N'+self.datestr[0:4]+'{:03d}'.format(self.doy)+'06'+met_postfix
                self.met2_name='N'+self.datestr[0:4]+'{:03d}'.format(self.doy)+'12'+met_postfix
                self.metdt=(self.l1btime-datetime.strptime(self.datestr[0:8]+'060000','%Y%m%d%H%M%S').timestamp())/3600/24/0.25
            elif(int(self.datestr[8:14])<=180000): #before 18:00
                self.met1_name='N'+self.datestr[0:4]+'{:03d}'.format(self.doy)+'12'+met_postfix
                self.met2_name='N'+self.datestr[0:4]+'{:03d}'.format(self.doy)+'18'+met_postfix
                self.metdt=(self.l1btime-datetime.strptime(self.datestr[0:8]+'120000','%Y%m%d%H%M%S').timestamp())/3600/24/0.25
            else: #after 18:00
                dt=(datetime.strptime(self.datestr[0:8]+'000000','%Y%m%d%H%M%S')+timedelta(days=1)).timetuple()
                self.met1_name='N'+self.datestr[0:4]+'{:03d}'.format(self.doy)+'18'+met_postfix
                self.met2_name='N'+str(dt.tm_year)+'{:03d}'.format(dt.tm_yday)+'00'+met_postfix
                self.metdt=(self.l1btime-datetime.strptime(self.datestr[0:8]+'180000','%Y%m%d%H%M%S').timestamp())/3600/24/0.25
            
            if exists(self.path+self.oz1_name):
                print('Ozone file {} located on local drive.'.format(self.oz1_name))
            elif exists(self.path+self.oz1_name[0:10]+oz_postfix2):
                self.oz1_name = self.oz1_name[0:10]+oz_postfix2
                print('Ozone file {} located on local drive.'.format(self.oz1_name))
            else:
                print('Ozone file {} not found on local drive, downloading from NASA OBPG ...'.format(self.oz1_name))        
                request='/getfile/'+self.oz1_name
                status = httpdl(server, request, localpath=self.path, uncompress=False)
                if status:
                    print('OMI ozone data unavailable, downloading TOMS ozone data ...')
                    self.oz1_name = self.oz1_name[0:10]+oz_postfix2 #rename oz file
                    request='/getfile/'+self.oz1_name
                    status = httpdl(server, request, localpath=self.path, uncompress=False)
                    if status:
                        print('Real time ozone data unavailable, using climatology data ...') 
                        self.oz1_name = 'ozone_climatology_v2014.hdf'#rename oz file to the climatology data                    
                   
            if exists(self.path+self.oz2_name):
                print('Ozone file {} located on local drive.'.format(self.oz2_name))
            elif exists(self.path+self.oz2_name[0:10]+oz_postfix2):
                self.oz2_name = self.oz2_name[0:10]+oz_postfix2
                print('Ozone file {} located on local drive.'.format(self.oz2_name))
            else:
                print('Ozone file {} not found on local drive, downloading from NASA OBPG ...'.format(self.oz2_name))
                request='/getfile/'+self.oz2_name
                status = httpdl(server, request, localpath=self.path, uncompress=False)
                if status:
                    print('OMI ozone data unavailable, downloading TOMS ozone data ...')
                    self.oz2_name = self.oz2_name[0:10]+oz_postfix2 #rename oz file
                    request='/getfile/'+self.oz2_name
                    status = httpdl(server, request, localpath=self.path, uncompress=False)
                    if status:
                        print('Real time ozone data unavailable, using climatology data ...') 
                        self.oz2_name = 'ozone_climatology_v2014.hdf'#rename oz file to the climatology data       
                    
            if exists(self.path+self.met1_name):
                print('MET file {} located on local drive.'.format(self.met1_name))
            else:
                print('MET file {} not found on local drive, downloading from NASA OBPG ...'.format(self.met1_name)) 
                request='/getfile/'+self.met1_name
                status = httpdl(server, request, localpath=self.path, uncompress=False)
                if status:
                    print('Real time MET data unavailable ...')
    #            else:
    #                f=open(self.path+self.met1_name,'rb')
    #                compdata=f.read()
    #                f.close()
    #                decompdata=bz2.decompress(compdata)
    #                f=open(self.path+self.met1_name,'wb')
    #                f.write(decompdata)
    #                f.flush()
    #                os.remove(self.path+self.met1_name+'.bz2')

                
            if exists(self.path+self.met2_name):
                print('MET file {} located on local drive.'.format(self.met2_name))
            else:
                print('MET file {} not found on local drive, downloading from NASA OBPG ...'.format(self.met2_name))  
                request='/getfile/'+self.met2_name
                status = httpdl(server, request, localpath=self.path, uncompress=False)
                if status:
                    print('Real time MET data unavailable ...')
    #            else:
    #                f=open(self.path+self.met2_name,'rb')
    #                compdata=f.read()
    #                f.close()
    #                decompdata=bz2.decompress(compdata)
    #                f=open(self.path+self.met2_name,'wb')
    #                f.write(decompdata)
    #                f.flush()
    #                os.remove(self.path+self.met2_name+'.bz2')
        elif self.metsource == 'MERRA2':            
            time1_str = self.datestr[0:8] + 'T' + self.datestr[8:10] + '0000'
            time1 = datetime.strptime(time1_str[0:8]+time1_str[9:15],'%Y%m%d%H%M%S').timestamp()
            time2 = datetime.fromtimestamp(self.l1btime+3600)
            time2_str = str(time2.year) + str(time2.month).zfill(2) + str(time2.day).zfill(2) + 'T' + str(time2.hour).zfill(2) + '0000'
            self.ozdt = (self.l1btime - time1)/3600
            self.metdt = self.ozdt
            self.oz1_name='GMAO_MERRA2.'+ time1_str + '.MET.nc'
            self.oz2_name='GMAO_MERRA2.'+ time2_str + '.MET.nc'
            self.met1_name='GMAO_MERRA2.'+ time1_str + '.MET.nc'
            self.met2_name='GMAO_MERRA2.'+ time2_str + '.MET.nc'
            
            #download
            if exists(self.path+self.oz1_name):
                print('Ozone & Met file {} located on local drive.'.format(self.oz1_name))
            else:
                print('Ozone & MET file {} not found on local drive, downloading from NASA OBPG ...'.format(self.met1_name)) 
                request='/getfile/'+self.met1_name
                status = httpdl(server, request, localpath=self.path, uncompress=False)
                if status:
                    print('MERRA2 refined ancillary data unavailable, downloading MERRA2 IT NRT ancillary data ... ') 
                    self.oz1_name='GMAO_IT.'+ time1_str + '.MET.NRT.nc'
                    self.met1_name='GMAO_IT.'+ time1_str + '.MET.NRT.nc'
                    if not exists(self.path+self.oz1_name):                        
                        request='/getfile/'+self.met1_name
                        status = httpdl(server, request, localpath=self.path, uncompress=False)
                        if status:
                            print('Real time MET data unavailable ...')
            
            if exists(self.path+self.oz2_name):
                print('Ozone & Met file {} located on local drive.'.format(self.oz2_name))
            else:
                print('Ozone & MET file {} not found on local drive, downloading from NASA OBPG ...'.format(self.met2_name)) 
                request='/getfile/'+self.met2_name
                status = httpdl(server, request, localpath=self.path, uncompress=False)
                if status:
                    print('MERRA2 refined ancillary data unavailable, downloading MERRA2 IT NRT ancillary data ... ')
                    self.oz2_name='GMAO_IT.'+ time2_str + '.MET.NRT.nc'
                    self.met2_name='GMAO_IT.'+ time2_str + '.MET.NRT.nc'
                    if not exists(self.path+self.oz2_name):
                        request='/getfile/'+self.met2_name
                        status = httpdl(server, request, localpath=self.path, uncompress=False)
                        if status:
                            print('Real time MET data unavailable ...')
            
    def read_ozone(self):
        print('Reading Ozone data ...')
        if self.metsource == 'NCEP':
            self.oz_lat=np.arange(90.5,-91.5,-1)
            self.oz_lon=np.arange(-180.5,181.5,1)
            
            if 'climatology' in self.oz1_name or 'climatology' in self.oz2_name:
                print('Warning: real time ozone data unavailable, using climatology ozone data ...')
            if 'climatology' in self.oz1_name:
                f1 = SD('./auxdata/common/'+self.oz1_name,SDC.READ)
                ozone1 = f1.select('ozone_mean_'+'{:03d}'.format(self.doy))[:,:]*0.001
            elif 'OMI' in self.oz1_name:
                f1 = SD(self.path+self.oz1_name,SDC.READ)
                ozone1 = f1.select('ozone')[:,:]*0.001
            # elif 'TOMS' in self.oz1_name:
            #     f1 = SD(self.path+self.oz1_name,SDC.READ)
            #     data = f1.select('ozone')[:,:]*0.001
            #     #TOMS ozone data use a different grid, interpolate to 1 degree grid point
            #     lat = np.arange(89.5,-90.5,-1)
            #     lon = np.arange(-180.625,181.875,1.25)
                
            #     nline = len(lat)
            #     npix = len(lon)
            #     ozone = np.zeros((nline,npix), dtype='float64')
            #     ozone[:,1:npix-1] = data
            #     ozone[:,0] = data[:,-1]
            #     ozone[:,-1] = data[:,0]
            #     grid_yt=np.arange(-179.5,180.5,1)
            #     func = interpolate.interp2d(lon,lat,ozone,kind='linear')
            #     ozone1 = np.flip(np.flip(func(grid_yt,lat),1))
            
            if 'climatology' in self.oz2_name:
                f2 = SD('./auxdata/common/'+self.oz2_name,SDC.READ)
                ozone2 = f2.select('ozone_mean_'+'{:03d}'.format(self.doy))[:,:]*0.001
            elif 'OMI' in self.oz2_name:
                f2 = SD(self.path+self.oz2_name,SDC.READ)
                ozone2 = f2.select('ozone')[:,:]*0.001
            # elif 'TOMS' in self.oz2_name:
            #     f2 = SD(self.path+self.oz2_name,SDC.READ)
            #     data = f2.select('ozone')[:,:]*0.001
            #     #TOMS ozone data use a different grid, interpolate to 1 degree grid point
            #     lat = np.arange(89.5,-90.5,-1)
            #     lon = np.arange(-180.625,181.875,1.25)
            #     nline = len(lat)
            #     npix = len(lon)
            #     ozone = np.zeros((nline,npix), dtype='float64')
            #     ozone[:,1:npix-1] = data
            #     ozone[:,0] = data[:,-1]
            #     ozone[:,-1] = data[:,0]
            #     grid_yt=np.arange(-179.5,180.5,1)
            #     func = interpolate.interp2d(lon,lat,ozone,kind='linear')
            #     ozone2 = np.flip(np.flip(func(grid_yt,lat),1))
                
        elif self.metsource == 'MERRA2':
            self.oz_lat = np.arange(-90.5,91.0,0.5)
            self.oz_lon = np.arange(-180.3125,180.625,0.625)
            if 'climatology' in self.oz1_name or 'climatology' in self.oz2_name:
                print('Warning: real time ozone data unavailable, using climatology ozone data ...')
            if 'climatology' in self.oz1_name:
                f1 = SD('./auxdata/common/'+self.oz1_name,SDC.READ)
                ozone1 = f1.select('ozone_mean_'+'{:03d}'.format(self.doy))[:,:]*0.001
            else:
                fd=Dataset(self.path+self.oz1_name,'r')                
                ozone1 = fd.variables['TO3'][:]*0.001
                fd.close()
            
            if 'climatology' in self.oz2_name:
                f2 = SD('./auxdata/common/'+self.oz2_name,SDC.READ)
                ozone2 = f2.select('ozone_mean_'+'{:03d}'.format(self.doy))[:,:]*0.001
            else:
                fd=Dataset(self.path+self.oz2_name,'r')
                ozone2 = fd.variables['TO3'][:]*0.001
                fd.close()
                
        oz_nline=len(self.oz_lat)
        oz_npixl=len(self.oz_lon)
        ozone=(ozone1*(1-self.ozdt)+ozone2*self.ozdt) # interpolate in time and convert unit   
        self.ozmap=np.zeros((oz_nline,oz_npixl), dtype='float64')
        self.ozmap[1:oz_nline-1,1:oz_npixl-1]=ozone
        self.ozmap[:,0]=self.ozmap[:,oz_npixl-2]
        self.ozmap[:,oz_npixl-1]=self.ozmap[:,1]
        self.ozmap[0,:]=self.ozmap[1,:]
        self.ozmap[oz_nline-1,:]=self.ozmap[oz_nline-2,:]
        
    def read_met(self):
        print('Reading windspeed, pressure and RH data ...')
        if self.metsource =='NCEP':
            f1=SD(self.path+self.met1_name,SDC.READ)
            zwind=f1.select('z_wind')[:,:]
            mwind=f1.select('m_wind')[:,:]
            ws1=np.power(np.power(zwind,2)+np.power(mwind,2),0.5)
            press1=f1.select('press')[:,:]
            rh1=f1.select('rel_hum')[:,:]
            
            f2=SD(self.path+self.met2_name,SDC.READ)
            zwind=f2.select('z_wind')[:,:]
            mwind=f2.select('m_wind')[:,:]
            ws2=np.power(np.power(zwind,2)+np.power(mwind,2),0.5)
            press2=f2.select('press')[:,:]
            rh2=f2.select('rel_hum')[:,:]
            
            ws=ws1*(1-self.metdt)+ws2*self.metdt
            press=press1*(1-self.metdt)+press2*self.metdt
            rh=rh1*(1-self.metdt)+rh2*self.metdt        
            
            self.met_lat=np.arange(91,-92,-1)
            self.met_lon=np.arange(-180.5,181.5,1)
            met_nline=len(self.met_lat)
            met_npixl=len(self.met_lon)
            
        elif self.metsource == 'MERRA2':
            fd1 = Dataset(self.path+self.met1_name,'r')
            zwind = fd1.variables['U10M'][:]
            mwind = fd1.variables['V10M'][:]
            ws1 = np.power(np.power(zwind,2)+np.power(mwind,2),0.5)
            press1 = fd1.variables['PS'][:]*0.01 # convert to millibar
            t10 = fd1.variables['T10M'][:] # 10m temperature [K]
            sh = fd1.variables['QV10M'][:] # 10m specific humidity [kg/kg]
            rh1 = 0.263*sh*press1*100/(np.exp(17.67*(t10-273.15)/(t10-29.65))) #convert specific humidity to relative humidity
            rh1[rh1<0]=0.0
            rh1[rh1>100]=100.0
            total_precipitable_water1 = fd1.variables['TQV'][:]/0.1 # Convert from kg m^-2 to g cm^-2
            fd1.close()
            
            fd2 = Dataset(self.path+self.met2_name,'r')
            zwind = fd2.variables['U10M'][:]
            mwind = fd2.variables['V10M'][:]
            ws2 = np.power(np.power(zwind,2)+np.power(mwind,2),0.5)
            press2 = fd2.variables['PS'][:]*0.01 # convert to millibar
            t10 = fd2.variables['T10M'][:] # 10m temperature [K]
            sh = fd2.variables['QV10M'][:] # 10m specific humidity [kg/kg]
            rh2 = 0.263*sh*press2*100/(np.exp(17.67*(t10-273.15)/(t10-29.65))) #convert specific humidity to relative humidity
            rh2[rh2<0]=0.0
            rh2[rh2>100]=100.0
            total_precipitable_water2 = fd2.variables['TQV'][:]/0.1 # Convert from kg m^-2 to g cm^-2
            fd2.close()
            
            ws = ws1*(1 - self.metdt) + ws2*self.metdt
            press = press1*(1 - self.metdt) + press2*self.metdt
            rh = rh1*(1 - self.metdt) + rh2*self.metdt 
            total_precipitable_water = total_precipitable_water1*(1 - self.metdt) + total_precipitable_water2*self.metdt 
            
            self.met_lat = np.arange(-90.5,91.0,0.5)
            self.met_lon = np.arange(-180.3125,180.625,0.625)
            met_nline=len(self.met_lat)
            met_npixl=len(self.met_lon)
        
        self.wsmap = np.zeros((met_nline,met_npixl), dtype='float64')
        self.pressmap = np.zeros((met_nline,met_npixl), dtype='float64')
        self.rhmap = np.zeros((met_nline,met_npixl), dtype='float64')
        self.precipitable_water_map = np.zeros((met_nline,met_npixl), dtype='float64')
        
        self.wsmap[1:met_nline-1,1:met_npixl-1]=ws
        self.wsmap[:,0]=self.wsmap[:,met_npixl-2]
        self.wsmap[:,met_npixl-1]=self.wsmap[:,1]
        self.wsmap[0,:]=self.wsmap[1,:]
        self.wsmap[met_nline-1,:]=self.wsmap[met_nline-2,:]
        
        self.pressmap[1:met_nline-1,1:met_npixl-1]=press
        self.pressmap[:,0]=self.pressmap[:,met_npixl-2]
        self.pressmap[:,met_npixl-1]=self.pressmap[:,1]
        self.pressmap[0,:]=self.pressmap[1,:]
        self.pressmap[met_nline-1,:]=self.pressmap[met_nline-2,:]
        
        self.rhmap[1:met_nline-1,1:met_npixl-1]=rh
        self.rhmap[:,0]=self.rhmap[:,met_npixl-2]
        self.rhmap[:,met_npixl-1]=self.rhmap[:,1]
        self.rhmap[0,:]=self.rhmap[1,:]
        self.rhmap[met_nline-1,:]=self.rhmap[met_nline-2,:]

        self.precipitable_water_map[1:met_nline-1,1:met_npixl-1] = total_precipitable_water
        self.precipitable_water_map[:,0] = self.precipitable_water_map[:,met_npixl-2]
        self.precipitable_water_map[:,met_npixl-1] = self.precipitable_water_map[:,1]
        self.precipitable_water_map[0,:] = self.precipitable_water_map[1,:]
        self.precipitable_water_map[met_nline-1,:] = self.precipitable_water_map[met_nline-2,:]


    def trans_ozone(self, ozone_absorption_cross_section, l1b_lat, l1b_lon, l1b_solz, l1b_senz):
        # Deprecated, replaced by compute_ozone_transmittance, which calls gas corrections library
        # Old version retained here just in case it is needed (e.g. for comparison with new version)

        #interpolate ozone map to the L1B grid and compute transmittance
        print('Compute Ozone transmittance ...')
        #npix=len(l1b_lat)
        func=interpolate.RegularGridInterpolator((np.flip(self.oz_lat),self.oz_lon),np.flip(self.ozmap,0))
        ozone_concentration=func(np.array([l1b_lat,l1b_lon]).transpose())
        self.ozone_concentration=ozone_concentration

        start_time = time.perf_counter()

        solar_zenith=np.exp(np.matmul(np.expand_dims(-ozone_concentration/np.cos(np.deg2rad(l1b_solz)), 1),[ozone_absorption_cross_section]))
        sensor_zenith=np.exp(np.matmul(np.expand_dims(-ozone_concentration/np.cos(np.deg2rad(l1b_senz)), 1),[ozone_absorption_cross_section]))

        end_time = time.perf_counter()
        python_time_elapsed = end_time - start_time
        print("Python Ozone time: ", python_time_elapsed)

        return solar_zenith, sensor_zenith
    
    
    def compute_ozone_transmittance(self, ozone_absorption_cross_section, l1b_lat, l1b_lon, l1b_solz, l1b_senz):
        print('Compute Ozone transmittance ...')

        cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
        cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))

        l1_data = gas_corrections.L1_Data()
        l1_data.cos_solar_zenith = cos_solar_zenith
        l1_data.cos_sensor_zenith = cos_sensor_zenith
        l1_data.num_pixels = len(cos_solar_zenith)
        l1_data.num_wavelengths = len(ozone_absorption_cross_section)

        # Interpolate ozone map to the L1B grid
        func = interpolate.RegularGridInterpolator((np.flip(self.oz_lat), self.oz_lon), np.flip(self.ozmap, 0))
        ozone_concentration = func(np.array([l1b_lat, l1b_lon]).transpose())
        self.ozone_concentration = ozone_concentration

        ancillary_data = gas_corrections.Ancillary_Data()
        ancillary_data.ozone_absorption_cross_section = ozone_absorption_cross_section
        ancillary_data.ozone_concentration = ozone_concentration

        # start_time = time.perf_counter()
        gas_transmittances = gas_corrections.ozone_transmittance(l1_data=l1_data, \
                                                                 ancillary_data=ancillary_data)
        # end_time = time.perf_counter()
        # cpp_time_elapsed = end_time - start_time
        # print("C++ Ozone time: ", cpp_time_elapsed)

        return gas_transmittances.solar_zenith, gas_transmittances.sensor_zenith
    

    def compute_co2_transmittance(self, l1b_solz, l1b_senz, sensor_wavelengths):
        print('Compute CO2 transmittance ...')

        l1_data = gas_corrections.L1_Data()
        l1_data.cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
        l1_data.cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))
        l1_data.num_pixels = len(l1_data.cos_solar_zenith)
        l1_data.wavelengths = sensor_wavelengths
        l1_data.num_wavelengths = len(l1_data.wavelengths)

        gas_transmittance_table = self.gas_transmittance_manager.read_gas_transmittance_table(self.gas_transmittance_filepath)

        # start_time = time.perf_counter()
        gas_transmittances = gas_corrections.co2_transmittance(l1_data=l1_data, \
                                                               gas_transmittance_table=gas_transmittance_table)
        # end_time = time.perf_counter()
        # cpp_time_elapsed = end_time - start_time
        # print("C++ CO2 time: ", cpp_time_elapsed)

        return gas_transmittances.solar_zenith, gas_transmittances.sensor_zenith
    

    def compute_co_transmittance(self, l1b_solz, l1b_senz, sensor_wavelengths):
        print('Compute CO transmittance ...')

        l1_data = gas_corrections.L1_Data()
        l1_data.cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
        l1_data.cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))
        l1_data.num_pixels = len(l1_data.cos_solar_zenith)
        l1_data.wavelengths = sensor_wavelengths
        l1_data.num_wavelengths = len(l1_data.wavelengths)

        gas_transmittance_table = self.gas_transmittance_manager.read_gas_transmittance_table(self.gas_transmittance_filepath)

        # start_time = time.perf_counter()
        gas_transmittances = gas_corrections.co_transmittance(l1_data=l1_data, \
                                                                             gas_transmittance_table=gas_transmittance_table)
        # end_time = time.perf_counter()
        # cpp_time_elapsed = end_time - start_time
        # print("C++ CO time: ", cpp_time_elapsed)

        return gas_transmittances.solar_zenith, gas_transmittances.sensor_zenith
    

    def compute_ch4_transmittance(self, l1b_solz, l1b_senz, sensor_wavelengths):
        print('Compute CH4 transmittance ...')

        l1_data = gas_corrections.L1_Data()
        l1_data.cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
        l1_data.cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))
        l1_data.num_pixels = len(l1_data.cos_solar_zenith)
        l1_data.wavelengths = sensor_wavelengths
        l1_data.num_wavelengths = len(l1_data.wavelengths)

        gas_transmittance_table = self.gas_transmittance_manager.read_gas_transmittance_table(self.gas_transmittance_filepath)

        # start_time = time.perf_counter()
        gas_transmittances = gas_corrections.ch4_transmittance(l1_data=l1_data, \
                                                                              gas_transmittance_table=gas_transmittance_table)
        # end_time = time.perf_counter()
        # cpp_time_elapsed = end_time - start_time
        # print("C++ CH4 time: ", cpp_time_elapsed)

        return gas_transmittances.solar_zenith, gas_transmittances.sensor_zenith
    
    
    def compute_n2o_transmittance(self, l1b_solz, l1b_senz, sensor_wavelengths):
        print('Compute N2O transmittance ...')

        l1_data = gas_corrections.L1_Data()
        l1_data.cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
        l1_data.cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))
        l1_data.num_pixels = len(l1_data.cos_solar_zenith)
        l1_data.wavelengths = sensor_wavelengths
        l1_data.num_wavelengths = len(l1_data.wavelengths)

        gas_transmittance_table = self.gas_transmittance_manager.read_gas_transmittance_table(self.gas_transmittance_filepath)

        # start_time = time.perf_counter()
        gas_transmittances = gas_corrections.n2o_transmittance(l1_data=l1_data, \
                                                                              gas_transmittance_table=gas_transmittance_table)
        # end_time = time.perf_counter()
        # cpp_time_elapsed = end_time - start_time
        # print("C++ N2O time: ", cpp_time_elapsed)

        return gas_transmittances.solar_zenith, gas_transmittances.sensor_zenith
    
        
    def compute_o2_transmittance(self, l1b_solz, l1b_senz, l1b_reflectance, sensor_wavelengths):
        print('Compute O2 transmittance ...')

        l1_data = gas_corrections.L1_Data()
        l1_data.cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))
        l1_data.cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
        l1_data.num_pixels = len(l1_data.cos_solar_zenith)
        l1_data.reflectance = l1b_reflectance
        l1_data.wavelengths = sensor_wavelengths
        l1_data.num_wavelengths = len(sensor_wavelengths)

        gas_transmittance_table = self.gas_transmittance_manager.read_gas_transmittance_table(self.gas_transmittance_filepath)

        # start_time = time.perf_counter()
        gas_transmittances = gas_corrections.o2_transmittance(l1_data=l1_data, \
                                                                             gas_transmittance_table=gas_transmittance_table, \
                                                                             oxygen_A_band_option=gas_corrections.Oxygen_A_Band_Option().TRANSMITTANCE_TABLE)
        # end_time = time.perf_counter()
        # cpp_time_elapsed = end_time - start_time
        # print("C++ O2 time: ", cpp_time_elapsed)

        return gas_transmittances.solar_zenith, gas_transmittances.sensor_zenith
    
    
    def trans_no2(self, no2_absorption_cross_section, month, l1b_lat, l1b_lon, l1b_solz, l1b_senz):
        # Deprecated, replaced by compute_no2_transmittance, which calls gas corrections library
        # Old version retained here just in case it is needed (e.g. for comparison with new version)

        #interpolate no2 map to the L1B grid and compute transmittance
        print('Compute NO2 transmittance ...')
        #npix=len(l1b_lat)
        a_285 = no2_absorption_cross_section * (1.0 - 0.003*(285.0-294.0))
        a_225 = no2_absorption_cross_section * (1.0 - 0.003*(225.0-294.0))
        func=interpolate.RegularGridInterpolator((np.flip(self.no2_frac_lat),self.no2_frac_lon),np.flip(self.no2_frac,0))
        fraction_tropospheric_no2_above_200m=func(np.array([l1b_lat,l1b_lon]).transpose())
        no2_strat=self.no2_strat[int(month)-1,:,:]
        func=interpolate.RegularGridInterpolator((np.flip(self.no2_lat),self.no2_lon),np.flip(no2_strat,0))
        stratospheric_no2_concentration=func(np.array([l1b_lat,l1b_lon]).transpose())
        no2_tropo=self.no2_tropo[int(month)-1,:,:]
        func=interpolate.RegularGridInterpolator((np.flip(self.no2_lat),self.no2_lon),np.flip(no2_tropo,0))
        tropospheric_no2_concentration=func(np.array([l1b_lat,l1b_lon]).transpose())

        start_time = time.perf_counter()

        l1b_no2_trop200=fraction_tropospheric_no2_above_200m*tropospheric_no2_concentration
        l1b_no2_trop200[l1b_no2_trop200<0]=0.
        solar_zenith=np.exp(-(np.matmul(np.expand_dims(l1b_no2_trop200/np.cos(np.deg2rad(l1b_solz)), 1),[a_285])+np.matmul(np.expand_dims(stratospheric_no2_concentration/np.cos(np.deg2rad(l1b_solz)),1),[a_225])))
        sensor_zenith=np.exp(-(np.matmul(np.expand_dims(l1b_no2_trop200/np.cos(np.deg2rad(l1b_senz)), 1),[a_285])+np.matmul(np.expand_dims(stratospheric_no2_concentration/np.cos(np.deg2rad(l1b_senz)),1),[a_225])))
        
        end_time = time.perf_counter()
        python_time_elapsed = end_time - start_time
        print("Python NO2 time: ", python_time_elapsed)

        self.fraction_tropospheric_no2_above_200m=fraction_tropospheric_no2_above_200m
        self.stratospheric_no2_concentration=stratospheric_no2_concentration
        self.tropospheric_no2_concentration=tropospheric_no2_concentration
        return solar_zenith, sensor_zenith
    
    
    def compute_no2_transmittance(self, no2_absorption_cross_section, month, l1b_lat, l1b_lon, l1b_solz, l1b_senz):
        print('Compute NO2 transmittance ...')

        cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
        cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))

        l1_data = gas_corrections.L1_Data()
        l1_data.cos_solar_zenith = cos_solar_zenith
        l1_data.cos_sensor_zenith = cos_sensor_zenith
        l1_data.num_pixels = len(cos_solar_zenith)
        l1_data.num_wavelengths = len(no2_absorption_cross_section)
        
        # Interpolate no2 map to the L1B grid
        func = interpolate.RegularGridInterpolator((np.flip(self.no2_frac_lat), self.no2_frac_lon), np.flip(self.no2_frac, 0))
        fraction_tropospheric_no2_above_200m = func(np.array([l1b_lat, l1b_lon]).transpose())

        no2_strat = self.no2_strat[int(month)-1, :, :]
        func = interpolate.RegularGridInterpolator((np.flip(self.no2_lat), self.no2_lon), np.flip(no2_strat, 0))
        stratospheric_no2_concentration = func(np.array([l1b_lat, l1b_lon]).transpose())

        no2_tropo = self.no2_tropo[int(month)-1, :, :]
        func = interpolate.RegularGridInterpolator((np.flip(self.no2_lat), self.no2_lon), np.flip(no2_tropo, 0))
        tropospheric_no2_concentration = func(np.array([l1b_lat, l1b_lon]).transpose())

        ancillary_data = gas_corrections.Ancillary_Data()
        ancillary_data.no2_absorption_cross_section = no2_absorption_cross_section
        ancillary_data.fraction_tropospheric_no2_above_200m = fraction_tropospheric_no2_above_200m
        ancillary_data.tropospheric_no2_concentration = tropospheric_no2_concentration
        ancillary_data.stratospheric_no2_concentration = stratospheric_no2_concentration

        # start_time = time.perf_counter()
        gas_transmittances = gas_corrections.no2_transmittance(l1_data=l1_data, \
                                                                              ancillary_data=ancillary_data)
        # end_time = time.perf_counter()
        # cpp_time_elapsed = end_time - start_time
        # print("C++ NO2 time: ", cpp_time_elapsed)        

        self.fraction_tropospheric_no2_above_200m=fraction_tropospheric_no2_above_200m
        self.stratospheric_no2_concentration=stratospheric_no2_concentration
        self.tropospheric_no2_concentration=tropospheric_no2_concentration

        return gas_transmittances.solar_zenith, gas_transmittances.sensor_zenith
    

    def compute_h2o_transmittance(self, l1b_solz, l1b_senz, l1b_reflectance, sensor_wavelengths):
        print('Compute Water Vapor transmittance ...')

        l1_data = gas_corrections.L1_Data()
        l1_data.cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))
        l1_data.cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
        l1_data.num_pixels = len(l1_data.cos_solar_zenith)
        l1_data.reflectance = l1b_reflectance
        l1_data.wavelengths = sensor_wavelengths
        l1_data.num_wavelengths = len(sensor_wavelengths)

        gas_transmittance_table = self.gas_transmittance_manager.read_gas_transmittance_table(self.gas_transmittance_filepath)

        # start_time = time.perf_counter()
        gas_transmittances = gas_corrections.h2o_transmittance(l1_data=l1_data, \
                                                               gas_transmittance_table=gas_transmittance_table)
        # end_time = time.perf_counter()
        # cpp_time_elapsed = end_time - start_time
        # print("C++ H2O time: ", cpp_time_elapsed)

        return gas_transmittances.solar_zenith, gas_transmittances.sensor_zenith
    
    
    def trans_o2_aer(self,l1b_solz, l1b_senz):
        #O2 transmittance for aerosols
        ao2=np.array([-1.0796, 9.0481e-2, -6.8452e-3])
        airmass=1/np.cos(np.deg2rad(l1b_solz))+1/np.cos(np.deg2rad(l1b_senz))
        t_o2=1.0+np.power(10,ao2[0]+ao2[1]*airmass+ao2[2]*airmass**2)
        return t_o2
           
    def trans_o2_ray(self,l1b_solz, l1b_senz):
        #O2 transmittance for Rayleigh
        ao2=np.array([-1.3491, 0.1155, -7.0218e-3])
        airmass=1/np.cos(np.deg2rad(l1b_solz))+1/np.cos(np.deg2rad(l1b_senz))
        ray_o2=1.0/(1.0+np.power(10,ao2[0]+ao2[1]*airmass+ao2[2]*airmass**2))
        return ray_o2     
                
    def get_metdata(self, l1b_lat, l1b_lon):
        #interpolate pressure, relative humidity (RH) and windspeed to the L1B grid
        func=interpolate.RegularGridInterpolator((np.flip(self.met_lat),self.met_lon),np.flip(self.pressmap,0))
        l1b_press=func(np.array([l1b_lat,l1b_lon]).transpose())

        func=interpolate.RegularGridInterpolator((np.flip(self.met_lat),self.met_lon),np.flip(self.rhmap,0))
        l1b_rh=func(np.array([l1b_lat,l1b_lon]).transpose())

        func=interpolate.RegularGridInterpolator((np.flip(self.met_lat),self.met_lon),np.flip(self.wsmap,0))
        l1b_ws=func(np.array([l1b_lat,l1b_lon]).transpose())
        
        return l1b_press, l1b_rh, l1b_ws
    
    def whitecaps(self, band, l1b_solz, l1b_senz, l1b_ws, l1b_pressure, taur):
        print('Compute whitecaps reflectance ...')
        awc_band=np.array([380,412,443,490,510,555,670,765,865,1000,5000])
        awc_tab=np.array([1.0,1.0,1.0,1.0,1.0,1.0,0.889225,0.760046,0.644950,0.0,0.0])
        wc_ws_min=6.33
        wc_ws_max=12
        p0=1013.25
        npix=len(l1b_ws)
        nband=len(band)
        l1b_tlf=np.zeros([npix,nband])
        func=interpolate.interp1d(awc_band,awc_tab,kind='linear')
        awc=func([band])
        idx1 = l1b_ws <= wc_ws_max 
        idx2 = l1b_ws >= wc_ws_min
        idx  = idx1 & idx2
        x=np.matmul(np.expand_dims(1.925e-5*np.power(l1b_ws[idx]-wc_ws_min,3),1),[awc])
        l1b_tlf[idx,:]=x
        solar_zenith=np.exp(-0.5*np.matmul(np.expand_dims(l1b_pressure/p0/np.cos(np.deg2rad(l1b_solz)),1),[taur]))
        sensor_zenith=np.exp(-0.5*np.matmul(np.expand_dims(l1b_pressure/p0/np.cos(np.deg2rad(l1b_senz)),1),[taur]))
        x1=np.zeros([npix,nband])
        for i in range(nband):
            x1[:,i]=np.cos(np.deg2rad(l1b_solz))
        l1b_tlf=l1b_tlf / np.pi * x1 * solar_zenith * sensor_zenith
        return l1b_tlf
    
     











     