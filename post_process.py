# -*- coding: utf-8 -*-
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import datetime
import dask
import domain
import glob
import itertools
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.animation
from matplotlib.colors import LogNorm
import matplotlib.colors as col
import math
import multiprocessing
import netCDF4
import numpy
from operator import itemgetter
import os
import pyproj
import pathlib
import pandas
import re
from scipy import stats
import sqlalchemy
from sklearn import metrics
import sys
import shutil
from sklearn.linear_model import LinearRegression
import stat
import time
from warnings import simplefilter
import wrf
import xarray

simplefilter(action='ignore', category=FutureWarning)

def calc_kernel(samp):
    return stats.gaussian_kde(samp)(samp)
airdb_engine = sqlalchemy.create_engine("dialect+driver://username:password@host:port/database")  #observation station set1
airdb_engine2 = sqlalchemy.create_engine("dialect+driver://username:password@host:port/database") #observation station set2


class Post_Process():
    def __init__(self,data_type,date_start,date_end,root_dir,source_dirs,extracted_layers,sites_includedprovinces,sites_includedcities):
        self.data_type = data_type
        self.date_start = datetime.datetime.strptime(date_start, "%Y-%m-%d")
        self.date_end = datetime.datetime.strptime(date_end, "%Y-%m-%d") + datetime.timedelta(hours=23)
        self.date_cover = pandas.date_range(self.date_start + datetime.timedelta(hours=8), self.date_end + datetime.timedelta(hours=8), freq='H')
        self.root_dir = root_dir
        self.source_dirs = source_dirs
        self.extracted_layers = extracted_layers
        self.sites_includedprovinces = sites_includedprovinces
        self.sites_includedcities = sites_includedcities
        self.heiti_font = matplotlib.font_manager.FontProperties\
            (fname=root_dir + "/utils/simhei.ttf", size=25)
        if self.data_type == 'WRF':
            self.factors = ['TC|°(C)|-10~40|ambient temperature','WS|m/s|0~10|ambient wind speed',\
        'WD|degree|0~360|ambient wind direction','Pressure|hPa|900~1020|ambient pressure',\
        'RH|%|40~100|ambient relative humidity','PBLH|m|0~1000|planetary boundary layer']

            if self.sites_includedprovinces != ['all']:
                where_str = 'station_code is not null and station_province in ' \
                    + str(self.sites_includedprovinces).replace('[','(').replace(']',')')
            else:
                where_str = 'station_code is not null'
            if self.sites_includedcities != ['all']:
                where_str = where_str + ' and city_name in ' \
                    + str(self.sites_includedcities).replace('[','(').replace(']',')')
            ds_station = pandas.read_sql_query('SELECT * FROM t_weather_data_station WHERE ' \
                + where_str + ' ORDER BY station_code', airdb_engine)
            self.siteIDs = ds_station['station_code'].to_numpy().astype(str)
            self.sitenames = ds_station['station_name'].to_numpy()
            self.sitelats = ds_station['station_lat'].to_numpy()
            self.sitelons = ds_station['station_lon'].to_numpy()
            self.sitecities = ds_station['city_name'].to_numpy()
            self.siteprovince = ds_station['station_province'].to_numpy()
            self.sitecols = []
            self.siterows = []
            # 筛选网格空间范围内的站点
            for source_index in range(len(self.source_dirs.split(';'))):
                nc_ds = netCDF4.Dataset(self.source_dirs.split(';')[source_index] + '/wrfout_d01_' + (self.date_start + datetime.timedelta(days = -1)).strftime('%Y-%m-%d_12:00:00'))
                Domain_NCOLS = int(nc_ds.getncattr('WEST-EAST_PATCH_END_UNSTAG')) # len(nc_ds.dimensions['west_east']) also works
                Domain_NROWS = int(nc_ds.getncattr('SOUTH-NORTH_PATCH_END_UNSTAG')) # len(nc_ds.dimensions['south_north']) also works
                Domain_CEN_LAT = float(nc_ds.getncattr('CEN_LAT'))
                Domain_CEN_LON = float(nc_ds.getncattr('CEN_LON'))
                Domain_TRUELAT1 = float(nc_ds.getncattr('TRUELAT1'))
                Domain_TRUELAT2 = float(nc_ds.getncattr('TRUELAT2'))
                Domain_XORIG = float(wrf.cartopy_xlim(wrfin=nc_ds)[0])
                Domain_YORIG = float(wrf.cartopy_ylim(wrfin=nc_ds)[0])
                Domain_DX = float(nc_ds.getncattr('DX'))
                Domain_DY = float(nc_ds.getncattr('DY'))
                LCCProj = pyproj.Proj(proj='lcc', lat_1=Domain_TRUELAT1,lat_2=Domain_TRUELAT2,lat_0=Domain_CEN_LAT,lon_0=Domain_CEN_LON,a=6370000,b=6370000)
                lcc_x, lcc_y = LCCProj(self.sitelons, self.sitelats)
                self.sitecols.append(numpy.trunc((lcc_x - Domain_XORIG) / Domain_DX).astype(int))
                self.siterows.append(numpy.trunc((lcc_y - Domain_YORIG) / Domain_DY).astype(int))
                condition = (self.sitecols[source_index] >= 0) & (self.sitecols[source_index] < Domain_NCOLS) \
                    & (self.siterows[source_index] >= 0) & (self.siterows[source_index] < Domain_NROWS)
                selected_indices = numpy.where(condition)[0]
                self.siteIDs = self.siteIDs[selected_indices]
                self.sitenames = self.sitenames[selected_indices]
                self.sitelats = self.sitelats[selected_indices]
                self.sitelons = self.sitelons[selected_indices]
                self.sitecities = self.sitecities[selected_indices]
                self.siteprovince = self.siteprovince[selected_indices]
                self.sitecols[source_index] = self.sitecols[source_index][selected_indices]
                self.siterows[source_index] = self.siterows[source_index][selected_indices]
             
        if self.data_type == 'CMAQ_regularsites':  # 针对全国空气质量国控点
        #     self.factors = ['CO|mg/m3|0.1~2|carbon monoxide','NO2|μg/m3|0~60|nitrogen dioxide','O3|μg/m3|60~260|Ozone',
        # 'SO2|μg/m3|0~60|sulfur dioxide','PM25|μg/m3|0~200|fine particulate matter','PM10|μg/m3|30~200|inhalable particulate matter']
            self.factors = ['NO2|μg/m3|0~60|nitrogen dioxide']
            #提取观测站点信息
            if self.sites_includedprovinces != ['all']:
                where_str = 'station_code is not null and station_province in ' \
                    + str(self.sites_includedprovinces).replace('[','(').replace(']',')')
            else:
                where_str = 'station_code is not null'
            if self.sites_includedcities != ['all']:
                where_str = where_str + ' and city_name in ' \
                    + str(self.sites_includedcities).replace('[','(').replace(']',')')
            ds_station = pandas.read_sql_query('SELECT * FROM t_na_stations WHERE ' \
                + where_str + ' ORDER BY station_code', airdb_engine)
            self.siteIDs = ds_station['station_code'].to_numpy().astype(str) 
            self.sitenames = ds_station['station_name'].to_numpy()
            self.sitelats = ds_station['station_lat'].to_numpy()
            self.sitelons = ds_station['station_lon'].to_numpy()
            self.sitecities = ds_station['city_name'].to_numpy()
            self.siteprovince = ds_station['station_province'].to_numpy()
            self.sitecols = []
            self.siterows = []
            # 筛选网格空间范围内的站点
            for source_index in range(len(self.source_dirs.split(';'))):
                init_file = self.source_dirs.split(';')[source_index] + '/ACONC_' + self.date_start.strftime('%Y-%m-%d') + '.nc'            
                if not os.path.exists(init_file):
                    init_file = self.source_dirs.split(';')[source_index] + '/CCTM_ACONC_' + self.date_start.strftime('%Y-%m-%d') + '.nc'
                if not os.path.exists(init_file):
                    print('the error of data sources!!')
                    sys.exit()
                else:
                    nc_ds = netCDF4.Dataset(init_file)
                Domain_CEN_LAT = float(nc_ds.getncattr('YCENT'))
                Domain_CEN_LON = float(nc_ds.getncattr('XCENT'))
                Domain_TRUELAT1 = float(nc_ds.getncattr('P_ALP'))
                Domain_TRUELAT2 = float(nc_ds.getncattr('P_BET'))
                Domain_XORIG = float(nc_ds.getncattr('XORIG'))
                Domain_YORIG = float(nc_ds.getncattr('YORIG'))
                Domain_DX = float(nc_ds.getncattr('XCELL'))
                Domain_DY = float(nc_ds.getncattr('YCELL'))
                Domain_NCOLS = int(nc_ds.getncattr('NCOLS'))
                Domain_NROWS = int(nc_ds.getncattr('NROWS'))
                LCCProj = pyproj.Proj(proj='lcc', lat_1=Domain_TRUELAT1,lat_2=Domain_TRUELAT2,lat_0=Domain_CEN_LAT,lon_0=Domain_CEN_LON,a=6370000,b=6370000)
                lcc_x, lcc_y = LCCProj(self.sitelons, self.sitelats)
                self.sitecols.append(numpy.trunc((lcc_x - Domain_XORIG) / Domain_DX).astype(int))
                self.siterows.append(numpy.trunc((lcc_y - Domain_YORIG) / Domain_DY).astype(int))
                condition = (self.sitecols[source_index] >= 0) & (self.sitecols[source_index] < Domain_NCOLS) \
                    & (self.siterows[source_index] >= 0) & (self.siterows[source_index] < Domain_NROWS)
                selected_indices = numpy.where(condition)[0]
                self.siteIDs = self.siteIDs[selected_indices]
                self.sitenames = self.sitenames[selected_indices]
                self.sitelats = self.sitelats[selected_indices]
                self.sitelons = self.sitelons[selected_indices]
                self.sitecities = self.sitecities[selected_indices]
                self.siteprovince = self.siteprovince[selected_indices]
                self.sitecols[source_index] = self.sitecols[source_index][selected_indices]
                self.siterows[source_index] = self.siterows[source_index][selected_indices]

        if self.data_type == 'CMAQ_supersites':
            self.factors = ['CO|mg/m3|0.1~2|carbon monoxide','NO2|μg/m3|0~60|nitrogen dioxide','O3|μg/m3|60~260|Ozone',
        'SO2|μg/m3|0~60|sulfur dioxide','NO|μg/m3|0~60|nitric oxide','NH3|μg/m3|0~10|ammonia','HNO3|μg/m3|0~30|nitric acid',
        'PM25|μg/m3|0~200|fine particulate matter','PM10|μg/m3|30~200|inhalable particulate matter',
        'PM25_NH4|μg/m3|0~30|ammonium ion in fine particles','PM25_NO3|μg/m3|0~30|nitrate ion in fine particles',
        'PM25_SO4|μg/m3|0~30|sulfate ion in fine particles','PM25_OC|μg/m3|0~50|organic carbon in fine particles',
        'PM25_EC|μg/m3|0~10|elemental carbon in fine particles','TS|μg/m3|0~30|Total Sulfur',\
        'TNN|μg/m3|0~100|Total Nitrate Nitrogen','TAN|μg/m3|0~30|Total Ammonia Nitrogen']
        #'PM25_Cl|μg/m3|0~5|chlorine ion in fine particles','PM25_SOIL|μg/m3|0~10|Soil matter in fine particles'\
        #'PM25_Na|μg/m3|0~5|sodium ion in fine particles','PM25_Mg|μg/m3|0~5|magnesium ion in fine particles',\
        #'PM25_K|μg/m3|0~5|potassium ion in fine particles','PM25_Ca|μg/m3|0~5|calcium ion in fine particles',\
        #'PM25_POC|μg/m3|0~30|primary organic carbon in fine particles','PM25_POM|μg/m3|0~30|primary organic matter in fine particles',\
        #'PM25_SOC|μg/m3|0~30|secondary organic carbon in fine particles','PM25_SOM|μg/m3|0~30|secondary organic matter in fine particles'
            #提取观测站点信息
            if self.sites_includedprovinces != ['all']:
                where_str = 'station_code is not null and station_province in ' \
                    + str(self.sites_includedprovinces).replace('[','(').replace(']',')')
            else:
                where_str = 'station_code is not null'
            if self.sites_includedcities != ['all']:
                where_str = where_str + ' and station_city in ' \
                    + str(self.sites_includedcities).replace('[','(').replace(']',')')
            ds_station = pandas.read_sql_query('SELECT * FROM Supersite_Sites WHERE ' \
                    + where_str + ' ORDER BY station_code', airdb_engine2)
            self.siteIDs = ds_station['station_code'].to_numpy().astype(str)
            self.sitenames = ds_station['station_name'].to_numpy()
            self.sitelats = ds_station['station_lat'].to_numpy()
            self.sitelons = ds_station['station_lon'].to_numpy()
            self.sitecities = ds_station['station_city'].to_numpy()
            self.siteprovince = ds_station['station_province'].to_numpy()
            self.sitecols = []
            self.siterows = []
            # 筛选网格空间范围内的站点
            for source_index in range(len(self.source_dirs.split(';'))):
                init_file = self.source_dirs.split(';')[source_index] + '/ACONC_' + self.date_start.strftime('%Y-%m-%d') + '.nc'            
                if not os.path.exists(init_file):
                    init_file = self.source_dirs.split(';')[source_index] + '/CCTM_ACONC_' + self.date_start.strftime('%Y-%m-%d') + '.nc'
                if not os.path.exists(init_file):
                    print('the error of data sources!!')
                    sys.exit()
                else:
                    nc_ds = netCDF4.Dataset(init_file)
                Domain_CEN_LAT = float(nc_ds.getncattr('YCENT'))
                Domain_CEN_LON = float(nc_ds.getncattr('XCENT'))
                Domain_TRUELAT1 = float(nc_ds.getncattr('P_ALP'))
                Domain_TRUELAT2 = float(nc_ds.getncattr('P_BET'))
                Domain_XORIG = float(nc_ds.getncattr('XORIG'))
                Domain_YORIG = float(nc_ds.getncattr('YORIG'))
                Domain_DX = float(nc_ds.getncattr('XCELL'))
                Domain_DY = float(nc_ds.getncattr('YCELL'))
                Domain_NCOLS = int(nc_ds.getncattr('NCOLS'))
                Domain_NROWS = int(nc_ds.getncattr('NROWS'))
                LCCProj = pyproj.Proj(proj='lcc', lat_1=Domain_TRUELAT1,lat_2=Domain_TRUELAT2,lat_0=Domain_CEN_LAT,lon_0=Domain_CEN_LON,a=6370000,b=6370000)
                lcc_x, lcc_y = LCCProj(self.sitelons, self.sitelats)
                self.sitecols.append(numpy.trunc((lcc_x - Domain_XORIG) / Domain_DX).astype(int))
                self.siterows.append(numpy.trunc((lcc_y - Domain_YORIG) / Domain_DY).astype(int))
                condition = (self.sitecols[source_index] >= 0) & (self.sitecols[source_index] < Domain_NCOLS) \
                    & (self.siterows[source_index] >= 0) & (self.siterows[source_index] < Domain_NROWS)
                selected_indices = numpy.where(condition)[0]
                self.siteIDs = self.siteIDs[selected_indices]
                self.sitenames = self.sitenames[selected_indices]
                self.sitelats = self.sitelats[selected_indices]
                self.sitelons = self.sitelons[selected_indices]
                self.sitecities = self.sitecities[selected_indices]
                self.siteprovince = self.siteprovince[selected_indices]
                self.sitecols[source_index] = self.sitecols[source_index][selected_indices]
                self.siterows[source_index] = self.siterows[source_index][selected_indices]


        if self.data_type == 'TROPOMI':
            self.start_day = datetime.datetime.strptime(date_start, "%Y-%m-%d") #开始日期，日期选择
            self.end_day = datetime.datetime.strptime(date_end, "%Y-%m-%d") #结束日期，日期选择
            #self.factors = ['NO2|mol/m2|0~0.0192|nitrogen_dioxide_total_column|v221']
            self.factors = ['NO2|mol/m2|0~30|nitrogen_dioxide_total_column|v221']
            #'CO|mol/m2|0~0.07|carbonmonoxide_total_column', 'O3|mol/m2|0~0.3|ozone_total_column','CH4|mol/m2|0~100|methane_total_column']
            self.mpijob_partition = "64c512g"
            self.serialjob_partition = "64c512g"
            self.mpijob_npn = 64
            self.domain = domain.domain_definition('China_12km',35.0,105,25.0,45.0,512,512,12000,12000,-3072000.0,-3072000.0)
            #self.domain = domain.domain_definition("CHINA_12km","lambert",35.0,105,258,258,25.0,45.0,512,512,12000,12000,-3072000.0,-3072000.0,45,5000)
            self.source_path = self.source_dirs
            # self.savepath_MP4 = self.output_dir + '/MP4'  # 图片输出文件夹
            # self.single_PNG = self.output_dir + '/MP4_PNG'
            #仅用于保持一致性
            nc_ds = netCDF4.Dataset(self.source_dirs[1]+'/ACONC_2019-01-31.nc')
            self.Domain_NCOLS = int(nc_ds.getncattr('NCOLS'))
            self.Domain_NROWS = int(nc_ds.getncattr('NROWS'))
            self.Domain_CEN_LAT = float(nc_ds.getncattr('YCENT'))
            self.Domain_CEN_LON = float(nc_ds.getncattr('XCENT'))
            self.Domain_TRUELAT1 = float(nc_ds.getncattr('P_ALP'))
            self.Domain_TRUELAT2 = float(nc_ds.getncattr('P_BET'))
            self.Domain_XORIG = float(nc_ds.getncattr('XORIG'))
            self.Domain_YORIG = float(nc_ds.getncattr('YORIG'))
            self.Domain_DX = float(nc_ds.getncattr('XCELL'))
            self.Domain_DY = float(nc_ds.getncattr('YCELL'))
        
        if self.data_type == 'adjoint_sensitivity':
            self.factors = ['CO|CF/(mol/s)|0.000005~0.001|carbon monoxide',\
                'SO2|CF/(mol/s)|0.000005~0.001|SO2', 'NO|CF/(mol/s)|0.000005~0.001|NO', 'NO2|CF/(mol/s)|0.000005~0.001|NO2',\
                'NH3|CF/(mol/s)|0.000005~0.001|NH3', 'VOCs|CF/(mol/s)|0.000005~0.001|VOCs', 'PNO3|CF/(g/s)|0.000005~0.001|primary NO3 emission',\
                'PMC|CF/(g/s)|0.000005~0.001|primary coarse emission','PSO4|CF/(g/s)|0.000005~0.001|primary SO4 emission',\
                'PEC|CF/(g/s)|0.000005~0.001|primary EC emission', 'POC|CF/(g/s)|0.000005~0.001|primary POC emission',\
                'PMFINE|CF/(g/s)|0.000005~0.001|primary ohter fine particles emission'] 
                #VOCs = MEOH+OLE+TERP+TOL+ETHA+ETOH+ETH+FORM+ALDX+IOLE+PAR+XYL+ISOP+ALD2
            nc_ds = netCDF4.Dataset(self.source_dirs.split(';')[0] + '/lgrid_em_' + self.date_start.strftime('%Y-%m-%d') + '.nc')
            self.Domain_NCOLS = int(nc_ds.getncattr('NCOLS'))
            self.Domain_NROWS = int(nc_ds.getncattr('NROWS'))
            self.Domain_CEN_LAT = float(nc_ds.getncattr('YCENT'))
            self.Domain_CEN_LON = float(nc_ds.getncattr('XCENT'))
            self.Domain_TRUELAT1 = float(nc_ds.getncattr('P_ALP'))
            self.Domain_TRUELAT2 = float(nc_ds.getncattr('P_BET'))
            self.Domain_XORIG = float(nc_ds.getncattr('XORIG'))
            self.Domain_YORIG = float(nc_ds.getncattr('YORIG'))
            self.Domain_DX = float(nc_ds.getncattr('XCELL'))
            self.Domain_DY = float(nc_ds.getncattr('YCELL'))

        if self.data_type == 'adjoint_semi_normalized_sensitivity':
            self.factors = ['CO|CF|0.000005~0.001|carbon monoxide',\
                'SO2|CF|0.000005~0.001|SO2', 'NO|CF|0.000005~0.001|NO', 'NO2|CF|0.000005~0.001|NO2',\
                'NH3|CF|0.000005~0.001|NH3', 'VOCs|CF|0.000005~0.001|VOCs', 'PNO3|CF|0.000005~0.001|primary NO3 emission',\
                'PMC|CF|0.000005~0.001|primary coarse emission','PSO4|CF|0.000005~0.001|primary SO4 emission',\
                'PEC|CF|0.000005~0.001|primary EC emission', 'POC|CF|0.000005~0.001|primary POC emission',\
                'PMFINE|CF|0.000005~0.001|primary ohter fine particles emission'] 
                #VOCs = MEOH+OLE+TERP+TOL+ETHA+ETOH+ETH+FORM+ALDX+IOLE+PAR+XYL+ISOP+ALD2
            nc_ds = netCDF4.Dataset(self.source_dirs.split(';')[0] + '/lgrid_em_' + self.date_start.strftime('%Y-%m-%d') + '.nc')
            self.Domain_NCOLS = int(nc_ds.getncattr('NCOLS'))
            self.Domain_NROWS = int(nc_ds.getncattr('NROWS'))
            self.Domain_CEN_LAT = float(nc_ds.getncattr('YCENT'))
            self.Domain_CEN_LON = float(nc_ds.getncattr('XCENT'))
            self.Domain_TRUELAT1 = float(nc_ds.getncattr('P_ALP'))
            self.Domain_TRUELAT2 = float(nc_ds.getncattr('P_BET'))
            self.Domain_XORIG = float(nc_ds.getncattr('XORIG'))
            self.Domain_YORIG = float(nc_ds.getncattr('YORIG'))
            self.Domain_DX = float(nc_ds.getncattr('XCELL'))
            self.Domain_DY = float(nc_ds.getncattr('YCELL'))

        if self.data_type == 'emission':
            self.factors = ['CO','NO2','O3','SO2']
            nc_ds = netCDF4.Dataset(self.source_dirs.split(';')[0] + '/emission_' + self.date_start.strftime('%Y-%m-%d') + '.nc')
            self.Domain_NCOLS = int(nc_ds.getncattr('NCOLS'))
            self.Domain_NROWS = int(nc_ds.getncattr('NROWS'))
            self.Domain_CEN_LAT = float(nc_ds.getncattr('YCENT'))
            self.Domain_CEN_LON = float(nc_ds.getncattr('XCENT'))
            self.Domain_TRUELAT1 = float(nc_ds.getncattr('P_ALP'))
            self.Domain_TRUELAT2 = float(nc_ds.getncattr('P_BET'))
            self.Domain_XORIG = float(nc_ds.getncattr('XORIG'))
            self.Domain_YORIG = float(nc_ds.getncattr('YORIG'))
            self.Domain_DX = float(nc_ds.getncattr('XCELL'))
            self.Domain_DY = float(nc_ds.getncattr('YCELL'))

        matplotlib.rcParams['agg.path.chunksize'] = 10000
        self.simhei_font = matplotlib.font_manager.FontProperties\
            (fname=self.root_dir + "/utils/simhei.ttf")
        self.timesnr_font = matplotlib.font_manager.FontProperties\
            (fname=self.root_dir + "/utils/TimesNewRoman.ttf")

    def makefiles_semi_normalized_sensitivity_hourly(self, bwdDir):
        with open(f'{bwdDir}/BWD.slurm') as f:
            bwdSlurm = f.read()
        fwdDir = re.search("FWDDIR.*", bwdSlurm).group().split("=")[-1].strip() + '/'
        emissionname = f"{fwdDir}/emission.nc"
        emission = netCDF4.Dataset(emissionname, mode='r', open=True)
        varLst = ['PMFINE', 'MEOH', 'POC', 'OLE', 'TERP', 'TOL', 'PEC', 'NO', 'ETHA', 'ETOH', 'ETH', 'NO2', 'NH3', 'FORM', 'PNO3', \
                    'ALDX', 'PMC', 'IOLE', 'PAR', 'XYL', 'PSO4', 'ISOP', 'CO', 'ALD2', 'SO2']
        for day_index in range((self.date_end - self.date_start).days + 1):
            current_day  = self.date_start + datetime.timedelta(days=day_index)
            filename = f"{self.output_dir}/semi_normalized_sensitivity_em_" + current_day.strftime('%Y-%m-%d') + ".nc"
            if os.path.exists(filename): os.system(f"rm {filename}")
            lgridPath = f"{bwdDir}/lgrid_em_" + current_day.strftime('%Y-%m-%d') + ".nc"
            lgrid = netCDF4.Dataset(lgridPath, mode='r', open=True)
            tflagArr = numpy.array(lgrid.variables["TFLAG"][:][1:25, 0, :])
            tflagArr = tflagArr.astype(int)
            conchr = lgrid.variables['PMFINE'][:][:][:][:]
            ltime, llay, lrow, lcol = conchr.shape
            ltime = 24
            # 每层都需要则设置为7，所有高度层累加则设置为1
            llay = 1
            ldattim = 2            
            sensfile = netCDF4.Dataset(filename, mode='w', format='NETCDF3_64BIT')
            attrs=["IOAPI_VERSION", "EXEC_ID", "FTYPE", "CDATE", "CTIME", "WDATE", "WTIME", "SDATE", "STIME", \
                            "TSTEP", "NTHIK", "NCOLS", "NROWS", "GDTYP", "P_ALP", "P_BET", "P_GAM", "XCENT", "YCENT", "XORIG",\
                            "YORIG", "XCELL", "YCELL", "VGTYP", "VGTOP", "VGLVLS", "GDNAM", "HISTORY"]
            for attr in attrs:
                if hasattr(lgrid, attr):
                    attrVal = getattr(lgrid, attr)
                    setattr(sensfile, attr, attrVal)
            setattr(sensfile,"NVARS", len(varLst)) # only one variable
            setattr(sensfile,"FILEDESC","sensitivity file: lgrid_em * emission")

            sensfile.createDimension("TSTEP", None)
            sensfile.createDimension("DATE-TIME", ldattim)
            sensfile.createDimension("LAY", llay)
            sensfile.createDimension("VAR", len(varLst))
            sensfile.createDimension("ROW", lrow)
            sensfile.createDimension("COL", lcol)
            sensfile.sync()
            dattim = numpy.zeros([ltime,len(varLst),ldattim])
            sens_tflag = sensfile.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
            for j in range(len(varLst)):
                dattim[:, j, :] = tflagArr
            sens_tflag[:] = dattim
            varattrs=["long_name","units","var_desc"]
            for varattr in varattrs:
                if hasattr(lgrid.variables['TFLAG'], varattr): 
                    varattrVal = getattr(lgrid.variables['TFLAG'], varattr)
                    setattr(sens_tflag, varattr, varattrVal)
            for var in varLst:
                x = int((current_day - self.date_start).days)
                emissionArr = emission.variables[var][:][:][:][:][x*24+1 : x*24+25,0:7,:,:]
                lgridArr = lgrid.variables[var][:][:][:][:][1:25, 0:7,:,:]
                val = numpy.multiply(emissionArr, lgridArr)
                # 每层都需要保留则注释此行，所有层全部累加则保留此行
                val = numpy.sum(val, axis = 1, keepdims=True)
                sens_species = sensfile.createVariable(var,'f4',('TSTEP','LAY','ROW','COL'))  
                sens_species[:] = numpy.zeros([ltime,llay,lrow,lcol])
                sens_species[:] = numpy.reshape(val,[ltime,llay,lrow,lcol])
                for varattr in varattrs:
                    varattrVal = getattr(lgrid.variables[var], varattr)
                    setattr(sens_species, varattr, varattrVal)
                setattr(sens_species, 'units', "CF")
            sensfile.close()
        emission.close()

    def makefiles_semi_normalized_sensitivity_daily(self):
        varLst = ['PMFINE', 'MEOH', 'POC', 'OLE', 'TERP', 'TOL', 'PEC', 'NO', 'ETHA', 'ETOH', 'ETH', 'NO2', 'NH3', 'FORM', 'PNO3', \
                'ALDX', 'PMC', 'IOLE', 'PAR', 'XYL', 'PSO4', 'ISOP', 'CO', 'ALD2', 'SO2']
        for day_index in range((self.date_end - self.date_start).days + 1):
            current_day  = self.date_start + datetime.timedelta(days=day_index)
            filename = f"{self.output_dir}/semi_normalized_sensitivity_em_Day_" + current_day.strftime('%Y-%m-%d') + ".nc"
            if os.path.exists(filename): os.system(f"rm {filename}")
            lgridPath = f"{self.output_dir}/semi_normalized_sensitivity_em_" + current_day.strftime('%Y-%m-%d') + ".nc"
            lgrid = netCDF4.Dataset(lgridPath, mode='r', open=True)
            tflagArr = numpy.array(lgrid.variables["TFLAG"][:][0:1, 0, :]).astype(int)
            conchr = lgrid.variables['PMFINE'][:][:][:][:]
            ltime, llay, lrow, lcol = conchr.shape
            ltime = 1
            ldattim = 2            
            sensfile = netCDF4.Dataset(filename, mode='w', format='NETCDF3_64BIT')
            attrs=["IOAPI_VERSION", "EXEC_ID", "FTYPE", "CDATE", "CTIME", "WDATE", "WTIME", "SDATE", "STIME", \
                            "TSTEP", "NTHIK", "NCOLS", "NROWS", "GDTYP", "P_ALP", "P_BET", "P_GAM", "XCENT", "YCENT", "XORIG",\
                            "YORIG", "XCELL", "YCELL", "VGTYP", "VGTOP", "VGLVLS", "GDNAM", "HISTORY"]
            for attr in attrs:
                if hasattr(lgrid, attr):
                    attrVal = getattr(lgrid, attr)
                    setattr(sensfile, attr, attrVal)
            setattr(sensfile,"NVARS", len(varLst)) # only one variable
            setattr(sensfile,"FILEDESC","sensitivity file: lgrid_em * emission")

            sensfile.createDimension("TSTEP", None)
            sensfile.createDimension("DATE-TIME", ldattim)
            sensfile.createDimension("LAY", llay)
            sensfile.createDimension("VAR", len(varLst))
            sensfile.createDimension("ROW", lrow)
            sensfile.createDimension("COL", lcol)
            sensfile.sync()
            dattim = numpy.zeros([ltime,len(varLst),ldattim])
            sens_tflag = sensfile.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
            for j in range(len(varLst)):
                dattim[:, j, :] = tflagArr
            sens_tflag[:] = dattim
            # sens_tflag[:] = numpy.reshape(dattim,[ltime,llay,ldattim])
            varattrs=["long_name","units","var_desc"]
            for varattr in varattrs:
                if hasattr(lgrid.variables['TFLAG'], varattr): 
                    varattrVal = getattr(lgrid.variables['TFLAG'], varattr)
                    setattr(sens_tflag, varattr, varattrVal)
            for var in varLst:
                lgridArr = lgrid.variables[var][:][:][:][:][:, :,:,:]
                val = numpy.sum(lgridArr, axis = 0, keepdims=True)
                sens_species = sensfile.createVariable(var,'f4',('TSTEP','LAY','ROW','COL'))  
                sens_species[:] = numpy.zeros([ltime,llay,lrow,lcol])
                sens_species[:] = numpy.reshape(val,[ltime,llay,lrow,lcol])
                for varattr in varattrs:
                    varattrVal = getattr(lgrid.variables[var], varattr)
                    setattr(sens_species, varattr, varattrVal)
                setattr(sens_species, 'units', "CF")
            sensfile.close()

    def makefiles_semi_normalized_sensitivity_period(self):
        filename = f"{self.output_dir}/semi_normalized_sensitivity_em_Tot_" + self.date_start.strftime('%Y-%m-%d') + "~" + self.date_end.strftime('%Y-%m-%d') + ".nc"
        if os.path.exists(filename): os.system(f"rm {filename}")
        sensfile = netCDF4.Dataset(filename, mode='w', format='NETCDF3_64BIT')    
        varLst = ['PMFINE', 'MEOH', 'POC', 'OLE', 'TERP', 'TOL', 'PEC', 'NO', 'ETHA', 'ETOH', 'ETH', 'NO2', 'NH3', 'FORM', 'PNO3', \
                'ALDX', 'PMC', 'IOLE', 'PAR', 'XYL', 'PSO4', 'ISOP', 'CO', 'ALD2', 'SO2']
        lgridPath = f"{self.output_dir}/semi_normalized_sensitivity_em_Day_" + self.date_start.strftime('%Y-%m-%d') + ".nc"
        lgrid = netCDF4.Dataset(lgridPath, mode='r', open=True)
        tflagArr = numpy.array(lgrid.variables["TFLAG"][:][0:1, 0, :]).astype(int)
        # shape : (25, 44, 106, 88)
        conchr = lgrid.variables['PMFINE'][:][:][:][:]
        ltime, llay, lrow, lcol = conchr.shape
        ldattim = 2
        attrs=["IOAPI_VERSION", "EXEC_ID", "FTYPE", "CDATE", "CTIME", "WDATE", "WTIME", "SDATE", "STIME", \
                        "TSTEP", "NTHIK", "NCOLS", "NROWS", "GDTYP", "P_ALP", "P_BET", "P_GAM", "XCENT", "YCENT", "XORIG",\
                        "YORIG", "XCELL", "YCELL", "VGTYP", "VGTOP", "VGLVLS", "GDNAM", "HISTORY"]
        for attr in attrs:
            if hasattr(lgrid, attr):
                attrVal = getattr(lgrid, attr)
                setattr(sensfile, attr, attrVal)
        setattr(sensfile,"NVARS", len(varLst)) # only one variable
        setattr(sensfile,"FILEDESC","sensitivity file: lgrid_em * emission")

        sensfile.createDimension("TSTEP", None)
        sensfile.createDimension("DATE-TIME", ldattim)
        sensfile.createDimension("LAY", llay)
        sensfile.createDimension("VAR", len(varLst))
        sensfile.createDimension("ROW", lrow)
        sensfile.createDimension("COL", lcol)
        sensfile.sync()

        sens_tflag = sensfile.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        dattim = numpy.zeros([ltime,len(varLst),ldattim])
        for j in range(len(varLst)):
            dattim[:, j, :] = tflagArr
        sens_tflag[:] = dattim
        varattrs=["long_name","units","var_desc"]
        for varattr in varattrs:
            if hasattr(lgrid.variables['TFLAG'], varattr): 
                varattrVal = getattr(lgrid.variables['TFLAG'], varattr)
                setattr(sens_tflag, varattr, varattrVal)
        ans = numpy.zeros([len(varLst), ltime, llay, lrow, lcol])
        for day_index in range((self.date_end - self.date_start).days + 1):
            current_day  = self.date_start + datetime.timedelta(days=day_index)
            lgridPath = f"{self.output_dir}/semi_normalized_sensitivity_em_Day_" + current_day.strftime('%Y-%m-%d') + ".nc"
            lgrid = netCDF4.Dataset(lgridPath, mode='r', open=True)
            for i, var in enumerate(varLst):
                lgridArr = lgrid.variables[var][:][:][:][:][:, :,:,:]
                ans[i, : :, :, :] = ans[i, : :, :, :] + lgridArr
        for i, var in enumerate(varLst):
            val = ans[i, :, :, :, :]
            sens_species = sensfile.createVariable(var,'f4',('TSTEP','LAY','ROW','COL'))  
            sens_species[:] = numpy.zeros([ltime,llay,lrow,lcol])
            sens_species[:] = numpy.reshape(val,[ltime,llay,lrow,lcol])
            for varattr in varattrs:
                varattrVal = getattr(lgrid.variables[var], varattr)
                setattr(sens_species, varattr, varattrVal)

    def makefiles_sensitivity_hourly(self, bwdDir):   
        varLst = ['PMFINE', 'MEOH', 'POC', 'OLE', 'TERP', 'TOL', 'PEC', 'NO', 'ETHA', 'ETOH', 'ETH', 'NO2', 'NH3', 'FORM', 'PNO3', \
                'ALDX', 'PMC', 'IOLE', 'PAR', 'XYL', 'PSO4', 'ISOP', 'CO', 'ALD2', 'SO2']
        for day_index in range((self.date_end - self.date_start).days + 1):
            current_day  = self.date_start + datetime.timedelta(days=day_index)
            filename = f"{self.output_dir}/sensitivity_em_" + current_day.strftime('%Y-%m-%d') + ".nc"
            if os.path.exists(filename): os.system(f"rm {filename}")
            lgridPath = f"{bwdDir}/lgrid_em_" + current_day.strftime('%Y-%m-%d') + ".nc"
            lgrid = netCDF4.Dataset(lgridPath, mode='r', open=True)
            # print(numpy.array(lgrid.variables["TFLAG"][:])[0:1, 0, :].shape)
            tflagArr = numpy.array(lgrid.variables["TFLAG"][:][1:, 0, :]).astype(int)
            # shape : (25, 44, 106, 88)
            conchr = lgrid.variables['PMFINE'][:][:][:][:]
            ltime, llay, lrow, lcol = conchr.shape
            ltime = 24
            llay = 1
            ldattim = 2
            sensfile = netCDF4.Dataset(filename, mode='w', format='NETCDF3_64BIT')
            attrs=["IOAPI_VERSION", "EXEC_ID", "FTYPE", "CDATE", "CTIME", "WDATE", "WTIME", "SDATE", "STIME", \
                            "TSTEP", "NTHIK", "NCOLS", "NROWS", "GDTYP", "P_ALP", "P_BET", "P_GAM", "XCENT", "YCENT", "XORIG",\
                            "YORIG", "XCELL", "YCELL", "VGTYP", "VGTOP", "VGLVLS", "GDNAM", "HISTORY"]
            for attr in attrs:
                if hasattr(lgrid, attr):
                    attrVal = getattr(lgrid, attr)
                    setattr(sensfile, attr, attrVal)
            setattr(sensfile,"NVARS", len(varLst)) # only one variable
            setattr(sensfile,"FILEDESC",getattr(lgrid, "FILEDESC"))

            sensfile.createDimension("TSTEP", None)
            sensfile.createDimension("DATE-TIME", ldattim)
            sensfile.createDimension("LAY", llay)
            sensfile.createDimension("VAR", len(varLst))
            sensfile.createDimension("ROW", lrow)
            sensfile.createDimension("COL", lcol)
            sensfile.sync()
            dattim = numpy.zeros([ltime,len(varLst),ldattim])
            sens_tflag = sensfile.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
            for j in range(len(varLst)):
                dattim[:, j, :] = tflagArr
            sens_tflag[:] = dattim
            varattrs=["long_name","units","var_desc"]
            for varattr in varattrs:
                if hasattr(lgrid.variables['TFLAG'], varattr): 
                    varattrVal = getattr(lgrid.variables['TFLAG'], varattr)
                    setattr(sens_tflag, varattr, varattrVal)
            for var in varLst:
                lgridArr = lgrid.variables[var][:][:][:][:][1:, :,:,:]
                val = numpy.sum(lgridArr, axis = 1, keepdims=True)
                sens_species = sensfile.createVariable(var,'f4',('TSTEP','LAY','ROW','COL'))  
                sens_species[:] = numpy.zeros([ltime,llay,lrow,lcol])
                sens_species[:] = numpy.reshape(val,[ltime,llay,lrow,lcol])
                for varattr in varattrs:
                    varattrVal = getattr(lgrid.variables[var], varattr)
                    setattr(sens_species, varattr, varattrVal)
                # setattr(sens_species, 'units', "CF")
            sensfile.close()

    def makefiles_sensitivity_daily(self): 
        varLst = ['PMFINE', 'MEOH', 'POC', 'OLE', 'TERP', 'TOL', 'PEC', 'NO', 'ETHA', 'ETOH', 'ETH', 'NO2', 'NH3', 'FORM', 'PNO3', \
                'ALDX', 'PMC', 'IOLE', 'PAR', 'XYL', 'PSO4', 'ISOP', 'CO', 'ALD2', 'SO2']
        for day_index in range((self.date_end - self.date_start).days + 1):
            current_day  = self.date_start + datetime.timedelta(days=day_index)
            filename = f"{self.output_dir}/sensitivity_em_Day_" + current_day.strftime('%Y-%m-%d') + ".nc"
            if os.path.exists(filename): os.system(f"rm {filename}")
            lgridPath = f"{self.output_dir}/sensitivity_em_" + current_day.strftime('%Y-%m-%d') + ".nc"
            lgrid = netCDF4.Dataset(lgridPath, mode='r', open=True)
            tflagArr = numpy.array(lgrid.variables["TFLAG"][:][0:1, 0, :]).astype(int)
            conchr = lgrid.variables['PMFINE'][:][:][:][:]
            ltime, llay, lrow, lcol = conchr.shape
            ltime = 1
            ldattim = 2
            sensfile = netCDF4.Dataset(filename, mode='w', format='NETCDF3_64BIT')
            attrs=["IOAPI_VERSION", "EXEC_ID", "FTYPE", "CDATE", "CTIME", "WDATE", "WTIME", "SDATE", "STIME", \
                            "TSTEP", "NTHIK", "NCOLS", "NROWS", "GDTYP", "P_ALP", "P_BET", "P_GAM", "XCENT", "YCENT", "XORIG",\
                            "YORIG", "XCELL", "YCELL", "VGTYP", "VGTOP", "VGLVLS", "GDNAM", "HISTORY"]
            for attr in attrs:
                if hasattr(lgrid, attr):
                    attrVal = getattr(lgrid, attr)
                    setattr(sensfile, attr, attrVal)
            setattr(sensfile,"NVARS", len(varLst)) # only one variable
            setattr(sensfile,"FILEDESC","lgrid_em file: sum of one day")

            sensfile.createDimension("TSTEP", None)
            sensfile.createDimension("DATE-TIME", ldattim)
            sensfile.createDimension("LAY", llay)
            sensfile.createDimension("VAR", len(varLst))
            sensfile.createDimension("ROW", lrow)
            sensfile.createDimension("COL", lcol)
            sensfile.sync()
            dattim = numpy.zeros([ltime,len(varLst),ldattim])
            sens_tflag = sensfile.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
            for j in range(len(varLst)):
                dattim[:, j, :] = tflagArr
            sens_tflag[:] = dattim
            varattrs=["long_name","units","var_desc"]
            for varattr in varattrs:
                if hasattr(lgrid.variables['TFLAG'], varattr): 
                    varattrVal = getattr(lgrid.variables['TFLAG'], varattr)
                    setattr(sens_tflag, varattr, varattrVal)
            for var in varLst:
                lgridArr = lgrid.variables[var][:][:][:][:][:, :,:,:]
                val = numpy.sum(lgridArr, axis = 0, keepdims=True)
                sens_species = sensfile.createVariable(var,'f4',('TSTEP','LAY','ROW','COL'))  
                sens_species[:] = numpy.zeros([ltime,llay,lrow,lcol])
                sens_species[:] = numpy.reshape(val,[ltime,llay,lrow,lcol])
                for varattr in varattrs:
                    varattrVal = getattr(lgrid.variables[var], varattr)
                    setattr(sens_species, varattr, varattrVal)
                # setattr(sens_species, 'units', "CF")
            sensfile.close()

    def makefiles_sensitivity_period(self):
        filename = f"{self.output_dir}/sensitivity_em_Tot_" + self.date_start.strftime('%Y-%m-%d') + "~" + self.date_end.strftime('%Y-%m-%d') + ".nc"
        if os.path.exists(filename): os.system(f"rm {filename}")
        sensfile = netCDF4.Dataset(filename, mode='w', format='NETCDF3_64BIT')  
        varLst = ['PMFINE', 'MEOH', 'POC', 'OLE', 'TERP', 'TOL', 'PEC', 'NO', 'ETHA', 'ETOH', 'ETH', 'NO2', 'NH3', 'FORM', 'PNO3', \
                'ALDX', 'PMC', 'IOLE', 'PAR', 'XYL', 'PSO4', 'ISOP', 'CO', 'ALD2', 'SO2']
        lgridPath = f"{self.output_dir}/sensitivity_em_Day_" + self.date_start.strftime('%Y-%m-%d') + ".nc"
        lgrid = netCDF4.Dataset(lgridPath, mode='r', open=True)
        tflagArr = numpy.array(lgrid.variables["TFLAG"][:][0:1, 0, :]).astype(int)
        conchr = lgrid.variables['PMFINE'][:][:][:][:]
        ltime, llay, lrow, lcol = conchr.shape
        ldattim = 2
        attrs=["IOAPI_VERSION", "EXEC_ID", "FTYPE", "CDATE", "CTIME", "WDATE", "WTIME", "SDATE", "STIME", \
                        "TSTEP", "NTHIK", "NCOLS", "NROWS", "GDTYP", "P_ALP", "P_BET", "P_GAM", "XCENT", "YCENT", "XORIG",\
                        "YORIG", "XCELL", "YCELL", "VGTYP", "VGTOP", "VGLVLS", "GDNAM", "HISTORY"]
        for attr in attrs:
            if hasattr(lgrid, attr):
                attrVal = getattr(lgrid, attr)
                setattr(sensfile, attr, attrVal)
        setattr(sensfile,"NVARS", len(varLst)) # only one variable
        setattr(sensfile,"FILEDESC","total sensitivity file")

        sensfile.createDimension("TSTEP", None)
        sensfile.createDimension("DATE-TIME", ldattim)
        sensfile.createDimension("LAY", llay)
        sensfile.createDimension("VAR", len(varLst))
        sensfile.createDimension("ROW", lrow)
        sensfile.createDimension("COL", lcol)
        sensfile.sync()
        sens_tflag = sensfile.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        dattim = numpy.zeros([ltime,len(varLst),ldattim])
        for j in range(len(varLst)):
            dattim[:, j, :] = tflagArr
        sens_tflag[:] = dattim
        varattrs=["long_name","units","var_desc"]
        for varattr in varattrs:
            if hasattr(lgrid.variables['TFLAG'], varattr): 
                varattrVal = getattr(lgrid.variables['TFLAG'], varattr)
                setattr(sens_tflag, varattr, varattrVal)
        ans = numpy.zeros([len(varLst), ltime, llay, lrow, lcol])
        for day_index in range((self.date_end - self.date_start).days + 1):
            current_day  = self.date_start + datetime.timedelta(days=day_index)
            lgridPath = f"{self.output_dir}/sensitivity_em_Day_" + current_day.strftime('%Y-%m-%d') + ".nc"
            lgrid = netCDF4.Dataset(lgridPath, mode='r', open=True)
            for i, var in enumerate(varLst):
                lgridArr = lgrid.variables[var][:][:][:][:][:, :,:,:]
                ans[i, : :, :, :] = ans[i, : :, :, :] + lgridArr
        for i, var in enumerate(varLst):
            val = ans[i, :, :, :, :]
            sens_species = sensfile.createVariable(var,'f4',('TSTEP','LAY','ROW','COL'))  
            sens_species[:] = numpy.zeros([ltime,llay,lrow,lcol])
            sens_species[:] = numpy.reshape(val,[ltime,llay,lrow,lcol])
            for varattr in varattrs:
                varattrVal = getattr(lgrid.variables[var], varattr)
                setattr(sens_species, varattr, varattrVal)

    def draw_Tropomi_PNG(self, dataset, datetime, factor_index):
        factor_name = self.factors[factor_index].split('|')[0]
        print('factor name in mp4_png', factor_name)
        print('type of dataset:', type(dataset))
        print('type of dataset[2]:', type(dataset[2]))        
        # 绘制图，设置投影
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(1, 1, 1, projection=self.domain.LCCProj_crs)
        # 设置显示范围
        ax.set_extent(self.domain.LCC_plotextent, crs=self.domain.LCCProj_crs)   
        # ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())     
        # 读取并处理地图数据
        dem_file = self.root_dir + "/data/ETOPO2v2g_f4.nc"
        f = xarray.open_dataset(dem_file)
        lon = f['x'].sel(x=slice(50, 155)).values
        lat = f['y'].sel(y=slice(0, 70)).values
        dem = f['z'].sel(x=slice(50, 155), y=slice(0, 70)).values
        lon_a, lat_a = numpy.meshgrid(lon, lat)
        ret = ax.projection.transform_points(ccrs.PlateCarree(), numpy.array(lon_a), numpy.array(lat_a))
        
        xx = ret[..., 0]
        yy = ret[..., 1]
        dem_levels = [-8000, -6000, -4000, -2000, -1000, -200, -50, 0, 50, 200, 500, 1000, 1500, 2000, 3000, 4000, 5000,
                      6000, 7000, 8000]
        dem_color = ['#084594', '#2171b5', '#4292c6', '#6baed6', '#9ecae1', '#c6dbef', '#deebf7', '#006837', '#31a354',
                     '#78c679', '#addd8e', \
                     '#d9f0a3', '#f7fcb9', '#c9bc87', '#a69165', '#856b49', '#664830', '#ad9591', '#d7ccca']
        # 绘制地图
        # print('xx:', xx.shape)
        # print('yy:', yy.shape)
        # print('dem:', dem.shape)

        ax.contourf(xx, yy, dem, levels=dem_levels, colors=dem_color, extend='both')
        # plot city's mark
        df = pandas.read_csv(self.root_dir + "/utils/asia_city.csv")
        df_admin = df[df['capital'] == 'admin']
        lon_admin = df_admin['lng'].to_list()
        lat_admin = df_admin['lat'].to_list()
        name_admin = df_admin['city'].to_list()
        zip_object = zip(lon_admin, lat_admin, name_admin)
        # zorder 控制绘制顺序，值越小越先绘制
        for (lg, lt, ne) in zip_object:
            ax.text(lg, lt, ne, color='dimgrey', va='center', ha='center', transform=ccrs.Geodetic(), zorder=4,
                    fontsize=8)
        # plot capital's mark
        df_primary = df[df['capital'] == 'primary']
        lon_primary = df_primary['lng'].to_list()
        lat_primary = df_primary['lat'].to_list()
        name_primary = df_primary['city'].to_list()
        zip_object = zip(lon_primary, lat_primary, name_primary)
        for (lg, lt, ne) in zip_object:
            ax.plot(lg, lt, "ro", zorder=5, markersize=1, transform=ccrs.Geodetic())
            ax.text(lg + 1, lt - 1, ne, va='center', color='k', ha='center', transform=ccrs.Geodetic(), zorder=4,
                    fontsize=8)
        lon_ad, lat_ad = numpy.meshgrid(dataset[0], dataset[1])
        ret = ax.projection.transform_points(ccrs.PlateCarree(), numpy.array(lon_ad), numpy.array(lat_ad))
        x = ret[..., 0]
        y = ret[..., 1]
        # 获取数据范围
        # 创建掩码，将NaN值标记为True
        mask = numpy.isnan(dataset[2])

        # 使用掩码将NaN值替换为灰色
        dataset_mask = numpy.where(mask, -999.0, dataset[2])
        minvalue = dataset_mask.min()
        maxvalue = dataset_mask.max()
        
        #minvalue = float(self.factors[factor_index].split('|')[2].split('~')[0])
        #maxvalue = float(self.factors[factor_index].split('|')[2].split('~')[1])
        # 自定义颜色条
        
        # conc_color = ['blue', 'cyan', 'yellow', 'orange', 'red']
        # mycmap = col.LinearSegmentedColormap.from_list('own', conc_color, N=len(conc_color))
        # my_cmap_aplpha = mycmap(numpy.arange(mycmap.N))
        # for i in range(my_cmap_aplpha.shape[0]):
        #     if i < 1:
        #         my_cmap_aplpha[i, -1] = 0.9
        #     else:
        #         my_cmap_aplpha[i, -1] = 0.9
        # my_cmap_aplpha = col.ListedColormap(my_cmap_aplpha, name='myowncolor')
        # 在地图上绘制数据
        conc_color = ['white', 'blue', 'cyan', 'yellow', 'red']
        # c_levels = [x for x in numpy.linspace( minvalue, maxvalue, num=len(conc_color) + 1)]
        # c_levels = [-999] + c_levels
        c_levels = [-999, 0, 6, 12, 18, 24, 30]
        my_norm = col.BoundaryNorm(c_levels, 6)
        my_cmap = col.ListedColormap(['grey', 'white', 'blue', 'cyan', 'yellow', 'red'])
        p = ax.pcolormesh(x, y, dataset_mask, norm=my_norm,\
            cmap=my_cmap)
        #p = ax.contourf(x, y, dataset_mask, norm=my_norm,\
        #    cmap=my_cmap)
        # 添加颜色条

        # 使用 plt.axes 添加独立的坐标轴，并创建颜色条
        cax = plt.axes([0.9, 0.15, 0.03, 0.7])  # [left, bottom, width, height]
        cbar = plt.colorbar(p, cax=cax, ticks=c_levels)

        # 在颜色条上添加NaN字符标签
        cbar.ax.set_yticklabels([f'{val}' if not val==-999 else 'NaN' for val in c_levels])

        plt.text(0.2, 0.95, datetime + ' (BeiJing Time)', fontsize=12, \
                 transform=ax.transAxes, color='r', bbox=dict(facecolor='w', pad=1, alpha=0.7, edgecolor='none'))
        # 绘制标题
        title_str = factor_name + ' (' + self.factors[factor_index].split('|')[3] + ',' + \
            self.factors[factor_index].split('|')[1] + ')'
        ax.set_title(title_str, fontdict={'size': 16, 'weight': 'bold'})
        # 设置坐标标签、网格线
        ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
        ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
        self.lon_tickrange = numpy.arange(80,140,10)  # the range of lon to be ticked
        self.lat_tickrange = numpy.arange(10,60,10) # the range of lat to be ticked
        gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5, alpha=1,
                          x_inline=False, y_inline=False, color='gray')
        gl.top_labels = False
        gl.right_labels = False
        gl.xpadding = 8
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlocator = matplotlib.ticker.FixedLocator(self.lon_tickrange)
        gl.ylocator = matplotlib.ticker.FixedLocator(self.lat_tickrange)
        gl.xlabel_style = {'rotation': 0, 'ha': 'center'}
        gl.ylabel_style = {'rotation': 0}
        ax.spines['right'].set_visible(True)
        ax.spines['left'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        ax.spines['top'].set_visible(True)
        # 保存图片
        savename = self.output_dir + '/' + 'Tropomi_' + factor_name + '_' + datetime + '.png'
        print('savename in mp4:',savename) 
        plt.savefig(savename, bbox_inches='tight', dpi=600, format='PNG')
        plt.close()

    def Plot_Tropomi_PNGs(self,factor_index): # The main function of plot 
        #num_process = ((self.end_day - self.start_day).days + 1) * 24
        num_process = (self.end_day - self.start_day).days + 1
        num_nodes = math.ceil(num_process/self.mpijob_npn)
        for i_node in range(num_nodes):
            job_tag = self.output_dir + '/plot_TROPOMI_' + str(i_node+1)
            python_file = job_tag + '.py'
            python_job = open(python_file, 'w')
            python_job.write("import os\n")
            python_job.write("import datetime\n")
            python_job.write("import cartopy.crs as ccrs\n")
            python_job.write("import PIL\n")
            python_job.write("import pandas\n")
            python_job.write("import xarray\n")
            python_job.write("import operator\n")
            python_job.write("import math\n")
            python_job.write("from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER\n")
            python_job.write("import stat\n")
            python_job.write("import matplotlib.animation\n")
            python_job.write("import matplotlib.colors as col\n")
            python_job.write("import matplotlib.pyplot as plt\n")
            python_job.write("import sys\n")
            python_job.write("import numpy\n")
            python_job.write('sys.path.append("'+ self.root_dir + '")\n')
            python_job.write("import post_process\n")
            python_job.write("import cv2\n")
            python_job.write("import multiprocessing\n")
            python_job.write("from osgeo import gdal\n")
            python_job.write("from osgeo import osr\n")
            # python_job.write("obj_TROPOMI = post_process.Post_Process()\n")
            python_job.write("obj_TROPOMI = post_process.Post_Process('" + self.data_type \
                + "','" + self.date_start.strftime('%Y-%m-%d') + "','" + self.date_end.strftime('%Y-%m-%d') \
                 + "','" + str(self.root_dir) + "','" + str(self.source_dirs) + "')\n")
            python_job.write("obj_TROPOMI.output_dir = '" + str(self.output_dir) + "'\n")
            python_job.write("print(obj_TROPOMI, obj_TROPOMI.factors, ' begin')\n")
            python_job.write("def call_draw_PNG(day_index,factor_index):\n")
            python_job.write("    print('===============*****')\n")  
            python_job.write("    factor_name = obj_TROPOMI.factors[factor_index].split('|')[0]\n")
            python_job.write("    factor_version =  obj_TROPOMI.factors[factor_index].split('|')[4]\n")
            python_job.write("    current_time = obj_TROPOMI.start_day + datetime.timedelta(days=day_index)\n")
            python_job.write("    dt_str = current_time.strftime('%Y-%m-%d')\n")                     
            python_job.write("    filename = 'Tropomi_' + 'col_' + factor_name + '_' + current_time.strftime('%Y%m%d') + '_' + factor_version +'_gridded.tif'\n")
            #python_job.write("    df = gdal.Open(filename)\n")             
            # python_job.write("    if os.path.isfile(obj_TROPOMI.source_path + '/' + factor_name + '/' + filename) :\n")
            # python_job.write("       df = gdal.Open(obj_TROPOMI.source_path + '/' + factor_name + '/' + filename,)\n")
            # python_job.write("    else:\n")
            # python_job.write("       df = gdal.Open(obj_TROPOMI.source_path + '/null_values.tif')\n")
            python_job.write("    df = gdal.Open(obj_TROPOMI.source_path + '/' + factor_name + '/' + filename,)\n")
            python_job.write("    geo_transform = df.GetGeoTransform()\n")  # 获取tif数据的基本信息
            python_job.write("    prosrs = osr.SpatialReference()\n")#获得给定数据的投影参考系和地理参考系
            python_job.write("    prosrs.ImportFromWkt(df.GetProjection())\n")#投影参考系
            python_job.write("    geosrs = prosrs.CloneGeogCS()\n")#地理参考系    
            python_job.write("    origin_x = geo_transform[0]\n")  # 左上角经度
            python_job.write("    origin_y = geo_transform[3]\n")  # 左上角纬度
            python_job.write("    pixel_width = geo_transform[1]\n")  #像元宽度          
            python_job.write("    pixel_height = geo_transform[5]\n") #像元高度             
            python_job.write("    n_rows = df.RasterYSize\n")  #行数          
            python_job.write("    n_cols = df.RasterXSize\n")  #列数            
            python_job.write("    in_band = df.GetRasterBand(1)\n")  #打开波段            
            python_job.write("    print(in_band)\n")  
            python_job.write("    ds = in_band.ReadAsArray()\n")   #读数据
            python_job.write("    print(ds)\n") 
            python_job.write("    lon = numpy.linspace(origin_x,(origin_x+pixel_width*n_cols-pixel_width),n_cols)\n")  
            python_job.write("    lat = numpy.linspace(origin_y,(origin_y + pixel_height*n_rows-pixel_height),n_rows)\n")               
            python_job.write("    scale_factor = 10000\n")
            python_job.write("    ds = ds * scale_factor\n")
            python_job.write("    print('***********')\n")      
            #python_job.write("    ds = numpy.mean(ds)\n")
            python_job.write("    print(obj_TROPOMI, obj_TROPOMI.factors)\n")        
            python_job.write("    for factor_index in range(len(obj_TROPOMI.factors)):\n")
            python_job.write("        factor_name = obj_TROPOMI.factors[factor_index].split('|')[0]\n")
            # python_job.write("        print(factor_name)\n")
            python_job.write("    obj_TROPOMI.draw_Tropomi_PNG((lon, lat, ds), dt_str, factor_index)\n")
            #python_job.write("    obj_TROPOMI.draw_upload_PNG((lon, lat, ds), dt_str, factor_index)\n")
            if i_node == (num_nodes-1):
                range_process = num_process - (i_node * self.mpijob_npn)
            else:
                range_process = self.mpijob_npn
            for process_node in range(range_process):
                index_process = i_node * self.mpijob_npn + process_node
                print(index_process,factor_index)
                python_job.write("multiprocessing.Process(target=call_draw_PNG, args=("+str(index_process)+", "+str(factor_index)+")).start()\n")
            python_job.close()
            os.chmod(python_file, stat.S_IRWXO+stat.S_IRWXG+stat.S_IRWXU)
            slurm_fname = job_tag + '.slurm'
            slurm_file = open(slurm_fname, 'w')
            slurm_file.write("#!/bin/bash\n")
            slurm_file.write("#SBATCH --job-name=plot_TROPOMI\n")
            slurm_file.write("#SBATCH --time=01:00:00\n")     # maximum running time of one hour
            slurm_file.write("#SBATCH --partition=" + self.serialjob_partition + "\n")
            slurm_file.write("#SBATCH --mail-type=end\n")
            slurm_file.write("#SBATCH --mail-user=\n")
            slurm_file.write("#SBATCH --output=" + self.output_dir + "/%j.out\n")
            slurm_file.write("#SBATCH --error=" + self.output_dir + "/%j.err\n")
            slurm_file.write("#SBATCH -N 1\n")
            slurm_file.write("#SBATCH --ntasks-per-node=" + str(self.mpijob_npn) + "\n")
            slurm_file.write('python' + ' ' + python_file + ' >& ' + python_file + '.log')
            slurm_file.close()
            os.system('sbatch ' + slurm_fname)          
            
    def Timeseries(self, factor_index, site_index, df_obs, df_simu_sites):  # 各站点各物种的时序对比图        
        factor_name = self.factors[factor_index].split('|')[0]
        factor_unit = self.factors[factor_index].split('|')[1]
        site_name = self.sitenames[site_index]
        site_province = self.siteprovince[site_index]
        site_city = self.sitecities[site_index]
        site_id = self.siteIDs[site_index]
        site_lon = self.sitelons[site_index]
        site_lat = self.sitelats[site_index]
        if site_province is None:
            str_province = 'null'
        else:
            str_province = str(site_province)
        if site_city is None:
            str_city = 'null'
        else:
            str_city = str(site_city)
        if site_name is None:
            str_site = 'null'
        else:
            str_site = str(site_name)
        df_obs = df_obs[df_obs['SiteCode']==site_id]
        x = df_obs.copy()
        x.loc[:,'TimePoint'] = pandas.to_datetime(df_obs.loc[:,'TimePoint'], format="%Y-%m-%d %H:%M:%S") 
        df_obs = x
        if factor_name in df_obs.columns:
            df_obs = df_obs.dropna(subset=[factor_name])
            df_obs = df_obs.drop(df_obs[(df_obs[factor_name]<-10000)|(df_obs[factor_name]>200000)].index) # 删除异常值的行
            if df_obs[[factor_name]].count().iloc[0] == 0 and len(self.source_dirs.split(';'))==1:
                return   # 如果观测站没有数据，且只有一个模拟情景，无需画图提前结束
            # fill the missing obs datatime with NAN values     
            if len(df_obs) != len(self.date_cover):
                dt_cover = pandas.DataFrame(self.date_cover)
                dt_cover.columns = ['TimePoint']
                df_obs = pandas.merge(dt_cover, df_obs, how='left', on='TimePoint')
            df_obs = df_obs[[factor_name]]
            df_obs.columns = [factor_name + "_obs"]
            df_obs.reset_index(drop=True,inplace=True)
        else: 
            if len(self.source_dirs.split(';'))==1: 
                return  # 如果观测站没有数据，且只有一个模拟情景，无需画图提前结束
        fig, ax = plt.subplots(figsize=(16,10))
        pandas.plotting.register_matplotlib_converters() 
        ax.set_title(factor_name + ',' + str_site +'(' + str_province + ','+ str_city + ',' + str(site_lon) \
            + 'E, ' + str(site_lat) + 'N)', fontsize=23, pad=(len(self.source_dirs.split(';'))+1)*35,fontproperties=self.heiti_font)        
        color_list = ['red','blue','green','black']
        linestyle = ['-',':','--','-.']        
        df_simu_sites = df_simu_sites[(df_simu_sites['factor_index'].astype(str) == str(factor_index)) & (df_simu_sites['SiteCode'].astype(str) == str(self.siteIDs[site_index]))]
        df_simu_sites = df_simu_sites.rename(columns={'factor_value':factor_name})
        df_simu_sites.sort_values(by=['TimePoint','source_index'],ascending=[True,True])
        df_simu_sites = df_simu_sites.drop(['SiteCode'], axis=1)
        df_simu_sites.reset_index(drop=True,inplace=True)            
        factor_min = df_simu_sites[factor_name].min() * 0.9
        factor_max = df_simu_sites[factor_name].max() * 1.1
        for source_index in range(len(self.source_dirs.split(';'))):
            df_simu_source = df_simu_sites[(df_simu_sites['source_index'] == source_index)]
            df_simu_source = df_simu_source.sort_values(by=['TimePoint'],ascending=[True])
            df_simu_source = df_simu_source.drop(['source_index'], axis=1)
            df_simu_source = df_simu_source.drop(['TimePoint'], axis=1)
            df_simu_source.reset_index(drop=True,inplace=True)      
            ax.plot(self.date_cover, df_simu_source[factor_name], linestyle[source_index], color=color_list[source_index],\
                label=self.source_dirs.split(';')[source_index].split('/')[-1], linewidth=3)
        if df_obs[[factor_name + "_obs"]].count().iloc[0] > 0:
            factor_min = min(df_obs[factor_name + "_obs"].min()*0.9, factor_min)
            factor_max = max(df_obs[factor_name + "_obs"].max()*1.1, factor_max)
        ax.plot_date(self.date_cover, df_obs, label='Observation', color=color_list[-1])
        ax.legend(fontsize=20, borderpad=0, ncol=1,frameon=0, facecolor='white', framealpha=0, \
            bbox_to_anchor=(0., 1, 1., (len(self.source_dirs.split(';'))+1)*0.06), loc='upper center')
        ax.tick_params(which='minor', length=3, width=1.5)
        ax.tick_params(labelsize=20, length=6, width=3)
        plt.xlabel("Date time", fontsize=20,labelpad=5,fontproperties=self.timesnr_font)
        plt.ylabel(factor_name + ' (' + factor_unit + ')',fontproperties=self.timesnr_font,fontsize=20,labelpad=5)
        ax.set_xlim(self.date_cover[0]-12*self.date_cover.freq, self.date_cover[len(self.date_cover)-1]+12*self.date_cover.freq)
        ax.set_ylim([factor_min, factor_max])
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['top'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.yaxis.set_major_locator(matplotlib.ticker.AutoLocator())
        ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        ax.xaxis.set_major_locator(matplotlib.dates.DayLocator(interval=math.ceil(len(self.date_cover)/7/24)))
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%d\n%Y'))
        ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        figname = factor_name.replace('_','-') + '_' +  str_province + '_' + str_city + '_' + str_site + '_' + str(site_lon) + '_' + str(site_lat) + '.png'
        fig.savefig(self.output_dir  + '/Timeseries/' + figname, bbox_inches='tight')        
        matplotlib.pyplot.close('all')        

    def Scatter(self, factor_index, source_index, df_obs, df_simu_sites): # 某情景某物种所有观测数据和模拟数据的散点对比图
        factor_name = self.factors[factor_index].split('|')[0]
        factor_unit = self.factors[factor_index].split('|')[1]
        factor_desc = self.factors[factor_index].split('|')[3]
        if factor_name not in df_obs.columns or df_obs[factor_name].count() < 1: # 如果没有观测数据吗，就跳出该物种
            return
        df_obs = df_obs.dropna(subset=[factor_name])
        df_obs = df_obs.drop(df_obs[(df_obs[factor_name]<-10000)|(df_obs[factor_name]>200000)].index) # 删除异常值的行
        df_obs = df_obs.rename(columns={factor_name:factor_name+'_obs'})
        df_obs = df_obs.sort_values(by=['TimePoint','SiteCode'],ascending=[True,True])
        x = df_obs.copy()
        x.loc[:,'TimePoint'] = pandas.to_datetime(df_obs.loc[:,'TimePoint'], format="%Y-%m-%d %H:%M:%S") 
        df_obs = x
        df_obs['SiteCode'] = df_obs['SiteCode'].astype(str)
        df_simu = df_simu_sites[(df_simu_sites['factor_index'].astype(str) == str(factor_index)) & (df_simu_sites['source_index'].astype(str) == str(source_index))]
        df_simu = df_simu.reset_index(drop=True)
        df_simu = df_simu.sort_values(by=['TimePoint','SiteCode'],ascending=[True,True])
        df_simu = df_simu.rename(columns={'factor_value':factor_name})
        df_simu = pandas.merge(df_simu, df_obs, on=['TimePoint','SiteCode'])  # 保证两个数据集大小一致
        df_simu.to_csv(self.output_dir + '/CSV/' + factor_name + '_source' + str(source_index) + '.csv', sep=',', index=False, header=True)
        df_obs = df_simu[[factor_name+'_obs']]
        df_simu = df_simu[[factor_name]]
        x = df_obs[factor_name + '_obs'].values
        y = df_simu[factor_name].values
        x_max = numpy.percentile(x, 99)
        y_max = numpy.percentile(y, 99)
        if x_max > y_max:
            max_value = math.ceil(x_max)
        else:
            max_value = math.ceil(y_max)
        x_min = numpy.percentile(x, 1)
        y_min = numpy.percentile(y, 1)
        if x_min < y_min:
            min_value = math.ceil(x_min)
        else:
            min_value = math.ceil(y_min)
        lr = LinearRegression(fit_intercept=True)
        lr.fit(x.reshape(-1, 1),y.reshape(-1, 1))
        K = lr.coef_
        b = lr.intercept_
        lr0 = LinearRegression(fit_intercept=False)
        lr0.fit(df_obs[[factor_name + '_obs']],df_simu[factor_name])
        K0 = lr0.coef_
        R2 = metrics.r2_score(y,lr.predict(x.reshape(-1, 1)))
        MAE = metrics.mean_absolute_error(y, lr.predict(x.reshape(-1, 1)))
        MSE = numpy.sqrt(metrics.mean_squared_error(y, lr.predict(x.reshape(-1, 1))))
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        # ===========Calculate the point density==========
        xy = numpy.vstack([x, y])
        if xy.size > 10000:
            torun = numpy.array_split(xy, 64, axis=1)
            pool = multiprocessing.Pool(processes = 64)
            results = pool.map(calc_kernel, torun)
            z = numpy.concatenate(results)
        else: z = calc_kernel(xy)
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
        dense=ax.scatter(x,y,8,marker='o',c=z*100,linewidths=0,cmap='jet') 
        clb=plt.colorbar(plt.cm.ScalarMappable(cmap='jet', norm=matplotlib.colors.Normalize(vmin=0, vmax=100)),ax=ax)
        plt.plot([min_value, max_value],[min_value, max_value],color='black',linestyle='--',lw=2)
        plt.plot(numpy.array([min_value, max_value]).reshape(-1, 1), lr.predict(numpy.array([min_value, max_value]).reshape(-1, 1)),'red',lw=2)
        plt.axis([0, max_value,0, max_value])
        plt.xlabel('Observation ('+factor_unit+')', weight="bold", fontsize=10, labelpad=5,fontproperties=self.timesnr_font)
        plt.ylabel(self.source_dirs.split(';')[source_index].split('/')[-1] + ' (' + factor_unit + ')', weight="bold", fontsize=10, labelpad=5,fontproperties=self.timesnr_font)
        plt.xticks(fontproperties=self.timesnr_font)
        plt.yticks(fontproperties=self.timesnr_font)
        plt.text(0.02, 0.96, '$N=%.f$' % len(y), fontsize=7, fontproperties=self.timesnr_font,transform=ax.transAxes,\
            bbox=dict(facecolor='w', alpha=1, edgecolor='none'))
        plt.text(0.02, 0.92, '$R^2=%.2f$' % R2, fontsize=7, fontproperties=self.timesnr_font,transform=ax.transAxes,\
            bbox=dict(facecolor='w', alpha=1, edgecolor='none'))
        if b>0:
            ope_str = '+'
        else:
            ope_str = ''
        plt.text(0.02, 0.88, '$Sim.=%.2f*Obs.$' % K + ope_str + '$%.2f$' % b, fontsize=7, fontproperties=self.timesnr_font,\
            bbox=dict(facecolor='w', alpha=1, edgecolor='none'), transform=ax.transAxes)
        plt.text(0.02, 0.84, '$MAE=%.2f$' % MAE + ' (' + factor_unit + ')', fontsize=7, fontproperties=self.timesnr_font,\
            bbox=dict(facecolor='w', alpha=1, edgecolor='none'), transform=ax.transAxes)
        plt.text(0.02, 0.80, '$MSE=%.2f$' % MSE + ' (' + factor_unit + ')', fontsize=7, fontproperties=self.timesnr_font,\
            bbox=dict(facecolor='w', alpha=1, edgecolor='none'), transform=ax.transAxes)
        plt.text(0.02, 0.76, '$Sim.=%.2f$*Obs.' % K0, fontsize=7, fontproperties=self.timesnr_font,\
            bbox=dict(facecolor='w', alpha=1, edgecolor='none'), transform=ax.transAxes)
        plt.xlim(min_value, max_value)
        plt.ylim(min_value, max_value)
        ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        ax.set_title(factor_desc + '(' + factor_name + '), (' + self.date_start.strftime('%Y-%m-%d') +' ~ '\
            + self.date_end.strftime('%Y-%m-%d') + ')', fontsize=12, fontproperties=self.timesnr_font, pad=5)
        figname = 'Scatter_' + factor_name.replace('_','-') + '_' + self.data_type + '_' + self.source_dirs.split(';')[source_index].split('/')[-1] + '.png'
        plt.savefig(self.output_dir + '/Scatter/' + figname,dpi=800,bbox_inches='tight')       
    
    def extract_2dsimu(self, source_index, layer_index, factor_index, day_index): # Extract the data record from corresponding nc file 
        current_day  = self.date_start + datetime.timedelta(days=day_index)
        factor_name = self.factors[factor_index].split('|')[0]
        layer = self.extracted_layers[layer_index]
        file_path = self.source_dirs.split(';')[int(source_index)]                   
        if self.data_type == 'WRF': # for WRF
            hour_start = 12 # start from 12 for our wrf running
            if ('-' not in str(source_index)):  #绝对值
                nc_file = file_path + '/wrfout_d01_' \
                    + (current_day + datetime.timedelta(days = -1)).strftime('%Y-%m-%d_12:00:00')
                nc_ds = netCDF4.Dataset(nc_file)
                Domain_NCOLS = int(nc_ds.getncattr('WEST-EAST_PATCH_END_UNSTAG')) # len(nc_ds.dimensions['west_east']) also works
                Domain_NROWS = int(nc_ds.getncattr('SOUTH-NORTH_PATCH_END_UNSTAG')) # len(nc_ds.dimensions['south_north']) also works
                ds_target = numpy.zeros([24, Domain_NROWS, Domain_NCOLS])
                ds = xarray.open_dataset(nc_file).isel({'Time': slice(13, None)}).sel(bottom_top=int(layer)-1)
                if factor_name == 'TC':
                    ds_target[:] = (ds['T']+300) * ((ds['P'] + ds['PB'])/100000)**(2/7) - 273.15
                if factor_name == 'RH':
                    TK = (ds['T']+300) * ((ds['P'] + ds['PB'])/100000)**(2/7)
                    ds_target[:] = 100 * ((ds['P'] + ds['PB']) * ds['QVAPOR']/(ds['QVAPOR'] * (1.-0.622) + 0.622))/(611.2 * numpy.exp(17.67 * (TK-273.15)/(TK-29.65)))
                if factor_name == 'Pressure':
                    ds_target[:]  = (ds['P'] + ds['PB']) / 100
                if factor_name == 'PBLH':
                    ds_target[:]  = ds['PBLH']
                if factor_name == 'WS':
                    ds = wrf.getvar(nc_ds, 'wspd_wdir', units = "m s-1", timeidx = wrf.ALL_TIMES)
                    ds_target[:]  = wrf.to_np(ds[0,hour_start:-1,int(layer)-1,:])
                if factor_name == 'WD':
                    ds = wrf.getvar(nc_ds, 'wspd_wdir', units = "m s-1", timeidx = wrf.ALL_TIMES)
                    ds_target[:]  = wrf.to_np(ds[1,hour_start:-1,int(layer)-1,:])  
                nc_ds.close()
                nc_ds = None
            if ('-' in str(source_index)):  #差值
                file_path1 = self.source_dirs.split(';')[int(source_index.split('-')[0])]
                file_path2 = self.source_dirs.split(';')[int(source_index.split('-')[1])]
                nc_file1 = file_path1 + '/wrfout_d01_' + (current_day + datetime.timedelta(days = -1)).strftime('%Y-%m-%d_12:00:00')
                nc_ds1 = netCDF4.Dataset(nc_file1)
                nc_file2 = file_path2 + '/wrfout_d01_' + (current_day + datetime.timedelta(days = -1)).strftime('%Y-%m-%d_12:00:00')
                nc_ds2 = netCDF4.Dataset(nc_file2)
                if factor_name in ['rh','tc','pressure']:
                    ds1 = wrf.getvar(nc_ds1, factor_name, timeidx = wrf.ALL_TIMES)
                    ds2 = wrf.getvar(nc_ds2, factor_name, timeidx = wrf.ALL_TIMES)
                    ds_target[:]  = wrf.to_np(ds1[hour_start:-1,int(layer)-1,:] - ds2[hour_start:-1,int(layer)-1,:])
                elif factor_name in ['PBLH']:
                    ds1 = wrf.getvar(nc_ds1, factor_name, timeidx = wrf.ALL_TIMES)
                    ds2 = wrf.getvar(nc_ds2, factor_name, timeidx = wrf.ALL_TIMES)
                    ds_target[:]  = wrf.to_np(ds1[hour_start:-1,:] - ds2[hour_start:-1,:])
                elif factor_name in ['ws','wd']:
                    ds1 = wrf.getvar(nc_ds1, 'wspd_wdir', units = "m s-1", timeidx = wrf.ALL_TIMES)
                    ds2 = wrf.getvar(nc_ds2, 'wspd_wdir', units = "m s-1", timeidx = wrf.ALL_TIMES)
                    if factor_name == 'ws':
                        ds_target[:]  = wrf.to_np(ds1[0,hour_start:-1,int(layer)-1,:] - ds2[0,hour_start:-1,int(layer)-1,:])
                    if factor_name == 'wd':
                        ds_target[:]  = wrf.to_np(ds1[1,hour_start:-1,int(layer)-1,:] - ds2[1,hour_start:-1,int(layer)-1,:])
                nc_ds1.close() 
                nc_ds2.close()
                nc_ds1 = None
                nc_ds2 = None

        if self.data_type == 'CMAQ_regularsites' or self.data_type == 'CMAQ_supersites':
            aconc_file = file_path + '/CCTM_ACONC_' + current_day.strftime('%Y-%m-%d') + '.nc'
            if not os.path.exists(aconc_file):
                aconc_file = file_path + '/ACONC_' + current_day.strftime('%Y-%m-%d') + '.nc' 
            apm_file = file_path + '/CCTM_APMDIAG_' + current_day.strftime('%Y-%m-%d') + '.nc'
            if not os.path.exists(apm_file): # CMAQ AE5的粒径分布文件名字较AE6有变化
                apm_file = file_path + '/AERODIAM_' + current_day.strftime('%Y-%m-%d') + '.nc'
            airdensity_file = file_path + '/METCRO3D_' + current_day.strftime('%Y-%m-%d') + '.nc'
            if not os.path.exists(airdensity_file):
                print('metcro3d file was not found in the current folder!!!')
                sys.exit()
            nc_ds = netCDF4.Dataset(aconc_file)
            Domain_NCOLS = int(nc_ds.getncattr('NCOLS'))
            Domain_NROWS = int(nc_ds.getncattr('NROWS'))
            ds_target = numpy.zeros([24, Domain_NROWS, Domain_NCOLS])
            ds_aconc = xarray.open_dataset(aconc_file).sel(LAY=int(layer)-1)
            ds_density = xarray.open_dataset(airdensity_file).sel(LAY=int(layer)-1)
            ds_apm = xarray.open_dataset(apm_file).sel(LAY=0)
            if factor_name == 'CO':
                ds_target[:] = ds_aconc['CO'][:] * 28 * ds_density['DENS'][1-24,:] / 29
            if factor_name == 'NO':
                ds_target[:] = 1000.0 * (ds_aconc['NO'][:] * 30 * ds_density['DENS'][1-24,:] / 29) 
            if factor_name == 'NO2':
                ds_target[:] = 1000.0 * (ds_aconc['NO2'][:] * 46 * ds_density['DENS'][1-24,:] / 29) 
            if factor_name == 'O3':
                ds_target[:] = 1000.0 * (ds_aconc['O3'][:] * 48 * ds_density['DENS'][1-24,:] / 29)
            if factor_name == 'SO2':
                ds_target[:] = 1000.0 * (ds_aconc['SO2'][:] * 64 * ds_density['DENS'][1-24,:] / 29) 
            if factor_name == 'HNO3':
                ds_target[:] = 1000.0 * (ds_aconc['HNO3'][:] * 63 * ds_density['DENS'][1-24,:] / 29) 
            if factor_name == 'NH3':
                ds_target[:] =  1000.0 * (ds_aconc['NH3'][:] * 17 * ds_density['DENS'][1-24,:] / 29)
            if factor_name == 'TS':
                ds_target[:] =  1000.0 * (ds_aconc['SO2'][:] * 32 * ds_density['DENS'][1-24,:] / 29)  \
                                + ds_aconc['ASO4I'][:] * 1 / 3 + ds_aconc['ASO4J'][:] * 1 / 3
            if factor_name == 'TNN':
                ds_target[:] =  1000.0 * (ds_aconc['NO'][:] * 14 * ds_density['DENS'][1-24,:] / 29)  \
                                + 1000.0 * (ds_aconc['NO2'][:] * 14 * ds_density['DENS'][1-24,:] / 29)  \
                                + 1000.0 * (ds_aconc['HNO3'][:] * 14 * ds_density['DENS'][1-24,:] / 29)  \
                                + ds_aconc['ANO3I'][:] * 14 / 62 + ds_aconc['ANO3J'][:] * 14 / 62 
            if factor_name == 'TAN':
                ds_target[:] =  1000.0 * (ds_aconc['NH3'][:] * 14 * ds_density['DENS'][1-24,:] / 29) + ds_aconc['ANH4I'][:] * 14 / 18 + ds_aconc['ANH4J'][:] * 14 / 18
            if 'PM' in factor_name:
                AE_version ='AE6'
                if 'AERODIAM' in apm_file: AE_version ='AE5'
                ds_target[:] = self.extract_pm(factor_name, ds_aconc, ds_density, ds_apm, AE_version)            
            
            if ('-' in str(source_index)):   #差值
                aconc_file1 = file_path1 + '/CCTM_ACONC_' + current_day.strftime('%Y-%m-%d') + '.nc'
                aconc_file2 = file_path2 + '/CCTM_ACONC_' + current_day.strftime('%Y-%m-%d') + '.nc'
                if not os.path.exists(aconc_file1):
                  aconc_file1 = file_path1 + '/ACONC_' + current_day.strftime('%Y-%m-%d') + '.nc'
                if not os.path.exists(aconc_file2):
                  aconc_file2 = file_path2 + '/ACONC_' + current_day.strftime('%Y-%m-%d') + '.nc'
                apm_file1 = file_path1 + '/CCTM_APMDIAG_' + current_day.strftime('%Y-%m-%d') + '.nc'
                apm_file2 = file_path2 + '/CCTM_APMDIAG_' + current_day.strftime('%Y-%m-%d') + '.nc'
                if not os.path.exists(apm_file1):
                   apm_file1 = file_path1 + '/AERODIAM_' + current_day.strftime('%Y-%m-%d') + '.nc'
                if not os.path.exists(apm_file2):
                   apm_file2 = file_path2 + '/AERODIAM_' + current_day.strftime('%Y-%m-%d') + '.nc'
                airdensity_file1 = file_path1 + '/METCRO3D_' + current_day.strftime('%Y-%m-%d') + '.nc'
                airdensity_file2 = file_path2 + '/METCRO3D_' + current_day.strftime('%Y-%m-%d') + '.nc'
                ds_aconc1 = xarray.open_dataset(aconc_file1).sel(LAY=int(layer)-1)
                ds_aconc2 = xarray.open_dataset(aconc_file2).sel(LAY=int(layer)-1)
                ds_apm1 = xarray.open_dataset(apm_file1).sel(LAY=0)
                ds_apm2 = xarray.open_dataset(apm_file2).sel(LAY=0)
                ds_density1 = xarray.open_dataset(airdensity_file1).sel(LAY=int(layer)-1)
                ds_density2 = xarray.open_dataset(airdensity_file2).sel(LAY=int(layer)-1)
                if factor_name == 'CO':
                    ds_target[:] = ds_aconc1['CO'][:] * 28 * ds_density1['DENS'][1-24,:] / 29 \
                            - ds_aconc2['CO'][:] * 28 * ds_density2['DENS'][1-24,:] / 29
                if factor_name == 'NO':
                    ds_target[:] = 1000.0 * (ds_aconc1['NO'][:] * 30 * ds_density1['DENS'][1-24,:] / 29) \
                            - 1000.0 * (ds_aconc2['NO'][:] * 30 * ds_density2['DENS'][1-24,:] / 29)
                if factor_name == 'NO2':
                    ds_target[:] = 1000.0 * (ds_aconc1['NO2'][:] * 46 * ds_density1['DENS'][1-24,:] / 29) \
                            - 1000.0 * (ds_aconc2['NO2'][:] * 46 * ds_density2['DENS'][1-24,:] / 29)
                if factor_name == 'O3':
                    ds_target[:] = 1000.0 * (ds_aconc1['O3'][:] * 48 * ds_density1['DENS'][1-24,:] / 29)\
                            - 1000.0 * (ds_aconc2['O3'][:] * 48 * ds_density2['DENS'][1-24,:] / 29) 
                if factor_name == 'SO2':
                    ds_target[:] = 1000.0 * (ds_aconc1['SO2'][:] * 64 * ds_density1['DENS'][1-24,:] / 29) \
                            - 1000.0 * (ds_aconc2['SO2'][:] * 64 * ds_density2['DENS'][1-24,:] / 29)
                if factor_name == 'HNO3':
                    ds_target[:] = 1000.0 * (ds_aconc1['HNO3'][:] * 63 * ds_density1['DENS'][1-24,:] / 29) \
                            - 1000.0 * (ds_aconc2['HNO3'][:] * 63 * ds_density2['DENS'][1-24,:] / 29)
                if factor_name == 'NH3':
                    ds_target[:] = 1000.0 * (ds_aconc1['NH3'][:] * 17 * ds_density1['DENS'][1-24,:] / 29) \
                            - 1000.0 * (ds_aconc2['NH3'][:] * 17 * ds_density2['DENS'][1-24,:] / 29)
                if 'PM' in factor_name:
                    AE_version1 ='AE6'
                    AE_version2 ='AE6'
                    if 'AERODIAM' in apm_file1: AE_version1 ='AE5'
                    if 'AERODIAM' in apm_file2: AE_version2 ='AE5'
                    ds_target[:] = self.extract_pm(factor_name, ds_aconc1, ds_density1, ds_apm1, AE_version1) \
                        - self.extract_pm(factor_name, ds_aconc2, ds_density2, ds_apm2, AE_version2)

        if self.data_type.lower() == 'adjoint_sensitivity':
            if ('-' not in str(source_index)):  #绝对值            
                adjoint_file = self.output_dir + '/sensitivity_em_' + current_day.strftime('%Y-%m-%d') + '.nc'
                ds_adjoint = xarray.open_dataset(adjoint_file).sel(LAY=0)
                if factor_name == 'VOCs':
                    ds_target[:] = ds_adjoint['MEOH'][:] + ds_adjoint['OLE'][:] + ds_adjoint['TERP'][:] + ds_adjoint['TOL'][:] + ds_adjoint['ETHA'][:]\
                        + ds_adjoint['ETOH'][:] + ds_adjoint['ETH'][:] + ds_adjoint['FORM'][:] + ds_adjoint['ALDX'][:] + ds_adjoint['IOLE'][:] + ds_adjoint['PAR'][:]\
                        + ds_adjoint['XYL'][:] + ds_adjoint['ISOP'][:] + ds_adjoint['ALD2'][:]
                else:
                    ds_target[:] = ds_adjoint[factor_name][:]  

        if self.data_type.lower() == 'adjoint_semi_normalized_sensitivity':
            if ('-' not in str(source_index)):  #绝对值                          
                adjoint_file = self.output_dir + '/semi_normalized_sensitivity_em_' + current_day.strftime('%Y-%m-%d') + '.nc'
                ds_adjoint = xarray.open_dataset(adjoint_file).sel(LAY=0)
                if factor_name == 'VOCs':
                    ds_target[:] = ds_adjoint['MEOH'][:] + ds_adjoint['OLE'][:] + ds_adjoint['TERP'][:] + ds_adjoint['TOL'][:] + ds_adjoint['ETHA'][:]\
                        + ds_adjoint['ETOH'][:] + ds_adjoint['ETH'][:] + ds_adjoint['FORM'][:] + ds_adjoint['ALDX'][:] + ds_adjoint['IOLE'][:] + ds_adjoint['PAR'][:]\
                        + ds_adjoint['XYL'][:] + ds_adjoint['ISOP'][:] + ds_adjoint['ALD2'][:]
                else:
                    ds_target[:] = ds_adjoint[factor_name][:]
        return ds_target
    
    def extract_pm(self, factor_name, ds_aconc, ds_density, ds_apm, AE_version): # Reconstruct the PM2.5 and PM10 mass
        if AE_version =='AE5':  
            ds_apm = ds_apm[['PM25AT','PM25AC','PM25CO']]
        if AE_version =='AE6':
            ds_apm = ds_apm[['PM25AT','PM25AC','PM25CO','PM10AT','PM10AC','PM10CO']] 
        if factor_name in ['PM25_Cl', 'PM25', 'PM10']:
            PM25_Cl = ds_aconc['ACLI'][:] * ds_apm['PM25AT'][:] + ds_aconc['ACLJ'][:] * ds_apm['PM25AC'][:] \
                    + ds_aconc['ACLK'][:] * ds_apm['PM25CO'][:]
        if factor_name in ['PM25_Na', 'PM25', 'PM10']:              
            if AE_version =='AE5': 
                PM25_Na = 0.78*(ds_aconc['ANAJ'][:] * ds_apm['PM25AC'][:] + ds_aconc['ANAK'][:] * ds_apm['PM25CO'][:])                            
            if AE_version =='AE6':
                PM25_Na = ds_aconc['ANAI'][:] * ds_apm['PM25AT'][:] + ds_aconc['ANAJ'][:] * ds_apm['PM25AC'][:] \
                    + (0.8373 * ds_aconc['ASEACAT'][:] + 0.0626 * ds_aconc['ASOIL'][:] \
                    + 0.0023 * ds_aconc['ACORS'][:]) * ds_apm['PM25CO'][:]
        if factor_name in ['PM25_Mg', 'PM25', 'PM10']:  
            if AE_version =='AE5': 
                PM25_Mg = 0                            
            if AE_version =='AE6':          
                PM25_Mg = ds_aconc['AMGJ'][:] * ds_apm['PM25AC'][:] + (0.0997 * ds_aconc['ASEACAT'][:] \
                    + 0.0170 * ds_aconc['ASOIL'][:] + 0.0032 * ds_aconc['ACORS'][:]) * ds_apm['PM25CO'][:]                
        if factor_name in ['PM25_K', 'PM25', 'PM10']:                                 
            if AE_version =='AE5': 
                PM25_Mg = 0 
            if AE_version =='AE6': 
                PM25_K = ds_aconc['AKJ'][:] * ds_apm['PM25AC'][:] + (0.0310 * ds_aconc['ASEACAT'][:] \
                    + 0.0242 * ds_aconc['ASOIL'][:] + 0.0176 * ds_aconc['ACORS'][:]) * ds_apm['PM25CO'][:]                
        if factor_name in ['PM25_Ca', 'PM25', 'PM10']: 
            if AE_version =='AE5': 
                PM25_Ca = 0
            if AE_version =='AE6':
                PM25_Ca = ds_aconc['ACAJ'][:] * ds_apm['PM25AC'][:] + (0.0320 * ds_aconc['ASEACAT'][:] \
                    + 0.0838 * ds_aconc['ASOIL'][:] + 0.0562 * ds_aconc['ACORS'][:]) * ds_apm['PM25CO'][:]                
        if factor_name in ['PM25_NH4', 'PM25', 'PM10']:
            PM25_NH4 = ds_aconc['ANH4I'][:] * ds_apm['PM25AT'][:] + ds_aconc['ANH4J'][:] * ds_apm['PM25AC'][:] \
                + ds_aconc['ANH4K'][:] * ds_apm['PM25CO'][:]
        if factor_name in ['PM25_NO3', 'PM25', 'PM10']:
            PM25_NO3 = ds_aconc['ANO3I'][:] * ds_apm['PM25AT'][:] + ds_aconc['ANO3J'][:] * ds_apm['PM25AC'][:] \
                + ds_aconc['ANO3K'][:] * ds_apm['PM25CO'][:]
        if factor_name in ['PM25_SO4', 'PM25', 'PM10']:
            PM25_SO4 = ds_aconc['ASO4I'][:] * ds_apm['PM25AT'][:] + ds_aconc['ASO4J'][:] * ds_apm['PM25AC'][:] \
                + ds_aconc['ASO4K'][:] * ds_apm['PM25CO'][:]
        if factor_name in ['PM25_EC', 'PM25', 'PM10']:
            PM25_EC = ds_aconc['AECI'][:]*ds_apm['PM25AT'][:] + ds_aconc['AECJ'][:]*ds_apm['PM25AC'][:]                
        if factor_name in ['PM25_SOIL', 'PM25', 'PM10']:
            if AE_version =='AE5': 
                PM25_SOIL = 0
            if AE_version =='AE6':                        
                PM25_SOIL = (2.20 * ds_aconc['AALJ'][:] +2.49 * ds_aconc['ASIJ'][:] + 1.63 * ds_aconc['ACAJ'][:] \
                    + 2.42 * ds_aconc['AFEJ'][:] + 1.94 * ds_aconc['ATIJ'][:]) * ds_apm['PM25AC'][:] \
                    + ds_aconc['ASOIL'][:] * ds_apm['PM25CO'][:]
        if factor_name in ['PM25_POC', 'PM25_OC']:
            if AE_version =='AE5': 
                PM25_POC = 0
            if AE_version =='AE6':
                APOCI =  ds_aconc['ALVPO1I'][:]/1.39 + ds_aconc['ASVPO1I'][:]/1.32 + ds_aconc['ASVPO2I'][:]/1.26 \
                    + ds_aconc['APOCI'][:]
                APOCJ =  ds_aconc['ALVPO1J'][:]/1.39 + ds_aconc['ASVPO1J'][:]/1.32 + ds_aconc['ASVPO2J'][:]/1.26 \
                    + ds_aconc['ASVPO3J'][:]/1.21 + ds_aconc['AIVPO1J'][:]/1.17 +  ds_aconc['APOCJ'][:]
                PM25_POC = APOCI * ds_apm['PM25AT'][:] + APOCJ * ds_apm['PM25AC'][:]                
        if factor_name  in ['PM25_SOC', 'PM25_OC', 'PM25', 'PM10']:
            if AE_version =='AE5': 
                PM25_SOC = 0
            if AE_version =='AE6':
                ASOCI = ds_aconc['ALVOO1I'][:]/2.27 + ds_aconc['ALVOO2I'][:]/2.06 + ds_aconc['ASVOO1I'][:]/1.88 \
                    + ds_aconc['ASVOO2I'][:]/1.73
                if  'ALK1' in ds_aconc:  # saprc mech.
                    ASOCJ = ds_aconc['AISO1J'][:]/2.20 + ds_aconc['AISO2J'][:]/2.23 + ds_aconc['AISO3J'][:]/2.80 \
                        + ds_aconc['ASQTJ'][:]/1.52 + ds_aconc['AORGCJ'][:]/2.00 + ds_aconc['AOLGBJ'][:]/2.10 \
                        + ds_aconc['AOLGAJ'][:]/2.50 + ds_aconc['AIETETJ'][:]/2.27 + ds_aconc['AIEOSJ'][:]/3.6 \
                        + ds_aconc['ADIMJ'][:]/2.07 + ds_aconc['AIMGAJ'][:]/2.5 + ds_aconc['AIMOSJ'][:]/4.17 \
                        + ds_aconc['AMT1J'][:]/1.67 + ds_aconc['AMT2J'][:]/1.67 + ds_aconc['AMT3J'][:]/1.72 \
                        + ds_aconc['AMT4J'][:]/1.53 + ds_aconc['AMT5J'][:]/1.57 + ds_aconc['AMT6J'][:]/1.40 \
                        + ds_aconc['AMTNO3J'][:]/1.90 + ds_aconc['AMTHYDJ'][:]/1.54 + ds_aconc['AISOPNNJ'][:]/3.8 \
                        + ds_aconc['AGLYJ'][:]/2.13 + ds_aconc['ALVOO1J'][:]/2.27 + ds_aconc['ALVOO2J'][:]/2.06 \
                        + ds_aconc['ASVOO1J'][:]/1.88 + ds_aconc['ASVOO2J'][:]/1.73 + ds_aconc['ASVOO3J'][:]/1.60 \
                        + ds_aconc['APCSOJ'][:]/2.00 + ds_aconc['AAVB1J'][:]/2.70 + ds_aconc['AAVB2J'][:]/2.35 \
                        + ds_aconc['AAVB3J'][:]/2.17 + ds_aconc['AAVB4J'][:]/1.99
                elif 'IOLE' in ds_aconc:  # carbon mech. 
                    ASOCJ = ds_aconc['AISO1J'][:]/2.20 + ds_aconc['AISO2J'][:]/2.23 + ds_aconc['AISO3J'][:]/2.80 \
                        + ds_aconc['ASQTJ'][:]/1.52 + ds_aconc['AORGCJ'][:]/2.00 + ds_aconc['AOLGBJ'][:]/2.10 \
                        + ds_aconc['AOLGAJ'][:]/2.50 + ds_aconc['AMT1J'][:]/1.67 + ds_aconc['AMT2J'][:]/1.67 \
                        + ds_aconc['AMT3J'][:]/1.72 + ds_aconc['AMT4J'][:]/1.53 + ds_aconc['AMT5J'][:]/1.57 \
                        + ds_aconc['AMT6J'][:]/1.40 + ds_aconc['AMTNO3J'][:]/1.90 + ds_aconc['AMTHYDJ'][:]/1.54 \
                        + ds_aconc['AGLYJ'][:]/2.13 + ds_aconc['ALVOO1J'][:]/2.27 + ds_aconc['ALVOO2J'][:]/2.06 \
                        + ds_aconc['ASVOO1J'][:]/1.88 + ds_aconc['ASVOO2J'][:]/1.73 + ds_aconc['ASVOO3J'][:]/1.60 \
                        + ds_aconc['APCSOJ'][:]/2.00 + ds_aconc['AAVB1J'][:]/2.70 + ds_aconc['AAVB2J'][:]/2.35 \
                        + ds_aconc['AAVB3J'][:]/2.17 + ds_aconc['AAVB4J'][:]/1.99                              
                PM25_SOC = ASOCI * ds_apm['PM25AT'][:] + ASOCJ * ds_apm['PM25AC'][:]                
        if factor_name in ['PM25_OC']:
            if AE_version =='AE5': 
                AOCT = (ds_aconc['AXYL1J'][:]+ds_aconc['AXYL2J'][:]+ds_aconc['AXYL3J'][:])/2.0 \
                        +(ds_aconc['ATOL1J'][:]+ds_aconc['ATOL2J'][:]+ds_aconc['ATOL3J'][:])/2.0 \
                        +(ds_aconc['ABNZ1J'][:]+ds_aconc['ABNZ2J'][:]+ds_aconc['ABNZ3J'][:])/2.0 \
                        +(ds_aconc['AISO1J'][:]+ds_aconc['AISO2J'][:])/1.6+ds_aconc['AISO3J'][:]/2.7+(ds_aconc['ATRP1J'][:]+ds_aconc['ATRP2J'][:])/1.4+ds_aconc['ASQTJ'][:]/2.1 \
                        +0.64*ds_aconc['AALKJ'][:]+ds_aconc['AORGCJ'][:]/2.0+(ds_aconc['AOLGBJ'][:]+ds_aconc['AOLGAJ'][:])/2.1+ds_aconc['AORGPAI'][:]+ds_aconc['AORGPAJ'][:]
                PM25_OC = ds_aconc['AORGPAI'][:]*ds_apm['PM25AT'][:] + (AOCT - ds_aconc['AORGPAI'][:])*ds_apm['PM25AC'][:]
            if AE_version =='AE6':
                PM25_OC = (APOCI+ASOCI)*ds_apm['PM25AT'][:] + (APOCJ+ASOCJ)*ds_apm['PM25AC'][:]                
        if factor_name in ['PM25_POM', 'PM25_OM', 'PM25', 'PM10']:
            if AE_version =='AE5': 
                PM25_POM = 0
            if AE_version =='AE6':                        
                APOMI = ds_aconc['ALVPO1I'][:] + ds_aconc['ASVPO1I'][:] + ds_aconc['ASVPO2I'][:] + ds_aconc['APOCI'][:] \
                    + ds_aconc['APNCOMI'][:]
                APOMJ = ds_aconc['ALVPO1J'][:] + ds_aconc['ASVPO1J'][:] + ds_aconc['ASVPO2J'][:] + ds_aconc['APOCJ'][:] \
                    + ds_aconc['ASVPO3J'][:] + ds_aconc['AIVPO1J'][:] + ds_aconc['APNCOMJ'][:]                        
                PM25_POM = APOMI * ds_apm['PM25AT'][:] + APOMJ * ds_apm['PM25AC'][:]                        
        if factor_name in ['PM25_SOM', 'PM25_OM', 'PM25', 'PM10']:
            if AE_version =='AE5':  
                PM25_SOM = 0
            if AE_version =='AE6':                          
                ASOMI = ds_aconc['ALVOO1I'][:] + ds_aconc['ALVOO2I'][:] + ds_aconc['ASVOO1I'][:] + ds_aconc['ASVOO2I'][:]
                if  'ALK1' in ds_aconc:  # saprc mech.
                    ASOMJ = ds_aconc['AISO1J'][:] + ds_aconc['AISO2J'][:] + ds_aconc['AISO3J'][:] + ds_aconc['ASQTJ'][:] \
                        + ds_aconc['AORGCJ'][:] + ds_aconc['AOLGBJ'][:] + ds_aconc['AOLGAJ'][:] + ds_aconc['AIETETJ'][:] \
                        + ds_aconc['AIEOSJ'][:] + ds_aconc['ADIMJ'][:] + ds_aconc['AIMGAJ'][:] + ds_aconc['AIMOSJ'][:] \
                        + ds_aconc['AMT1J'][:] + ds_aconc['AMT2J'][:] + ds_aconc['AMT3J'][:] + ds_aconc['AMT4J'][:] \
                        + ds_aconc['AMT5J'][:] + ds_aconc['AMT6J'][:] + ds_aconc['AMTNO3J'][:] + ds_aconc['AMTHYDJ'][:] \
                        + ds_aconc['AGLYJ'][:] + ds_aconc['AISOPNNJ'][:] + ds_aconc['ALVOO1J'][:] + ds_aconc['ALVOO2J'][:] \
                        + ds_aconc['ASVOO1J'][:] + ds_aconc['ASVOO2J'][:] + ds_aconc['ASVOO3J'][:] + ds_aconc['APCSOJ'][:] \
                        + ds_aconc['AAVB1J'][:] + ds_aconc['AAVB2J'][:] + ds_aconc['AAVB3J'][:] + ds_aconc['AAVB4J'][:]
                elif 'IOLE' in ds_aconc:  # carbon mech. 
                    ASOMJ = ds_aconc['AISO1J'][:] + ds_aconc['AISO2J'][:] + ds_aconc['AISO3J'][:] + ds_aconc['ASQTJ'][:] \
                        + ds_aconc['AORGCJ'][:] + ds_aconc['AOLGBJ'][:] + ds_aconc['AOLGAJ'][:] + ds_aconc['AMT1J'][:] \
                        + ds_aconc['AMT2J'][:] + ds_aconc['AMT3J'][:] + ds_aconc['AMT4J'][:] + ds_aconc['AMT5J'][:] \
                        + ds_aconc['AMT6J'][:] + ds_aconc['AMTNO3J'][:] + ds_aconc['AMTHYDJ'][:] + ds_aconc['AGLYJ'][:] \
                        + ds_aconc['ALVOO1J'][:] + ds_aconc['ALVOO2J'][:] + ds_aconc['ASVOO1J'][:] + ds_aconc['ASVOO2J'][:] \
                        + ds_aconc['ASVOO3J'][:] + ds_aconc['APCSOJ'][:] + ds_aconc['AAVB1J'][:] + ds_aconc['AAVB2J'][:] \
                        + ds_aconc['AAVB3J'][:] + ds_aconc['AAVB4J'][:]
                PM25_SOM = ASOMI * ds_apm['PM25AT'][:] + ASOMJ * ds_apm['PM25AC'][:]
        if factor_name in ['PM25_OM']:
            if AE_version =='AE5':  
                AOCT = (ds_aconc['AXYL1J'][:]+ds_aconc['AXYL2J'][:]+ds_aconc['AXYL3J'][:])/2.0 \
                        +(ds_aconc['ATOL1J'][:]+ds_aconc['ATOL2J'][:]+ds_aconc['ATOL3J'][:])/2.0 \
                        +(ds_aconc['ABNZ1J'][:]+ds_aconc['ABNZ2J'][:]+ds_aconc['ABNZ3J'][:])/2.0 \
                        +(ds_aconc['AISO1J'][:]+ds_aconc['AISO2J'][:])/1.6+ds_aconc['AISO3J'][:]/2.7+(ds_aconc['ATRP1J'][:]+ds_aconc['ATRP2J'][:])/1.4+ds_aconc['ASQTJ'][:]/2.1 \
                        +0.64*ds_aconc['AALKJ'][:]+ds_aconc['AORGCJ'][:]/2.0+(ds_aconc['AOLGBJ'][:]+ds_aconc['AOLGAJ'][:])/2.1+ds_aconc['AORGPAI'][:]+ds_aconc['AORGPAJ'][:]
                PM25_OM = 1.4* (ds_aconc['AORGPAI'][:]*ds_apm['PM25AT'][:] + (AOCT - ds_aconc['AORGPAI'][:])*ds_apm['PM25AC'][:])
            if AE_version =='AE6':
                PM25_OM = (APOMI+ASOMI)*ds_apm['PM25AT'][:] + (APOMJ+ASOMJ)*ds_apm['PM25AC'][:]     
        if factor_name in ['PM25', 'PM10']:
            if AE_version =='AE5':
                AORGAT = ds_aconc['AXYL1J'][:]+ds_aconc['AXYL2J'][:]+ds_aconc['AXYL3J'][:]+ds_aconc['ATOL1J'][:]\
                    + ds_aconc['ATOL2J'][:]+ds_aconc['ATOL3J'][:]+ds_aconc['ABNZ1J'][:]+ds_aconc['ABNZ2J'][:]\
                    + ds_aconc['ABNZ3J'][:]+ds_aconc['AALKJ'][:]+ds_aconc['AOLGAJ'][:] 
                AORGBT = ds_aconc['AISO1J'][:]+ds_aconc['AISO2J'][:]+ds_aconc['AISO3J'][:]+ds_aconc['ATRP1J'][:]\
                    + ds_aconc['ATRP2J'][:]+ds_aconc['ASQTJ'][:]+ds_aconc['AOLGBJ'][:]
                PM25 = (ds_aconc['ASO4I'][:] + ds_aconc['ANO3I'][:] + ds_aconc['ANH4I'][:] + \
                    + ds_aconc['ACLI'][:] + ds_aconc['AECI'][:] + 1.4*ds_aconc['AORGPAI'][:]) * ds_apm['PM25AT'][:] \
                    + (ds_aconc['ASO4J'][:] + ds_aconc['ANO3J'][:] + ds_aconc['ANH4J'][:] \
                    + 1.4*ds_aconc['AORGPAJ'][:] + ds_aconc['ANAJ'][:] + ds_aconc['ACLJ'][:] \
                    + ds_aconc['AECJ'][:] + ds_aconc['A25J'][:] - 0.2*ds_aconc['AORGPAJ'][:] \
                    + AORGAT + AORGBT + ds_aconc['AORGCJ'][:]) * ds_apm['PM25AC'][:] \
                    + (ds_aconc['ASOIL'][:] + ds_aconc['ACORS'][:] + ds_aconc['ANAK'][:] + ds_aconc['ACLK'][:] \
                    + ds_aconc['ASO4K'][:] + ds_aconc['ANO3K'][:] + ds_aconc['ANH4K'][:]) * ds_apm['PM25CO'][:]   
            if AE_version =='AE6':
                PM25 = (ds_aconc['ASO4I'][:] + ds_aconc['ANO3I'][:] + ds_aconc['ANH4I'][:] + ds_aconc['ANAI'][:] \
                    + ds_aconc['ACLI'][:] + ds_aconc['AECI'][:] + (APOMI + ASOMI) + ds_aconc['AOTHRI'][:]) * ds_apm['PM25AT'][:] \
                    + (ds_aconc['ASO4J'][:] + ds_aconc['ANO3J'][:] + ds_aconc['ANH4J'][:] + ds_aconc['ANAJ'][:] + ds_aconc['ACLJ'][:] \
                    + ds_aconc['AECJ'][:] + (APOMJ + ASOMJ) + ds_aconc['AOTHRJ'][:] + ds_aconc['AFEJ'][:] \
                    + ds_aconc['ASIJ'][:] + ds_aconc['ATIJ'][:] + ds_aconc['ACAJ'][:] + ds_aconc['AMGJ'][:] \
                    + ds_aconc['AMNJ'][:] + ds_aconc['AALJ'][:] + ds_aconc['AKJ'][:]) * ds_apm['PM25AC'][:] \
                    + (ds_aconc['ASOIL'][:] + ds_aconc['ACORS'][:] + ds_aconc['ASEACAT'][:] + ds_aconc['ACLK'][:] \
                    + ds_aconc['ASO4K'][:] + ds_aconc['ANO3K'][:] + ds_aconc['ANH4K'][:]) * ds_apm['PM25CO'][:]               
        if factor_name == 'PM10':
            if AE_version == 'AE5': 
                PM10 = PM25 + ds_aconc['ASO4I'][:] + ds_aconc['ASO4J'][:] + ds_aconc['ASO4K'][:] - PM25_SO4 \
                    + ds_aconc['ANO3I'][:] + ds_aconc['ANO3J'][:] + ds_aconc['ANO3K'][:] - PM25_NO3 \
                    + ds_aconc['ANH4I'][:] + ds_aconc['ANH4J'][:] + ds_aconc['ANH4K'][:] - PM25_NH4 \
                    + ds_aconc['ACLI'][:] + ds_aconc['ACLJ'][:] + ds_aconc['ACLK'][:] - PM25_Cl \
                    + (ds_aconc['ANAJ'][:] + ds_aconc['ANAK'][:]) * 0.78 - PM25_Na \
                    + (ds_aconc['ASOIL'][:] + ds_aconc['ACORS'][:]+ 0.22*ds_aconc['ANAK']) * (1.0-ds_apm['PM25CO'][:])
            if AE_version =='AE6':
                PM10 = (ds_aconc['ASO4I'][:] + ds_aconc['ANO3I'][:] + ds_aconc['ANH4I'][:] + ds_aconc['ANAI'][:] \
                    + ds_aconc['ACLI'][:] + ds_aconc['AECI'][:] + (APOMI + ASOMI) + ds_aconc['AOTHRI'][:]) * ds_apm['PM10AT'][:] \
                    + (ds_aconc['ASO4J'][:] + ds_aconc['ANO3J'][:] + ds_aconc['ANH4J'][:] + ds_aconc['ANAJ'][:] \
                    + ds_aconc['ACLJ'][:] + ds_aconc['AECJ'][:] + (APOMJ + ASOMJ) + ds_aconc['AOTHRJ'][:] \
                    + ds_aconc['AFEJ'][:] + ds_aconc['ASIJ'][:] + ds_aconc['ATIJ'][:] + ds_aconc['ACAJ'][:] \
                    + ds_aconc['AMGJ'][:] + ds_aconc['AMNJ'][:] + ds_aconc['AALJ'][:] + ds_aconc['AKJ'][:]) * ds_apm['PM10AC'][:] \
                    + (ds_aconc['ASOIL'][:] + ds_aconc['ACORS'][:] + ds_aconc['ASEACAT'][:] + ds_aconc['ACLK'][:] \
                    + ds_aconc['ASO4K'][:] + ds_aconc['ANO3K'][:] + ds_aconc['ANH4K'][:]) * ds_apm['PM10CO'][:]

        return locals()[factor_name]

    def draw_offline_PNG(self, source_index, layer_index, factor_index, day_index, dataset):
        factor_name = self.factors[factor_index].split('|')[0]
        factor_unit = self.factors[factor_index].split('|')[1]
        factor_desc = self.factors[factor_index].split('|')[3]
        layer = self.extracted_layers[layer_index]
        current_day  = self.date_start + datetime.timedelta(days=day_index)
        f = xarray.open_dataset(self.root_dir + "/data/ETOPO2v2g_f4.nc")
        lon = f['x'].sel(x=slice(50, 155)).values
        lat = f['y'].sel(y=slice(0, 70)).values
        self.dem = f['z'].sel(x=slice(50, 155),y=slice(0, 70)).values
        lon_a, lat_a = numpy.meshgrid(lon, lat)
        self.dem_levels = [-8000, -6000, -4000, -2000, -1000, -200, -50, 0, 50, 200, 500, 1000, 1500, 2000, 3000, 4000,5000, 6000, 7000, 8000]
        self.dem_color = ['#084594', '#2171b5', '#4292c6', '#6baed6', '#9ecae1', '#c6dbef', '#deebf7', '#006837', '#31a354', '#78c679', '#addd8e', \
                '#d9f0a3', '#f7fcb9', '#c9bc87', '#a69165', '#856b49', '#664830', '#ad9591', '#d7ccca']
        # plot city's mark
        df = pandas.read_csv(self.root_dir + "/utils/asia_city.csv")
        df_admin = df[df['capital'] == 'admin']
        self.lon_admin = df_admin['lng'].to_list()
        self.lat_admin = df_admin['lat'].to_list()
        if self.data_type == 'WRF':
            nc_ds = netCDF4.Dataset(self.source_dirs.split(';')[source_index] + '/wrfout_d01_' + (self.date_start + datetime.timedelta(days = -1)).strftime('%Y-%m-%d_12:00:00'))
            Domain_NCOLS = int(nc_ds.getncattr('WEST-EAST_PATCH_END_UNSTAG')) # len(nc_ds.dimensions['west_east']) also works
            Domain_NROWS = int(nc_ds.getncattr('SOUTH-NORTH_PATCH_END_UNSTAG')) # len(nc_ds.dimensions['south_north']) also works
            Domain_CEN_LAT = float(nc_ds.getncattr('CEN_LAT'))
            Domain_CEN_LON = float(nc_ds.getncattr('CEN_LON'))
            Domain_TRUELAT1 = float(nc_ds.getncattr('TRUELAT1'))
            Domain_TRUELAT2 = float(nc_ds.getncattr('TRUELAT2'))
            Domain_XORIG = float(wrf.cartopy_xlim(wrfin=nc_ds)[0])
            Domain_YORIG = float(wrf.cartopy_ylim(wrfin=nc_ds)[0])
            Domain_DX = float(nc_ds.getncattr('DX'))
            Domain_DY = float(nc_ds.getncattr('DY'))
        if self.data_type == 'CMAQ_regularsites' or self.data_type == 'CMAQ_supersites': 
            init_file = self.source_dirs.split(';')[source_index] + '/ACONC_' + self.date_start.strftime('%Y-%m-%d') + '.nc'            
            if not os.path.exists(init_file):
                init_file = self.source_dirs.split(';')[source_index] + '/CCTM_ACONC_' + self.date_start.strftime('%Y-%m-%d') + '.nc'
            if not os.path.exists(init_file):
                print('the error of data sources!!')
                sys.exit()
            else:
                nc_ds = netCDF4.Dataset(init_file)
            Domain_CEN_LAT = float(nc_ds.getncattr('YCENT'))
            Domain_CEN_LON = float(nc_ds.getncattr('XCENT'))
            Domain_TRUELAT1 = float(nc_ds.getncattr('P_ALP'))
            Domain_TRUELAT2 = float(nc_ds.getncattr('P_BET'))
            Domain_XORIG = float(nc_ds.getncattr('XORIG'))
            Domain_YORIG = float(nc_ds.getncattr('YORIG'))
            Domain_DX = float(nc_ds.getncattr('XCELL'))
            Domain_DY = float(nc_ds.getncattr('YCELL'))
            Domain_NCOLS = int(nc_ds.getncattr('NCOLS'))
            Domain_NROWS = int(nc_ds.getncattr('NROWS'))
        LCCProj = pyproj.Proj(proj='lcc', lat_1=Domain_TRUELAT1, lat_2=Domain_TRUELAT2,lat_0=Domain_CEN_LAT, lon_0=Domain_CEN_LON, a=6370000, b=6370000)
        lcc_x_admin, lcc_y_admin = LCCProj(df_admin['lng'].to_list(), df_admin['lat'].to_list())
        self.name_admin = df_admin['city'].to_list()
        # plot capital's mark
        df_primary = df[df['capital'] == 'primary']
        self.lon_primary = df_primary['lng'].to_list()
        self.lat_primary = df_primary['lat'].to_list()
        lcc_x_primary, lcc_y_primary = LCCProj(df_primary['lng'].to_list(), df_primary['lat'].to_list())
        self.name_primary = df_primary['city'].to_list()            
        self.minvalue = float(self.factors[factor_index].split('|')[2].split('~')[0])  
        self.maxvalue = float(self.factors[factor_index].split('|')[2].split('~')[1])
        if '-' in str(source_index):  # 如果比较差异，需要调整图例范围为原来的0.1      
            self.minvalue = self.minvalue * 0.1
            self.maxvalue = self.maxvalue * 0.1
        self.conc_color = ['white', 'blue', 'cyan', 'yellow', 'red', 'darkred']
        mycmap = col.LinearSegmentedColormap.from_list('own', self.conc_color, N=len(self.conc_color))
        my_cmap_aplpha = mycmap(numpy.arange(mycmap.N))
        for i in range(my_cmap_aplpha.shape[0]):
            if i < 1:
                my_cmap_aplpha[i,-1] = 0.1
            else:
                my_cmap_aplpha[i,-1] = 0.8
        self.my_cmap_aplpha = col.ListedColormap(my_cmap_aplpha,name='myowncolor')
        if factor_name in ['PBLH','ws'] and ('-' not in str(source_index)): self.my_cmap_aplpha = self.my_cmap_aplpha.reversed()

        lcc_proj = ccrs.LambertConformal(central_longitude=Domain_CEN_LON, central_latitude=Domain_CEN_LAT,\
                    standard_parallels=(Domain_TRUELAT1,Domain_TRUELAT2))
        lcc_extent = [Domain_XORIG, Domain_XORIG + Domain_DX*Domain_NCOLS,Domain_YORIG, Domain_YORIG + Domain_DY*Domain_NROWS]
        
        for hour_index in range(24):
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(1, 1, 1, projection = lcc_proj)
            self.ax.set_extent(lcc_extent, crs=lcc_proj)
            ret = self.ax.projection.transform_points(ccrs.PlateCarree(),numpy.array(lon_a),numpy.array(lat_a))
            self.xx = ret[..., 0]
            self.yy = ret[..., 1]
            self.ax.contourf(self.xx,self.yy,self.dem,levels=self.dem_levels,colors=self.dem_color,extend='both')
            zip_object = zip(lcc_x_admin, lcc_y_admin, self.name_admin)
            for (x, y, ne) in zip_object:
                if(x>Domain_XORIG and x<Domain_XORIG+Domain_DX*Domain_NCOLS and\
                y>Domain_YORIG and y<Domain_YORIG+Domain_DY*Domain_NROWS):
                    self.ax.text(x,y,ne,color='white',va='center',ha='center',zorder=4,fontsize=4)
            zip_object = zip(lcc_x_primary, lcc_y_primary, self.name_primary)
            for (x, y, ne) in zip_object:
                if(x>Domain_XORIG and x<Domain_XORIG+Domain_DX*Domain_NCOLS and\
                y>Domain_YORIG and y<Domain_YORIG+Domain_DY*Domain_NROWS):
                    self.ax.plot(x, y, "ro", zorder=5,  markersize=1)
                    self.ax.text(x+1000,y-1000,ne,va='center',color='k',ha='center', zorder=4,fontsize=5)
            p = self.ax.imshow(dataset[hour_index,:], cmap=self.my_cmap_aplpha, interpolation=None, alpha=1, vmin=self.minvalue, \
                        vmax=self.maxvalue, origin='lower', extent=lcc_extent)
            p.set_zorder(2)
            plt.colorbar(p, ticks=[x for x in numpy.linspace(self.minvalue,self.maxvalue,num=len(self.conc_color)+1)])
            dt_str = (current_day + datetime.timedelta(hours = hour_index + 8))\
                .strftime('%Y-%m-%d %H:00')
            plt.text(0.25,0.97, dt_str + ' (UTC/GMT +08)', fontproperties=self.timesnr_font, fontsize=9,\
                transform = self.ax.transAxes, color='k', bbox=dict(facecolor='w', pad=1, alpha=0.8, edgecolor='none'))
            title_str = 'Spatial distribution of ' + factor_desc + ' for Layer '+str(layer) + '  (' + factor_unit + ')'
            self.ax.set_title(title_str, fontdict={'size':20, 'weight':'bold'}, fontproperties=self.timesnr_font)        
            self.ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
            self.ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
            gl=self.ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.3, alpha=1, \
                x_inline=False, y_inline=False, color='gray')
            gl.top_labels = False                                
            gl.right_labels = False
            gl.xpadding = 8
            gl.xformatter = LONGITUDE_FORMATTER         
            gl.yformatter = LATITUDE_FORMATTER
            self.lon_tickrange = numpy.arange(80,140,10)  # the range of lon to be ticked
            self.lat_tickrange = numpy.arange(10,60,10) # the range of lat to be ticked
                                                    
            gl.xlocator=matplotlib.ticker.FixedLocator(self.lon_tickrange)                       
            gl.ylocator=matplotlib.ticker.FixedLocator(self.lat_tickrange)
            gl.xlabel_style = {'rotation': 0,'ha':'center'}
            gl.ylabel_style = {'rotation': 0}
            self.ax.spines['right'].set_visible(True)
            self.ax.spines['left'].set_visible(True)
            self.ax.spines['bottom'].set_visible(True)
            self.ax.spines['top'].set_visible(True)
            hours_count = (current_day - self.date_start).days*24 + hour_index + 1
            if '-' not in str(source_index):  #绝对值
                source_name = str(self.source_dirs.split(';')[int(source_index)].split('/')[-1])
            if '-' in str(source_index): #差值
                source_name = str(self.source_dirs.split(';')[int(source_index.split('-')[0])].split('/')[-1]) \
                    + '-' + str(self.source_dirs.split(';')[int(source_index.split('-')[1])].split('/')[-1])
            savename = self.output_dir + '/' + factor_name + '_' + source_name \
                +'_Layer' + str(layer) + '_' + str('%04d'%hours_count) + '.png'
            plt.savefig(savename, pad_inches=0, bbox_inches='tight', dpi=600, format='PNG',transparent=True)
            plt.close('all')

    def draw_online_PNG(self, source_index, layer_index, factor_index, day_index, dataset):
        factor_name = self.factors[factor_index].split('|')[0]
        factor_unit = self.factors[factor_index].split('|')[1]
        layer = self.extracted_layers[layer_index]
        current_day  = self.date_start + datetime.timedelta(days=day_index)          
        minvalue = float(self.factors[factor_index].split('|')[2].split('~')[0])  
        maxvalue = float(self.factors[factor_index].split('|')[2].split('~')[1])
        if ('-' in str(source_index)):  # 如果比较差异，需要调整图例范围为原来的0.1    
            minvalue = minvalue * 0.1
            maxvalue = maxvalue * 0.1
        conc_color = ['white', 'blue', 'cyan', 'yellow', 'red', 'darkred']
        mycmap = col.LinearSegmentedColormap.from_list('own', conc_color, N=len(conc_color))
        my_cmap_aplpha = mycmap(numpy.arange(mycmap.N))
        for i in range(my_cmap_aplpha.shape[0]):
            if i < 1:
                my_cmap_aplpha[i,-1] = 0.1
            else:
                my_cmap_aplpha[i,-1] = 0.8
        self.my_cmap_aplpha = col.ListedColormap(my_cmap_aplpha,name='myowncolor')
        if factor_name in ['PBLH','ws'] and (not '-' in str(source_index)): self.my_cmap_aplpha = self.my_cmap_aplpha.reversed()
        if self.data_type == 'WRF':
            nc_ds = netCDF4.Dataset(self.source_dirs.split(';')[source_index] + '/wrfout_d01_' + (self.date_start + datetime.timedelta(days = -1)).strftime('%Y-%m-%d_12:00:00'))
            Domain_NCOLS = int(nc_ds.getncattr('WEST-EAST_PATCH_END_UNSTAG')) # len(nc_ds.dimensions['west_east']) also works
            Domain_NROWS = int(nc_ds.getncattr('SOUTH-NORTH_PATCH_END_UNSTAG')) # len(nc_ds.dimensions['south_north']) also works
            Domain_CEN_LAT = float(nc_ds.getncattr('CEN_LAT'))
            Domain_CEN_LON = float(nc_ds.getncattr('CEN_LON'))
            Domain_TRUELAT1 = float(nc_ds.getncattr('TRUELAT1'))
            Domain_TRUELAT2 = float(nc_ds.getncattr('TRUELAT2'))
            Domain_XORIG = float(wrf.cartopy_xlim(wrfin=nc_ds)[0])
            Domain_YORIG = float(wrf.cartopy_ylim(wrfin=nc_ds)[0])
            Domain_DX = float(nc_ds.getncattr('DX'))
            Domain_DY = float(nc_ds.getncattr('DY'))
        if self.data_type == 'CMAQ_regularsites' or self.data_type == 'CMAQ_supersites':         
            init_file = self.source_dirs.split(';')[source_index] + '/ACONC_' + self.date_start.strftime('%Y-%m-%d') + '.nc'            
            if not os.path.exists(init_file):
                init_file = self.source_dirs.split(';')[source_index] + '/CCTM_ACONC_' + self.date_start.strftime('%Y-%m-%d') + '.nc'
            if not os.path.exists(init_file):
                print('the error of data sources!!')
                sys.exit()
            else:
                nc_ds = netCDF4.Dataset(init_file)
            Domain_CEN_LAT = float(nc_ds.getncattr('YCENT'))
            Domain_CEN_LON = float(nc_ds.getncattr('XCENT'))
            Domain_TRUELAT1 = float(nc_ds.getncattr('P_ALP'))
            Domain_TRUELAT2 = float(nc_ds.getncattr('P_BET'))
            Domain_XORIG = float(nc_ds.getncattr('XORIG'))
            Domain_YORIG = float(nc_ds.getncattr('YORIG'))
            Domain_DX = float(nc_ds.getncattr('XCELL'))
            Domain_DY = float(nc_ds.getncattr('YCELL'))
            Domain_NCOLS = int(nc_ds.getncattr('NCOLS'))
            Domain_NROWS = int(nc_ds.getncattr('NROWS'))
        lcc_proj = ccrs.LambertConformal(central_longitude=Domain_CEN_LON, central_latitude=Domain_CEN_LAT,\
                            standard_parallels=(Domain_TRUELAT1,Domain_TRUELAT2))
        lcc_extent = [Domain_XORIG, Domain_XORIG + Domain_DX*Domain_NCOLS,\
                                    Domain_YORIG, Domain_YORIG + Domain_DY*Domain_NROWS]
        # Create grid points
        LccX, LccY = numpy.meshgrid(numpy.arange(Domain_XORIG + Domain_DX / 2, Domain_XORIG + Domain_NCOLS * Domain_DX, Domain_DX),
                                numpy.arange(Domain_YORIG + Domain_DY/2, Domain_YORIG + Domain_NROWS * Domain_DY, Domain_DY))
        combine_list = numpy.vstack([LccX.flatten(), LccY.flatten()]).T
        # Transform points to Web Mercator
        # Pre-calculate projection transformation function
        transform_func = pyproj.Transformer.from_crs(lcc_proj, "epsg:3857").transform
        web_Mercator_List = numpy.array([list(transform_func(x, y)) for x, y in combine_list])
        x_min = numpy.min(web_Mercator_List[:, 0])
        y_min = numpy.min(web_Mercator_List[:, 1])
        x_max = numpy.max(web_Mercator_List[:, 0])
        y_max = numpy.max(web_Mercator_List[:, 1])
        mercator_extent = [x_min, x_max, y_min, y_max]
        for hour_index in range(24):
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1, frameon=False, projection=ccrs.epsg(3857))
            ax.set_extent(mercator_extent, crs=ccrs.epsg(3857))
            p = ax.imshow(dataset[hour_index,:], extent=lcc_extent, transform=lcc_proj, cmap=self.my_cmap_aplpha,\
                        interpolation=None, alpha=1, vmin=minvalue, vmax=maxvalue, origin='lower')        
            p.set_zorder(2)
            plt.axis('off')
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            dt_str = (current_day + datetime.timedelta(hours = hour_index + 8))\
                .strftime('%Y-%m-%d %H_00')
            if '-' not in str(source_index):  #绝对值
                source_name = str(self.source_dirs.split(';')[int(source_index)].split('/')[-1])
            if '-' in str(source_index): #差值
                source_name = str(self.source_dirs.split(';')[int(source_index.split('-')[0])].split('/')[-1]) \
                    + '-' + str(self.source_dirs.split(';')[int(source_index.split('-')[1])].split('/')[-1])            
            savedir = self.output_dir + '/PNG/online_' + source_name + '_' + factor_name + '_Layer' + str(layer)
            savename = savedir + '/' +  dt_str + '_lon_min(' + str(round(mercator_extent[0]/20037508.34*180,3))\
                +'),lon_max(' + str(round(mercator_extent[1]/20037508.34*180, 3)) + '),lat_min(' \
                + str(round(180 / math.pi * (2 * math.atan(math.exp(mercator_extent[2]/20037508.34*180 * math.pi / 180)) - math.pi / 2), 3))\
                + '),lat_max(' + str(round(180 / math.pi * (2 * math.atan(math.exp(mercator_extent[3]/20037508.34*180 * math.pi / 180)) - math.pi / 2), 3)) +').png'
            plt.savefig(savename, pad_inches=0, bbox_inches='tight', dpi=600, format='PNG',transparent=True)
            # 输出图例
            fig2, c_map_ax = plt.subplots(figsize=(12, 1.5))  # 设置colorbar的尺寸
            fig2.subplots_adjust(bottom=0.5, top=0.9,left=0.5,right=0.9)  # 调整图的上下边距
            sm = plt.cm.ScalarMappable(cmap=self.my_cmap_aplpha, norm=plt.Normalize(vmin=minvalue, vmax=maxvalue))
            cb = plt.colorbar(sm, cax=c_map_ax, ticks=[x for x in numpy.linspace(minvalue,maxvalue,num=len(conc_color)+1)], orientation = 'horizontal')
            #cb.set_label('label', fontsize=14) 
            plt.title(factor_name + '(' + factor_unit + ')', fontsize=14)
            fig2.savefig(savedir + '/legend.png', bbox_inches='tight', pad_inches=0)
            plt.close('all')

    def submit_spatial(self,num_nodes): # 提交空间图处理任务
        job_tag = 'spatial'
        python_file = self.output_dir + '/' + job_tag + '.py'
        python_job = open(python_file, 'w')
        python_job.write("import os\n")
        python_job.write("from mpi4py import MPI\n")
        python_job.write("import glob\n")
        python_job.write("import sys\n") 
        python_job.write('sys.path.append("'+ self.root_dir + '")\n')
        python_job.write("import post_process\n")
        python_job.write("comm = MPI.COMM_WORLD\n")
        python_job.write("rank = comm.Get_rank()\n")
        python_job.write("size = comm.Get_size()\n")
        python_job.write("if rank == 0:\n") 
        python_job.write("    objProcess = post_process.Post_Process('" + self.data_type \
            + "','" + self.date_start.strftime('%Y-%m-%d') + "','" + self.date_end.strftime('%Y-%m-%d') + "','" + str(self.root_dir) \
            + "','" + str(self.source_dirs) + "'," + str(self.extracted_layers)+ "," + str(self.sites_includedprovinces)+ "," + str(self.sites_includedcities) + ")\n")
        python_job.write("    objProcess.factors = " + str(self.factors) + "\n")
        python_job.write("    objProcess.output_dir = '" + str(self.output_dir) + "'\n")
        python_job.write("    comm.bcast(objProcess, root=0)\n") 
        python_job.write("else:\n")
        python_job.write("    objProcess = comm.bcast(None, root=0)\n")
        python_job.write("total_iterations = len(objProcess.source_dirs.split(';')) * len(objProcess.extracted_layers) * len(objProcess.factors) * ((objProcess.date_end - objProcess.date_start).days + 1)\n")
        python_job.write("iterations_per_process = total_iterations // size\n")
        python_job.write("remainder = total_iterations % size\n")
        python_job.write("start_index = rank * iterations_per_process + min(rank, remainder)\n")
        python_job.write("end_index = start_index + iterations_per_process + (1 if rank < remainder else 0)\n")
        python_job.write("for global_index in range(start_index, end_index):\n")
        python_job.write("    source_index = global_index % len(objProcess.source_dirs.split(';'))\n")
        python_job.write("    layer_index = (global_index // len(objProcess.source_dirs.split(';'))) % len(objProcess.extracted_layers)\n")
        python_job.write("    factor_index = (global_index // (len(objProcess.source_dirs.split(';')) * len(objProcess.extracted_layers))) % len(objProcess.factors)\n")
        python_job.write("    day_index = (global_index // (len(objProcess.source_dirs.split(';')) * len(objProcess.extracted_layers) * len(objProcess.factors))) % ((objProcess.date_end - objProcess.date_start).days + 1)\n")
        python_job.write("    ds_target = objProcess.extract_2dsimu(source_index, layer_index, factor_index, day_index)\n")
        python_job.write("    objProcess.draw_offline_PNG(source_index, layer_index, factor_index, day_index, ds_target)\n")
        python_job.write("    objProcess.draw_online_PNG(source_index, layer_index, factor_index, day_index, ds_target)\n") 
        python_job.write("comm.Barrier()\n")
        python_job.write("total_iterations = len(objProcess.source_dirs.split(';')) * len(objProcess.extracted_layers) * len(objProcess.factors)\n")
        python_job.write("iterations_per_process = total_iterations // size\n")
        python_job.write("remainder = total_iterations % size\n")
        python_job.write("start_index = rank * iterations_per_process + min(rank, remainder)\n")
        python_job.write("end_index = start_index + iterations_per_process + (1 if rank < remainder else 0)\n")
        python_job.write("for global_index in range(start_index, end_index):\n")
        python_job.write("    source_index = global_index % len(objProcess.source_dirs.split(';'))\n")
        python_job.write("    layer_index = (global_index // len(objProcess.source_dirs.split(';'))) % len(objProcess.extracted_layers)\n")
        python_job.write("    factor_index = (global_index // (len(objProcess.source_dirs.split(';')) * len(objProcess.extracted_layers))) % len(objProcess.factors)\n")
        python_job.write("    key_str = objProcess.factors[factor_index].split('|')[0]  + '_' + str(objProcess.source_dirs.split(';')[int(source_index)].split('/')[-1]) + '_Layer' + str(objProcess.extracted_layers[layer_index])\n")
        python_job.write("    make_mp4 = 'ffmpeg -y -r 2 -f image2 -i ' + objProcess.output_dir + '/' + key_str + '_%4d.png -pix_fmt yuv420p  -c:v libx264 -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" '+ objProcess.output_dir + '/MP4' + '/' + key_str + '.mp4'\n")
        python_job.write("    os.system(make_mp4)\n")
        python_job.write("comm.Barrier()\n")
        python_job.write("if rank == 0:\n")
        python_job.write("    for file in os.listdir(objProcess.output_dir):\n")
        python_job.write("        file_str = os.path.join(objProcess.output_dir, file)\n")
        python_job.write("        if  os.path.isfile(file_str) and file_str.endswith('.png'):\n")
        python_job.write("            os.remove(file_str)\n")
        python_job.write("    print('well done!')\n")
        python_job.close()
        os.chmod(python_file, stat.S_IRWXO+stat.S_IRWXG+stat.S_IRWXU)
        slurm_file = self.output_dir + '/' + job_tag + '.slurm'
        slurm_job = open(slurm_file, 'w')
        slurm_job.write("#!/bin/bash\n")
        slurm_job.write("#SBATCH --job-name=spatial_mpi4py\n")
        slurm_job.write("#SBATCH --time=01:00:00\n")     # maximum running time of two hours
        slurm_job.write("#SBATCH --partition=64c512g\n")
        slurm_job.write("#SBATCH --mail-type=end\n")
        slurm_job.write("#SBATCH --mail-user=\n")
        slurm_job.write("#SBATCH --output=" + self.output_dir + "/" + job_tag + ".log\n")
        slurm_job.write("#SBATCH --error=" + self.output_dir + "/" + job_tag + ".log\n")
        slurm_job.write("#SBATCH -N " + str(num_nodes) + " \n")
        slurm_job.write("#SBATCH --exclusive\n")
        slurm_job.write('mpirun -n '+ str(num_nodes*64) + ' python ' + self.output_dir + '/' + job_tag + '.py' \
            + ' >& ' + self.output_dir + '/' + job_tag + '.log\n')
        slurm_job.close()
        os.system('sbatch ' + slurm_file)       

    def submit_sites(self,num_nodes): # 提取站点数据并进行散点图及时序图绘制
        job_tag = self.data_type
        python_file = self.output_dir + '/' + job_tag + '.py'
        python_job = open(python_file, 'w')
        python_job.write("import os\n")
        python_job.write("import numpy\n")
        python_job.write("import pandas\n")
        python_job.write("import datetime\n")
        python_job.write("from mpi4py import MPI\n")
        python_job.write("import sqlalchemy\n")
        python_job.write("import sys\n") 
        python_job.write('sys.path.append("'+ self.root_dir + '")\n')
        python_job.write("import post_process\n")
        if (self.data_type == 'WRF') or (self.data_type == 'CMAQ_regularsites'):
            python_job.write("airdb_engine = sqlalchemy.create_engine('dialect+driver://username:password@host:port/database')\n")
        if (self.data_type == 'CMAQ_supersites'):
            python_job.write("airdb_engine = sqlalchemy.create_engine('dialect+driver://username:password@host:port/database')\n")
        python_job.write("comm = MPI.COMM_WORLD\n")
        python_job.write("rank = comm.Get_rank()\n")
        python_job.write("size = comm.Get_size()\n")
        python_job.write("if rank == 0:\n")
        python_job.write("    objProcess = post_process.Post_Process('" + self.data_type \
            + "','" + self.date_start.strftime('%Y-%m-%d') + "','" + self.date_end.strftime('%Y-%m-%d') + "','" + str(self.root_dir) \
            + "','" + str(self.source_dirs) + "'," + str(self.extracted_layers)+ "," + str(self.sites_includedprovinces)+ "," + str(self.sites_includedcities) + ")\n")
        python_job.write("    objProcess.output_dir = '" + str(self.output_dir) + "'\n")
        python_job.write("    list_obs = []\n")
        python_job.write("    for factor_index in range(len(objProcess.factors)):\n")
        python_job.write("        factor_name = objProcess.factors[factor_index].split('|')[0]\n")   
        if self.data_type == 'WRF':
            python_job.write("        factor_axis = factor_name\n")
            python_job.write("        if factor_name.upper()=='PBLH':\n")
            python_job.write("            list_obs.append(pandas.DataFrame())\n")
            python_job.write("            continue\n")
            python_job.write("        if factor_name.upper()=='TC':\n")
            python_job.write("            factor_axis = 'TEM'\n")
            python_job.write("        if factor_name.upper()=='RH':\n")
            python_job.write("            factor_axis = 'RHU'\n")
            python_job.write("        if factor_name=='Pressure':\n")
            python_job.write("            factor_axis = 'PRS'\n")
            python_job.write("        if factor_name.upper()=='WS':\n")
            python_job.write("            factor_axis = 'CASE WHEN WIN_S_Avg_10mi IS NULL THEN WIN_S_Avg_2mi ELSE WIN_S_Avg_10mi END AS wind_s'\n")
            python_job.write("        if factor_name.upper()=='WD':\n")
            python_job.write("            factor_axis = 'CASE WHEN WIN_D_Avg_10mi IS NULL THEN WIN_D_Avg_2mi ELSE WIN_D_Avg_10mi END AS wind_D'\n")
            python_job.write("        ds_sql = \"SELECT \" + factor_axis + \", station_code,TimePoint \"\
            + \"FROM t_weather_data WHERE station_code IN ('\" + \"','\".join(objProcess.siteIDs) \
            + \"') AND TimePoint BETWEEN '\" + (objProcess.date_start + datetime.timedelta(hours=8)).strftime(\"%Y-%m-%d %H:00\") + \"' AND '\" \
            + (objProcess.date_end + datetime.timedelta(hours=8)).strftime(\"%Y-%m-%d %H:00\") + \"' ORDER BY TimePoint\"\n")
        if self.data_type == 'CMAQ_regularsites':
            python_job.write("        ds_sql = \"SELECT \" + factor_name.lower() + \",station_code,pubtime FROM t_na_station_realtime \" \
            + \"WHERE station_code IN ('\" + \"','\".join(objProcess.siteIDs) + \"') AND pubtime BETWEEN '\" \
            + (objProcess.date_start + datetime.timedelta(hours=8)).strftime(\"%Y-%m-%d %H:00\") + \"' AND '\" + (objProcess.date_end + datetime.timedelta(hours=8)).strftime(\"%Y-%m-%d %H:00\") \
            + \"' ORDER BY pubtime\"\n")
        if self.data_type == 'CMAQ_supersites':
            python_job.write("        factor_axis = factor_name + '_ug_m3'\n")
            python_job.write("        if 'PM25_' in factor_name.upper():\n")
            python_job.write("            factor_axis = factor_name.replace('PM25_','') + '_ug_m3'\n")            
            python_job.write("        if factor_name.upper()=='CO':\n")
            python_job.write("            factor_axis = 'CO_mg_m3'\n")
            python_job.write("        if factor_name.upper()=='TS':\n")
            python_job.write("            factor_axis = '(SO2_ug_m3*1/2+SO4_ug_m3*1/3) AS TS'\n")
            python_job.write("        if factor_name.upper()=='TNN':\n")
            python_job.write("            factor_axis = '(NO_ug_m3*14/(14+16)+NO2_ug_m3*14/(14+32)+HNO3_ug_m3*14/(14+1+48)+NO3_ug_m3*14/(14+48)) AS TNN'\n")
            python_job.write("        if factor_name.upper()=='TAN':\n")
            python_job.write("            factor_axis = '(NH3_ug_m3*14/17+NH4_ug_m3*14/18) AS TAN'\n")
            python_job.write("        ds_sql = \"SELECT \" + factor_axis + \", station_code,obs_time FROM Supersite_PM_GAS \" \
                + \"WHERE station_code IN ('\" + \"','\".join(objProcess.siteIDs) + \"') AND obs_time BETWEEN '\" \
                + (objProcess.date_start + datetime.timedelta(hours=8)).strftime(\"%Y-%m-%d %H:00\") + \"' AND '\" + (objProcess.date_end + datetime.timedelta(hours=8)).strftime(\"%Y-%m-%d %H:00\") \
                + \"' ORDER BY station_code, obs_time\"\n")
        python_job.write("        df_obs = pandas.DataFrame(numpy.array(pandas.read_sql_query(ds_sql, airdb_engine)))\n")
        python_job.write("        df_obs.columns=[factor_name,'SiteCode','TimePoint']\n")
        python_job.write("        df_obs.drop_duplicates(subset=['SiteCode','TimePoint'],keep='first',inplace=True)\n")
        python_job.write("        dtypes = {factor_name: float, 'SiteCode': str}\n")
        python_job.write("        df_obs = df_obs.astype(dtypes)\n")
        python_job.write("        df_obs['TimePoint'] = pandas.to_datetime(df_obs['TimePoint'])\n")
        python_job.write("        list_obs.append(df_obs)\n")
        python_job.write("    comm.bcast(objProcess, root=0)\n") 
        python_job.write("    comm.bcast(list_obs, root=0)\n") 
        python_job.write("else:\n")
        python_job.write("    objProcess = comm.bcast(None, root=0)\n")
        python_job.write("    list_obs = comm.bcast(None, root=0)\n")
        #加载模拟数据
        python_job.write("df_simu_sites = pandas.DataFrame()\n")
        python_job.write("total_iterations = len(objProcess.source_dirs.split(';')) * len(objProcess.factors) * ((objProcess.date_end - objProcess.date_start).days + 1)\n")
        python_job.write("iterations_per_process = total_iterations // size\n")
        python_job.write("remainder = total_iterations % size\n")
        python_job.write("start_index = rank * iterations_per_process + min(rank, remainder)\n")
        python_job.write("end_index = start_index + iterations_per_process + (1 if rank < remainder else 0)\n")
        python_job.write("for global_index in range(start_index, end_index):\n")
        python_job.write("    source_index = global_index % len(objProcess.source_dirs.split(';'))\n")
        python_job.write("    factor_index = (global_index // (len(objProcess.source_dirs.split(';')))) % len(objProcess.factors)\n")
        python_job.write("    day_index = (global_index // (len(objProcess.source_dirs.split(';')) * len(objProcess.factors))) % ((objProcess.date_end - objProcess.date_start).days + 1)\n")
        python_job.write("    ds_target = objProcess.extract_2dsimu(source_index, 0, factor_index, day_index)\n")        
        python_job.write("    if 'df_simu_list' not in locals():\n")
        python_job.write("        df_simu_list = []\n")
        python_job.write("    site_indices = numpy.arange(objProcess.siteIDs.shape[0])\n")
        python_job.write("    factor_values = ds_target[:, objProcess.siterows[source_index][site_indices], objProcess.sitecols[source_index][site_indices]]\n")
        python_job.write("    for hour_index in range(24):\n")
        python_job.write("        current_time = objProcess.date_start + datetime.timedelta(days=day_index, hours=hour_index + 8)\n") # utc 00> +8
        python_job.write("        df_simu_single = pandas.DataFrame({\n")
        python_job.write("              'source_index': source_index,\n")
        python_job.write("              'SiteCode': objProcess.siteIDs[site_indices],\n")
        python_job.write("              'TimePoint': current_time,\n")
        python_job.write("              'factor_index':factor_index,\n")
        python_job.write("              'factor_value': factor_values[hour_index]\n")
        python_job.write("              })\n")
        python_job.write("        df_simu_list.append(df_simu_single)\n")
        python_job.write("    df_simu_sites = pandas.concat(df_simu_list, ignore_index=True)\n")
        python_job.write("df_simu_sites = comm.allgather(df_simu_sites)\n")
        python_job.write("df_simu_sites = pandas.concat(df_simu_sites, axis=0, ignore_index=True)\n")
        # 开始散点图        
        python_job.write("comm.Barrier()\n")
        python_job.write("total_iterations = len(objProcess.source_dirs.split(';')) * len(objProcess.factors)\n")
        python_job.write("iterations_per_process = total_iterations // size\n")
        python_job.write("remainder = total_iterations % size\n")
        python_job.write("start_index = rank * iterations_per_process + min(rank, remainder)\n")
        python_job.write("end_index = start_index + iterations_per_process + (1 if rank < remainder else 0)\n")
        python_job.write("for global_index in range(start_index, end_index):\n")
        python_job.write("    factor_index = global_index  % len(objProcess.factors)\n")
        python_job.write("    source_index = (global_index // (len(objProcess.factors))) % len(objProcess.source_dirs.split(';'))\n")
        python_job.write("    if not list_obs[factor_index].empty:\n")
        python_job.write("        objProcess.Scatter(factor_index,source_index,list_obs[factor_index],df_simu_sites)\n")
        # 开始时序图
        python_job.write("comm.Barrier()\n")
        python_job.write("total_iterations = len(objProcess.factors) * objProcess.siteIDs.shape[0]\n")
        python_job.write("iterations_per_process = total_iterations // size\n")
        python_job.write("remainder = total_iterations % size\n")
        python_job.write("start_index = rank * iterations_per_process + min(rank, remainder)\n")
        python_job.write("end_index = start_index + iterations_per_process + (1 if rank < remainder else 0)\n")
        python_job.write("for global_index in range(start_index, end_index):\n")
        python_job.write("    factor_index = global_index  % (len(objProcess.factors))\n")
        python_job.write("    site_index = (global_index // (len(objProcess.factors))) % (objProcess.siteIDs.shape[0])\n")
        python_job.write("    if not list_obs[factor_index].empty:\n")
        python_job.write("        objProcess.Timeseries(factor_index,site_index,list_obs[factor_index],df_simu_sites)\n")
        python_job.write("comm.Barrier()\n")
        python_job.write("if rank == 0:  print('well done!')\n")
        python_job.write("MPI.Finalize()\n")
        python_job.close()
        os.chmod(python_file, stat.S_IRWXO+stat.S_IRWXG+stat.S_IRWXU)
        slurm_file = self.output_dir + '/' + job_tag + '.slurm'
        slurm_job = open(slurm_file, 'w')
        slurm_job.write("#!/bin/bash\n")
        slurm_job.write("#SBATCH --job-name=" + job_tag + "\n")
        slurm_job.write("#SBATCH --time=01:00:00\n")     # maximum running time of two hours
        slurm_job.write("#SBATCH --partition=64c512g\n")
        slurm_job.write("#SBATCH --mail-type=end\n")
        slurm_job.write("#SBATCH --mail-user=\n")
        slurm_job.write("#SBATCH --output=" + self.output_dir + "/" + job_tag + ".log\n")
        slurm_job.write("#SBATCH --error=" + self.output_dir + "/" + job_tag + ".log\n")
        slurm_job.write("#SBATCH -N " + str(num_nodes) + " \n")
        slurm_job.write("#SBATCH --exclusive\n")
        slurm_job.write('mpirun -n '+ str(num_nodes*64) + ' python ' + self.output_dir + '/' + job_tag + '.py' \
            + ' >& ' + self.output_dir + '/' + job_tag + '.log\n')
        slurm_job.close()
        os.system('sbatch ' + slurm_file)        

    def run(self):
        if len(set(self.source_dirs.split(';')))!=len(self.source_dirs.split(';')):
            print('数据源有重复输入，请检查！')
            sys.exit ()
        print('Post processing start at: ' + datetime.datetime.now().strftime('%m-%d %H:%M'));
        self.output_dir = self.root_dir + '/PostProcess_' + self.data_type + '_' + self.date_start.strftime('%Y-%m-%d') \
            + '~' + self.date_end.strftime('%Y-%m-%d') # 后处理结果保存在家目录下
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir, ignore_errors=True)
        os.makedirs(self.output_dir)
        os.makedirs(self.output_dir + '/MP4')
        os.makedirs(self.output_dir + '/PNG')
        os.makedirs(self.output_dir + '/Scatter')
        os.makedirs(self.output_dir + '/CSV')
        os.makedirs(self.output_dir + '/Timeseries')
        for source_index in range(len(self.source_dirs.split(';'))):
            for factor_index in range(len(self.factors)):
                for layer_index in range(len(self.extracted_layers)):                
                    savedir = self.output_dir + '/PNG/online_' + self.source_dirs.split(';')[int(source_index)].split('/')[-1] \
                        + '_' + self.factors[factor_index].split('|')[0] + '_Layer' + str(self.extracted_layers[layer_index])
                    os.makedirs(savedir)

        self.submit_spatial(2)  # 提交画图任务
        self.submit_sites(2)  # 提交站点图任务

        if self.data_type.lower() == 'adjoint_sensitivity':
            os.makedirs(self.output_dir + '/hourly')
            os.makedirs(self.output_dir + '/daily')
            os.makedirs(self.output_dir + '/period')

        if self.data_type.lower() == 'adjoint_semi_normalized_sensitivity':               
            os.makedirs(self.output_dir + '/hourly')
            os.makedirs(self.output_dir + '/daily')
            os.makedirs(self.output_dir + '/period')

        for source_index in range(len(self.source_dirs.split(';'))): # 针对每个数据源进行处理
            if self.data_type.lower() == 'adjoint_sensitivity':
                self.makefiles_sensitivity_hourly(self.source_dirs.split(';')[int(source_index)])
                self.makefiles_sensitivity_daily()
                self.makefiles_sensitivity_period()
            if self.data_type.lower() == 'adjoint_semi_normalized_sensitivity':
                self.makefiles_semi_normalized_sensitivity_hourly(self.source_dirs.split(';')[int(source_index)])
                self.makefiles_semi_normalized_sensitivity_daily()
                self.makefiles_semi_normalized_sensitivity_period()            

        finish_flag = False
        while not finish_flag:
            finish_flag = True # 所有log文件都包含well done表明工作正常结束
            if len(glob.glob(self.output_dir + '/*.log')) < 1: finish_flag = False
            for file in glob.glob(self.output_dir + '/*.log'):
                if (not 'well done' in open(file, 'r').read()): finish_flag = False
            time.sleep(10)
        print('The post processing ended successfully at ' + datetime.datetime.now().strftime('%m-%d %H:%M'))