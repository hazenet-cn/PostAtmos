import cartopy
import numpy
import pyproj
import geopandas
import pandas as pd
from shapely.geometry import Point
import folium
from folium.features import DivIcon


class domain_definition():
    def __init__(self,domain_name,ref_latitude,ref_longtitude,truelat_first,truelat_second,NCOLS,NROWS,X_CELL,Y_CELL,X_ORIG,Y_ORIG):
        self.projection = "lambert"
        self.domain_name= domain_name
        # Gain the input parameters
        self.ref_latitude = ref_latitude
        self.ref_longtitude = ref_longtitude
        # ref_latitude, ref_longtitude: A real value specifying the latitude and longitude part of a (latitude, longitude) location 
        # whose (i, j) location in the simulation domain is known. For ARW, ref_lon gives the longitude
        # of the center-point of the coarse domain by default (i.e., when ref_x and ref_y are not specified). 
        # West longitudes are negative, and the value of ref_lon should be in the range [-180, 180]. No default value.     
        # REF_X, REF_Y : A real value specifying the i,j part of an (i, j) location whose (latitude,longitude) location
        # in the simulation domain is known. The (i, j) location is always given with respect to the mass-staggered grid,
        # whose dimensions are one less than the dimensions of the unstaggered grid. 
        # Default value is REF_X: (((E_WE-1.)+1.)/2.) = (E_WE/2.);  REF_Y: (((E_SN-1.)+1.)/2.) = (E_SN/2.). 
        self.truelat_first = truelat_first
        self.truelat_second = truelat_second
        self.stand_lontitude  =  self.ref_longtitude   
        # stand_lontitude: A real value specifying, for ARW, the longitude that is parallel with the y-axis 
        # in the Lambert conformal and polar stereographic projections.    
        self.NCOLS = NCOLS # the number of columns (dimensionality in the X direction).
        self.NROWS = NROWS # the number of rows (dimensionality in the Y direction)
        self.X_CELL = X_CELL # the cell dimension parallel to the X coordinate axis, meters
        self.Y_CELL = Y_CELL # the cell dimension parallel to the Y coordinate axis, meters 
        self.X_ORIG = X_ORIG # X coordinate of the grid origin (lower left corner of the cell at column=row=1, meters)
        self.Y_ORIG = Y_ORIG # Y coordinate of the grid origin (lower left corner of the cell at column=row=1, meters)

        # Define projection information based on parameters
        self.LCCProj = pyproj.Proj(proj='lcc', lat_1=self.truelat_first, lat_2=self.truelat_second,\
             lat_0=self.ref_latitude, lon_0=self.ref_longtitude, a=6370000, b=6370000)
        self.LCCProj_crs = cartopy.crs.LambertConformal(central_longitude=self.ref_longtitude, \
            central_latitude=self.ref_latitude, standard_parallels =(self.truelat_first, self.truelat_second))
        self.LCC_plotextent = [self.X_ORIG, self.X_ORIG + self.X_CELL * self.NCOLS, self.Y_ORIG,
                    self.Y_ORIG + self.Y_CELL * self.NROWS]
        
        # Construct a few important coordinate arrays
        # Transform the target LCC to WGS84 and Web Mercator EPSG:3857 
        self.LccXList = numpy.tile(numpy.arange(self.X_ORIG + self.X_CELL / 2, self.X_ORIG \
            + self.NCOLS * self.X_CELL,self.X_CELL), self.NROWS) # Inner lOOP X FIRST, THEN OUTER LOOP Y
        self.LccYList = numpy.repeat(numpy.arange(self.Y_ORIG + self.Y_CELL/2, self.Y_ORIG \
            + self.NROWS * self.Y_CELL,self.Y_CELL), self.NCOLS, axis=None) # Inner lOOP X FIRST, THEN OUTER LOOP Y 
        self.LccPXList = numpy.tile(numpy.arange(self.X_ORIG , self.X_ORIG \
            + (self.NCOLS + 1) * self.X_CELL,self.X_CELL), self.NROWS + 1)
        self.LccPYList = numpy.repeat(numpy.arange(self.Y_ORIG , self.Y_ORIG \
            + (self.NROWS + 1) * self.Y_CELL,self.Y_CELL), self.NCOLS + 1, axis=None)
        
        self.dstLonList, self.dstLatList = self.LCCProj(self.LccXList, self.LccYList, inverse=True)
        self.PdstLonList, self.PdstLatList = self.LCCProj(self.LccPXList, self.LccPYList, inverse=True)
        self.Grid_LonArray = numpy.zeros((self.NROWS,self.NCOLS))
        self.Grid_LatArray = numpy.zeros((self.NROWS,self.NCOLS))
        self.Points_LonArray = numpy.zeros((self.NROWS + 1,self.NCOLS + 1))
        self.Points_LatArray = numpy.zeros((self.NROWS + 1,self.NCOLS + 1))
        self.Grid_LonArray = self.dstLonList.reshape(self.NROWS,self.NCOLS)
        self.Grid_LatArray = self.dstLatList.reshape(self.NROWS,self.NCOLS)
        self.Points_LonArray = self.PdstLonList.reshape(self.NROWS + 1,self.NCOLS + 1)
        self.Points_LatArray = self.PdstLatList.reshape(self.NROWS + 1,self.NCOLS + 1)
        combine_list =  [list(t) for t in zip(self.LccXList,self.LccYList)]
        combine_list = [tuple(combine_list[n]) for n in  range(len(combine_list))]        
        self.web_Mercator_List = pyproj.itransform(pyproj.crs.CRS.from_proj4(self.LCCProj.to_proj4()), "epsg:3857", combine_list)
        self.web_Mercator_List = [list(pt) for pt in self.web_Mercator_List]