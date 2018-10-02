"""
make_CMG_products.py

This makes compataible CMG products
for the BA products...
"""
import os
import glob
import gdal
from utils import *
import numpy as np
import sys
import datetime
#from sklearn.covariance import MinCovDet

global NDAYS
NDAYS = 7

# constant to convert 500m data to km2
global AREA_CONST
AREA_CONST = 0.25

def make_CMG_prod(loadFunction, tile, N=120, start_year=2005, end_year=2012):
    """
    I think loading the products into
    memory won't take too much space?
    """
    # get tile geostuff
    data_, refds = loadMCD45(2008, tile, 0, 20, 0, 20)


    # create storage
    NYEARS = end_year - start_year + 1
    nT = 52 * NYEARS
    nS = int(2400/N)
    grid = -999*np.ones((nT, nS, nS))

    # run over areas
    for i, x0 in enumerate(xrange(0, 2400, N)):
        for j, y0 in enumerate(xrange(0, 2400, N)):
            x1 = x0 + N
            y1 = y0 + N
            x1 = np.minimum(x1, 2400)
            y1 = np.minimum(y1, 2400)
            #print "Doing %i %i " %(y0, x0)
            xmin, xmax, ymin, ymax = x0, x1, y0, y1
            """
            run function to get each product and BA field
            for the whole tile
            """
            for year in xrange(start_year, end_year+1):
                if loadFunction == loadMCD45:
                    data_, ds = loadFunction(year, tile, ymin, ymax, xmin, xmax)
                else:
                    data_ = loadFunction(year, tile, ymin, ymax, xmin, xmax)
                """
                Aggregate the burnt area..
                """
                dataH, _ = np.histogram(data_, range(1, 366, NDAYS))
                #
                t0 = (year-start_year)*52
                t1 = t0 + 52
                grid[t0:t1, j, i] = dataH * AREA_CONST
                 # save dates
                #dates.append([datetime.date(year, 1, 1) + datetime.timedelta(n) for n in range(1,366, NDAYS)] )
    return grid, refds


def saveCMG_product(product, tile, grid, tile_ds):
    outdir = "/home/users/jbrennan01/DATA2/TCol2/CMG_products/" + product + "/"
    # stuff stores the outputs
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    """
    set up the output file
    """
    mem_drv = gdal.GetDriverByName( 'GTIff' )
    # The size of the raster is given the new projection and pixel spacing
    # Using the values we calculated above. Also, setting it to store one band
    # and to use Float32 data type.
    nBands = grid.shape[0]
    int_nSize = grid.shape[1]
    int_nPxls =  int(2400/int_nSize)
    dest = mem_drv.Create(outdir+'%s_%i_%s.tif' % (product, int_nPxls, tile), int_nSize,
                                 int_nSize, nBands, gdal.GDT_Float32)
    # Calculate the new geotransform
    # Set the geotransform
    geo_t = list(tile_ds.GetGeoTransform())
    # do pixel re-sizing
    geo_t[1] =  (geo_t[1] * 2400.0) / int_nSize
    geo_t[-1] =  (geo_t[-1] * 2400.0) / int_nSize
    geo_t = tuple(geo_t)
    dest.SetGeoTransform( geo_t )
    dest.SetProjection ( tile_ds.GetProjectionRef() )
    # and write
    for band in xrange(1, nBands+1):
        dest.GetRasterBand(band).WriteArray(grid[band-1])
    # set nodata
    for b in xrange(1, nBands+1):
        band = dest.GetRasterBand(b)
        band.SetNoDataValue(-999)
    # Write to disk.
    dest.FlushCache()
    dest = None

if __name__ == "__main__":

    tile = sys.argv[1]
    """
    mcd64 2002-2015
    """
    grid, refds = make_CMG_prod(loadMCD45, tile=tile, N=120, start_year=2002, end_year=2014)
    saveCMG_product("MCD45", tile, grid, refds)
    grid, refds = make_CMG_prod(loadMCD64, tile=tile, N=120, start_year=2002, end_year=2014)
    saveCMG_product("MCD64", tile, grid, refds)

    grid, refds = make_CMG_prod(loadFireCCIv50, tile=tile, N=120, start_year=2002, end_year=2014)
    saveCMG_product("FireCCI50", tile, grid, refds)

    grid, refds = make_CMG_prod(loadFireCCI41, tile=tile, N=120, start_year=2005, end_year=2011)
    saveCMG_product("FireCCI41", tile, grid, refds)
