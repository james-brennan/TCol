import os
import glob
import gdal
import numpy as np
from collections import namedtuple
import json


"""
*-- IO stuff --*
"""
def loadFireCCIv50(year, tile, ymin, ymax, xmin, xmax):
    """
    *-- Load FireCCI50 --*
    """
    ysize = ymax - ymin
    xsize = xmax - xmin
    fccidir = "/home/users/jbrennan01/DATA2/TColBA/input_products/fire_cci_v5/%s/" % tile
    store = []
    for month in xrange(1, 13):
        #import pdb; pdb.set_trace()
        filename = "fcci_v5_%s_%02i_%i.tif" % (tile, month, year)
        dest = gdal.Open(fccidir+filename)
        """
        Load the data and find burnt pixels
        """
        mapB = dest.GetRasterBand(1).ReadAsArray(xoff=xmin, yoff=ymin, win_xsize=xsize, win_ysize=ysize)
        burnt_m = mapB[np.logical_and(mapB>0, mapB < 366)]
        store.append(burnt_m)
    store = np.hstack(store)
    return store

def loadFireCCI41(year, tile, ymin, ymax, xmin, xmax):
    """
    *-- Load FireCCI41 --*
    """
    ysize = ymax - ymin
    xsize = xmax - xmin
    merisdir = "/group_workspaces/cems2/nceo_generic/users/jbrennan01/TColBA/input_products/fireCCI_v4.1/%i/" % year
    """
    get the modis info
    """
    mcd45dir ="/group_workspaces/cems2/nceo_generic/users/jlgomezdans/%s/%i/MCD45A1/" % (tile, year)
    store = []
    files = glob.glob(mcd45dir+"*.hdf")
    files.sort()
    tmp = 'HDF4_EOS:EOS_GRID:"%s":MOD_GRID_Monthly_500km_BA:burndate' % files[0]
    refModis = gdal.Open(tmp)
    """
    to make things simpler we make a vrt
    to put the different meris areas together
    """
    store = []
    for month in xrange(1, 13):
        """
        save it because we re-use it later...
        """
        tmpdir = '/home/users/jbrennan01/DATA2/TColBA/tmp/'
        tmpfile = tmpdir + 'fcci_tmp_%02i_%i_%s.tif' % (month, year, tile)
        print tile, os.path.isfile(tmpfile)
        if not os.path.isfile(tmpfile):
            """
            want to produce a modis style intermediate product...
            """
            files = glob.glob(merisdir+"%i%02i*tif" % (year, month))
            # make a vrt to combine these things...
            f = ' '.join(files)
            vrt_file = 'fcci_tmp_%02i_%i_%s.vrt' % (month, year, tile)
            os.system("gdalbuildvrt %s/%s %s &> /dev/null" % (merisdir, vrt_file, f) )
            # load it
            g = gdal.Open("%s/%s" % (merisdir, vrt_file) )
            geo_t = g.GetGeoTransform ()
            """
            Now we want to limit the whole file
            just down to the extent of the MODIS tile
            """
            # Now, we create an in-memory raster
            mem_drv = gdal.GetDriverByName( 'GTiff' )
            dest = mem_drv.Create(tmpfile, 2400, 2400, 3, gdal.GDT_Int16,
                                    options = [ 'COMPRESS=LZW' ])
            # Calculate the new geotransform
            # Set the geotransform
            dest.SetGeoTransform( refModis.GetGeoTransform() )
            dest.SetProjection ( refModis.GetProjectionRef() )
            #import pdb; pdb.set_trace()
            # Perform the projection/resampling
            res = gdal.ReprojectImage( g, dest, \
                g.GetProjectionRef(), refModis.GetProjectionRef(), \
                gdal.GRA_NearestNeighbour, )
            print 'creating tmpfile'
        else:
            dest = gdal.Open(tmpfile)
        """
        Load the data and find burnt pixels
        """
        mapB = dest.GetRasterBand(1).ReadAsArray(xoff=xmin, yoff=ymin, win_xsize=xsize, win_ysize=ysize)
        burnt_m = mapB[np.logical_and(mapB>0, mapB < 366)]
        store.append(burnt_m)
    store = np.hstack(store)
    return store


def loadMCD45(year, tile, ymin, ymax, xmin, xmax):
    """
    *-- Load MCd45 --*
    """
    ysize = ymax - ymin
    xsize = xmax - xmin
    mcd45dir ="/group_workspaces/cems2/nceo_generic/users/jlgomezdans/%s/%i/MCD45A1/" % (tile, year)
    store = []
    F = []
    files = glob.glob(mcd45dir+"*.hdf")
    for f in files:
        tmp = 'HDF4_EOS:EOS_GRID:"%s":MOD_GRID_Monthly_500km_BA:burndate' % f
        mapB = gdal.Open(tmp).ReadAsArray(xoff=xmin, yoff=ymin, xsize=xsize, ysize=ysize)
        burnt_m = mapB[np.logical_and(mapB>0, mapB < 366)]
        store.append(burnt_m)
        F.append(mapB)
    store = np.hstack(store)
    return store, gdal.Open(tmp)

def loadMCD64(year, tile,  ymin, ymax, xmin, xmax):
    """
    *-- Load MCD64 --*
    """
    ysize = ymax - ymin
    xsize = xmax - xmin
    #mcd64dir ="/group_workspaces/cems2/nceo_generic/users/jlgomezdans/%s/%i/MCD64A1/" % (tile, year)
    # new data!
    mcd64dir = "/group_workspaces/cems2/nceo_generic/users/jbrennan01/TColBA/input_products/MCD64/%s/%i/" % (tile, year)
    store = []
    files = glob.glob(mcd64dir+"*.hdf")
    for f in files:
        tmp = 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_Monthly_500m_DB_BA:Burn Date' % f
        mapB = gdal.Open(tmp).ReadAsArray(xoff=xmin, yoff=ymin, xsize=xsize, ysize=ysize)
        burnt_m = mapB[np.logical_and(mapB>0, mapB < 366)]
        store.append(burnt_m)
    store = np.hstack(store)
    return store



# the stuff we'd like to save...
OutputDataset = namedtuple("OutputDataset", ['sig_MCD64', 'sig_FCCI', 'sig_MCD45', 'start_year', 'end_year', 'tile', 'FireCCI50', 'nObs']  )

def save_experiment(exp_name, output_tuple, tile, tile_ds):
    outdir = "/home/users/jbrennan01/DATA2/TCol2/outputs/" + exp_name + "/"
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

    int_nSize = output_tuple.sig_MCD64.shape[0]
    int_nPxls =  int(2400/int_nSize)
    dest = mem_drv.Create(outdir+'output_%s_%i.tif' % (tile, int_nPxls), int_nSize,
                                 int_nSize, 4, gdal.GDT_Float32)
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
    dest.GetRasterBand(1).WriteArray(output_tuple.sig_MCD64)
    dest.GetRasterBand(2).WriteArray(output_tuple.sig_MCD45)
    dest.GetRasterBand(3).WriteArray(output_tuple.sig_FCCI)
    dest.GetRasterBand(4).WriteArray(output_tuple.nObs)

    # set nodata
    for b in xrange(1, 4):
        band = dest.GetRasterBand(b)
        band.SetNoDataValue(-999)
    # Write to disk.
    dest.FlushCache()
    dest = None
    """
    make a file to save experiment stuff
    """
    options_dict = { 'start_year': output_tuple.start_year,
                     'end_year':   output_tuple.end_year,
                     'FireCCI50': output_tuple.FireCCI50,
                     'Nx': int_nPxls,}
    # save these
    with open(outdir+'%s.txt' % exp_name, 'w') as file:
        file.write(json.dumps(options_dict))



"""
*-- Other utilities --*
"""
def clean_up_tmpfiles(year, tile):
    """
    run this file at the end to remove the
    tempoary rasters for the firecci product
    """
    for month in xrange(1, 13):
        """
        save it because we re-use it later...
        """
        tmpdir = '/home/users/jbrennan01/DATA2/TColBA/tmp/'
        tmpfile = tmpdir + 'fcci_tmp_%02i_%i_%s.tif' % (month, year, tile)
        if os.path.isfile(tmpfile):
            # delete the file
            pass
        else:
            pass
