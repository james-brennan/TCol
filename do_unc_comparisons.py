"""
do_unc_comparisions.py
"""
import gdal
import os
import sys
import numpy as np
import glob
import skimage.measure
"""
*-- GFED uncertainties --*
    as we are in the MODIS era these
    are scalings by GFED region...
"""
def make_GFED_yearly_uncs(tile, start_year=2002, end_year=2014):
    """
    This produces a GFED4 unc tiff
    with the uncertainties at the resolution
    we use...

    """
    ODIR = '/group_workspaces/cems2/nceo_generic/users/jbrennan01/TCol2/outputs/FULL_RECORD/'
    GFED_DIR = "/home/users/jbrennan01/DATA2/TCol2/GFED4/"
    my_tile_ref = gdal.Open(ODIR+"output_%s_120.tif" % tile)

    # limit to the tile
    geo_t = my_tile_ref.GetGeoTransform()
    year_uncs = []

    for year in xrange(start_year, end_year+1):
        files = glob.glob(GFED_DIR+"*%i*hdf" % year)
        files.sort()
        # file template
        tmpl = 'HDF4_SDS:UNKNOWN:"%s":1'
        # save each month
        month_unc = []
        for m in files:
            a = gdal.Open(tmpl%m)
            """
            Reproject to the modis grid...
            and then limit to the tile extent...

            we need to ascribe
            the projection and geo_transform of the GFED data...

            we have this from the regions file...
            """
            modis_warp = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '
            regions = gdal.Open("/group_workspaces/cems2/nceo_generic/users/jbrennan01/TCol2/GFED_regions.tif")
            gfed_warp = regions.GetProjection()
            a.SetGeoTransform(regions.GetGeoTransform())
            a.SetProjection(gfed_warp)
            # do the warp...
            tile_gfed = gdal.Warp('', a, format='MEM', resampleAlg=gdal.GRA_Average,
                                     dstSRS=modis_warp,
                                        outputBounds=[geo_t[0], geo_t[3]+geo_t[-1]*my_tile_ref.RasterYSize,
                                        geo_t[0]+geo_t[1]*my_tile_ref.RasterXSize,
                                        geo_t[3] ])
            """
            load it and apply the scaling factor
            """
            file_scale = 0.009999999776
            arr_full_res = tile_gfed.ReadAsArray() * file_scale
            """
            downsampling from 0.25 to ~0.5
            -- eg sum up for the aggregation
            """
            arr = skimage.measure.block_reduce(arr_full_res, (2, 2), np.sum)
            """
            convert from hectares to km2
            """
            to_km2 = 0.01
            arr = arr * to_km2
            #import ipdb; ipdb.set_trace()
            # convert to variance...
            unc = arr**2
            month_unc.append(unc)
        month_unc = np.array(month_unc)
        # save
        #import pdb; pdb.set_trace()
        year_uncs.append(month_unc.sum(axis=0))
    year_uncs = np.array(year_uncs)
    # return
    return year_uncs, my_tile_ref

def save_GFED_uncs(tile, uncs, tile_ds):
    # save these
    outdir = '/group_workspaces/cems2/nceo_generic/users/jbrennan01/TCol2/GFEDUncs/'
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
    nBands = uncs.shape[0]
    int_nSize = uncs.shape[1]
    int_nPxls =  int(2400/int_nSize)
    dest = mem_drv.Create(outdir+'gfed4_yearlyUncertainties_%i_%s.tif' % (int_nPxls, tile), int_nSize,
                                 int_nSize, nBands, gdal.GDT_Float32)
    # Calculate the geotransform
    dest.SetGeoTransform( tile_ds.GetGeoTransform() )
    #import pdb; pdb.set_trace()
    dest.SetProjection ( tile_ds.GetProjectionRef() )
    # and write
    for band in xrange(1, nBands+1):
        dest.GetRasterBand(band).WriteArray(uncs[band-1])
    # set nodata
    for b in xrange(1, nBands+1):
        band = dest.GetRasterBand(b)
        band.SetNoDataValue(-999)
    # Write to disk.
    dest.FlushCache()
    dest = None


"""
*-- FireCCI50 uncertainties --*
"""
def make_FIRECCI50_yearly_uncs(tile, start_year=2002, end_year=2014):
    datadir = '/datacentre/archvol2/pan80/archive/spot-2243-esacci_fire/data/burned_area/MODIS/grid/v5.0/'

    ODIR = '/group_workspaces/cems2/nceo_generic/users/jbrennan01/TCol2/outputs/FULL_RECORD/'
    GFED_DIR = "/home/users/jbrennan01/DATA2/TCol2/GFED4/"
    my_tile_ref = gdal.Open(ODIR+"output_%s_120.tif" % tile)

    # limit to the tile
    geo_t = my_tile_ref.GetGeoTransform()
    year_uncs = []

    for year in xrange(start_year, end_year+1):
        files_ = glob.glob(datadir + "/%i/"%year + "*nc")
        files_.sort()
        """
        each file covers a 2-week period...
        """
        week_unc = []
        for fi in files_:
            tmpl = 'NETCDF:"%s":standard_error'
            a = gdal.Open(tmpl % fi)

            """
            reproject... and cut to
            the extent of the modis tile
            """
            regions = gdal.Open("/group_workspaces/cems2/nceo_generic/users/jbrennan01/TCol2/GFED_regions.tif")
            gfed_warp = regions.GetProjection()
            a.SetGeoTransform(regions.GetGeoTransform())
            a.SetProjection(gfed_warp)

            modis_warp = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '
            # do the warp...
            tile_fcci = gdal.Warp('', a, format='MEM', resampleAlg=gdal.GRA_Average,
                                     dstSRS=modis_warp,
                                        outputBounds=[geo_t[0], geo_t[3]+geo_t[-1]*my_tile_ref.RasterYSize,
                                        geo_t[0]+geo_t[1]*my_tile_ref.RasterXSize,
                                        geo_t[3] ])

            """
            load it and apply the scaling factor
            """
            arr_full_res = tile_fcci.ReadAsArray()
            """
            downsampling from 0.25 to ~0.5
            -- eg sum up for the aggregation
            """
            arr = skimage.measure.block_reduce(arr_full_res, (2, 2), np.sum)
            """
            convert from m2 to km2
            """
            to_km2 = 1e-6
            arr = arr * to_km2
            #import ipdb; ipdb.set_trace()
            # convert standard errors to variances...
            unc = arr**2
            week_unc.append(unc)
        week_unc  = np.array(week_unc)
        # save to yearly
        year_uncs.append(week_unc.sum(axis=0))
    year_uncs = np.array(year_uncs)
    # return
    return year_uncs, my_tile_ref

def save_FireCCI50_uncs(tile, uncs, tile_ds):
    # save these
    outdir = '/group_workspaces/cems2/nceo_generic/users/jbrennan01/TCol2/FCCI50Uncs/'
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
    nBands = uncs.shape[0]
    int_nSize = uncs.shape[1]
    int_nPxls =  int(2400/int_nSize)
    dest = mem_drv.Create(outdir+'fcci50_yearlyUncertainties_%i_%s.tif' % (int_nPxls, tile), int_nSize,
                                 int_nSize, nBands, gdal.GDT_Float32)
    # Calculate the geotransform
    dest.SetGeoTransform( tile_ds.GetGeoTransform() )
    #import pdb; pdb.set_trace()
    dest.SetProjection ( tile_ds.GetProjectionRef() )
    # and write
    for band in xrange(1, nBands+1):
        dest.GetRasterBand(band).WriteArray(uncs[band-1])
    # set nodata
    for b in xrange(1, nBands+1):
        band = dest.GetRasterBand(b)
        band.SetNoDataValue(-999)
    # Write to disk.
    dest.FlushCache()
    dest = None




"""
*-- Annual TC uncertainties --*
"""
def do_annual_TC_uncs(product='MCD64', tile='h30v10'):
    """
    This produces annual uncertainties
    for the TC method
    """
    start_year = 2002
    end_year = 2014
    dataDIR = '/home/users/jbrennan01/DATA2/TCol2/CMG_products/%s/' % product
    # get the CMG-esque file
    file_ = dataDIR + "%s_120_%s.tif" % (product, tile)
    g  = gdal.Open(file_)
    arr = g.ReadAsArray()
    """
    Year-splits
    are every 52 weeks...

    we want to store the number of
    weeks per year with burning
    """
    annual_weeks = []
    NWEEKS = 52
    for k in xrange(start_year-2002, end_year+1-2002):
        i0 = NWEEKS * k
        i1 = i0 + NWEEKS
        wk = (arr[i0:i1]>0)
        nwk = wk.sum(axis=0)
        annual_weeks.append(nwk)
    annual_weeks =np.array(annual_weeks)
    """
    Now we need the product sigmas...
    """
    uncDIR = '/group_workspaces/cems2/nceo_generic/users/jbrennan01/TCol2/outputs/FULL_RECORD/'
    file_2 = uncDIR + "output_%s_120.tif" % (tile)
    product_dict = {'MCD64':0, 'MCD45':1, 'FireCCI50':2, 'FireCCI41':2}
    g = gdal.Open(file_2)
    uncs = g.GetRasterBand(product_dict[product] +1).ReadAsArray()
    annual_uncs = uncs * annual_weeks
    return annual_uncs, g


def save_TC_annual_uncs(product, tile, uncs,  tile_ds):
    outdir = '//group_workspaces/cems2/nceo_generic/users/jbrennan01/TCol2/derived/%s/' % product
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
    nBands = uncs.shape[0]
    int_nSize = uncs.shape[1]
    int_nPxls =  int(2400/int_nSize)
    dest = mem_drv.Create(outdir+'%s_annualUncertainties_%i_%s.tif' % (product, int_nPxls, tile), int_nSize,
                                 int_nSize, nBands, gdal.GDT_Float32)
    # Calculate the geotransform
    dest.SetGeoTransform( tile_ds.GetGeoTransform() )
    #import pdb; pdb.set_trace()
    dest.SetProjection ( tile_ds.GetProjectionRef() )
    # and write
    for band in xrange(1, nBands+1):
        dest.GetRasterBand(band).WriteArray(uncs[band-1])
    # set nodata
    for b in xrange(1, nBands+1):
        band = dest.GetRasterBand(b)
        band.SetNoDataValue(-999)
    # Write to disk.
    dest.FlushCache()
    dest = None





if __name__ == "__main__":


    tiles = np.genfromtxt("world_tiles.txt", dtype=str)
    for tile in tiles:
        try:
            """
            GFED annual uncertainties
            """
            u, tile_ds = make_GFED_yearly_uncs(tile=tile)
            save_GFED_uncs(tile, u, tile_ds)
            """
            FireCCI50 product annual uncertainties
            """
            u, tile_ds = make_FIRECCI50_yearly_uncs(tile=tile)
            save_FireCCI50_uncs(tile, u, tile_ds)
            """
            MCD64 annual uncertainties
            """
            u, tile_ds = do_annual_TC_uncs("MCD64", tile)
            save_TC_annual_uncs("MCD64", tile, u, tile_ds)
            """
            MCD45 annual uncertainties
            """
            u, tile_ds = do_annual_TC_uncs("MCD45", tile)
            save_TC_annual_uncs("MCD45", tile, u, tile_ds)
            """
            FireCCI50 TC annual uncertainties
            """
            u, tile_ds = do_annual_TC_uncs("FireCCI50", tile)
            save_TC_annual_uncs("FireCCI50", tile, u, tile_ds)
            print tile
        except:
            print("failed: ", tile)
