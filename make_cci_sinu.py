"""
Produce a MODIS like fire_cci product...
"""

import os
import glob
import gdal
import numpy as np


def make_vrts():
    """
    this seems to take awhile too
    """
    for year in xrange(2001, 2015):
        merisdir = "/neodc/esacci/fire/data/burned_area/MODIS/pixel/v5.0/compressed/%i/" % year
        storeDir = "/home/users/jbrennan01/DATA2/TColBA/input_products/vrts/"
        """
        to make things simpler we make a vrt
        to put the different meris areas together
        """
        store = []
        for month in xrange(1, 13):
            """
            want to produce a modis style intermediate product...
            """
            files = glob.glob(merisdir+"%i%02i*-fv5.0.tar.gz" % (year, month))
            vsitar_strings = []
            for f in files:
                _dir = f.strip(f.split("/")[-1]) # lazy...
                filename = f.split("/")[-1]
                # construct the true tif filename
                tif_ = filename.strip(".tar.gz") + "-JD.tif"
                """
                Use GDAL vsitar function!

                -- Idea is to make a global vrt of each area...
                """ 
                strr = "/vsitar//%s/%s" % (f, tif_)
                vsitar_strings.append(strr)
            """
            Make a VRT of each area to produce a global file
            """
            #import pdb; pdb.set_trace()
            vsitar_files = ' '.join(vsitar_strings)
            vrt_file = 'fcci_v5_%02i_%i.vrt' % (month, year)
            os.system("gdalbuildvrt %s/%s %s &> /dev/null" % (storeDir, vrt_file, vsitar_files) )
        print month, year




def make_modis_tiles(year, month, tile):
    """
    From the VRT files extract and produce a modis cut tile of these...
    """
    mcd45dir ="/group_workspaces/cems2/nceo_generic/users/jlgomezdans/%s/%i/MCD45A1/" % (tile, year)
    store = []
    files = glob.glob(mcd45dir+"*.hdf")
    files.sort()
    tmp = 'HDF4_EOS:EOS_GRID:"%s":MOD_GRID_Monthly_500km_BA:burndate' % files[0]
    refModis = gdal.Open(tmp)
    # sort output info
    out_dir = "/home/users/jbrennan01/DATA2/TColBA/input_products/fire_cci_v5/%s/" % tile
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    out_file = "fcci_v5_%s_%02i_%i.tif" % (tile, month, year)

    # FIND the relevant global VRT
    mdir = "/home/users/jbrennan01/DATA2/TColBA/input_products/vrts/"
    file_vrt = glob.glob(mdir + "fcci_v5_%02i_%i*vrt" % (month, year))[0]
    # load it
    g = gdal.Open(file_vrt )
    geo_t = g.GetGeoTransform ()
    """
    Now we want to limit the whole file
    just down to the extent of the MODIS tile
    """
    # Now, we create an in-memory raster
    mem_drv = gdal.GetDriverByName( 'GTiff' )
    dest = mem_drv.Create(out_dir + out_file, 2400, 2400, 1, gdal.GDT_Int16,
                            options = [ 'COMPRESS=LZW' ])
    # Calculate the new geotransform
    # Set the geotransform
    dest.SetGeoTransform( refModis.GetGeoTransform() )
    dest.SetProjection ( refModis.GetProjectionRef() )
    #import pdb; pdb.set_trace()
    # Perform the projection/resampling
    # this i think will take awhile?
    res = gdal.ReprojectImage( g, dest, \
        g.GetProjectionRef(), refModis.GetProjectionRef(), \
        gdal.GRA_NearestNeighbour )
    # flush
    dest = None


if __name__ == "__main__":

    """
    make the vrt files...
    """
    #make_vrts()



    #get modis tiles
    tiles = np.genfromtxt("world_tiles.txt", dtype=str)
    for modis in tiles:
        for year in xrange(2001, 2015):
            for month in xrange(1, 13):
                try:
                    make_modis_tiles(year, month, modis)
                except:
                    pass
                print year, month, modis


