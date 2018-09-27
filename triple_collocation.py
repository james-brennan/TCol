"""
triple_collocation.py

for a given modis tile estimate the TriCol

-- Can use either log model or additive model
"""
import os
import glob
import gdal
from utils import *
import numpy as np
import sys
import datetime
#from sklearn.covariance import MinCovDet
import pytesmo.scaling as scaling
import pytesmo.metrics as metrics

global NDAYS
NDAYS = 7


def triple_collocation(mcd64, mcd45, fcci, dates):
    """
    first try at Tcol
    """
    #nObs = mcd64.shape[0]
    st = np.vstack((mcd64, mcd45, fcci)).T
    """
    Want to isolate to the fire season
    for now...
    Use percentile to mask out where fires are rare
    """

    """
    remove months with no fire for all products...
    """
    all0 = np.all(st==0, axis=1)
    perc10 = np.percentile(st[~all0], 5, axis=0)
    #ff = np.all(st>perc10, axis=1)
    mask = ~all0 #& ~ff
    st = st[mask]
    nObs = st.shape[0]
    x = st[:, 0]
    y = st[:, 1]
    z = st[:, 2]
    # re-scale to x
    # see if re-scaling to y
    # changes results...
    y_sca = scaling.mean_std(y, x)
    z_sca = scaling.mean_std(z, x)
    """
    if all nan something went wrong
    """
    y_sca[~np.isfinite(y_sca)]=0
    z_sca[~np.isfinite(z_sca)]=0
    # and if negative? -- not sure it matters?
    #y_sca[y_sca<0]=0
    #z_sca[z_sca<0]=0
    # perform triple collocation
    e_x, e_y, e_z = metrics.tcol_error(x, y_sca, z_sca)
    errors = np.array([e_x, e_y, e_z])
    """
    return what we need
    """
    return (errors,
            nObs)







def do_TriCol(tile, N=240, start_year=2005, end_year=2012, fireCCI5=True):
    """

    make versatile enough to vary experiments...

    Arguments:
        tile [string] -- modis tile to use
        N [int] -- effective grid size.. go for 150 for now?
                N = 240 corresponds to around 1degree lat lon...
    """
    nS = int(2400/N)
    sig_MCD64 =   -999*np.ones((nS, nS))
    sig_MCD45 =   -999*np.ones((nS, nS))
    sig_FCCI  =   -999*np.ones((nS, nS))
    nObs      =    -999*np.ones((nS, nS))

    # which cci product?
    if fireCCI5:
        loadFireCCI = loadFireCCIv50
    else:
        loadFireCCI = loadFireCCI41

    for i, x0 in enumerate(xrange(0, 2400, N)):
        for j, y0 in enumerate(xrange(0, 2400, N)):
            x1 = x0 + N
            y1 = y0 + N
            x1 = np.minimum(x1, 2400)
            y1 = np.minimum(y1, 2400)
            #print "Doing %i %i " %(y0, x0)
            xmin, xmax, ymin, ymax = x0, x1, y0, y1
            """
            get data for each year
            """
            dates = []
            mcd64 = []
            mcd45 = []
            fcci = []
            for year in xrange(start_year, end_year+1):
                """
                run function to get each product and BA field
                for the whole tile
                """
                mcd45_, ds = loadMCD45(year, tile, ymin, ymax, xmin, xmax)
                mcd64_ = loadMCD64(year, tile, ymin, ymax, xmin, xmax)
                fcci_ = loadFireCCI(year, tile, ymin, ymax, xmin, xmax)
                """
                Aggregate the burnt area..
                """
                mcd45H, _ = np.histogram(mcd45_, range(1, 366, NDAYS))
                mcd64H, _ = np.histogram(mcd64_, range(1, 366, NDAYS))
                fcciH, _ = np.histogram(fcci_, range(1, 366, NDAYS))
                #
                mcd45.append(mcd45H)
                mcd64.append(mcd64H)
                fcci.append(fcciH)
                # save dates
                dates.append([datetime.date(year, 1, 1) + datetime.timedelta(n) for n in range(1,366, NDAYS)] )
            dates = np.hstack(dates)
            mcd64 = np.hstack(mcd64)
            mcd45 = np.hstack(mcd45)
            fcci = np.hstack(fcci)
            try:
                """
                Do the actual calculation
                """
                (   _sigmas,
                    _nObs) = triple_collocation(mcd64, mcd45, fcci, dates)
                #print _sigmas
                #import pdb; pdb.set_trace()
                _sigMCD64 = _sigmas[0]
                _sigMCD45 = _sigmas[1]
                _sigFCCI =  _sigmas[2]
                #if np.sum(mcd64>0)>1:
                #    import pdb; pdb.set_trace()
                # keep it
                sig_MCD64[j, i]     =  _sigMCD64
                sig_MCD45[j, i]     =  _sigMCD45
                sig_FCCI[j, i]      =  _sigFCCI
                nObs[j, i]          =  _nObs
            except:
                pass
    out = OutputsDataset(  sig_MCD64=sig_MCD64,
                    sig_FCCI=sig_FCCI,
                    sig_MCD45=sig_MCD45,
                    start_year=start_year,
                    end_year=end_year,
                    tile=tile,
                    FireCCI50=fireCCI5,
                    nObs=nObs)


    # make an output dataset for this
    return out
