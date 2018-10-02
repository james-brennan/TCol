"""
plots.py
"""


import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kde

sns.set_style("dark")
# create data
sns.set_context("paper")


fig, ax = plt.subplots(figsize=(10, 4.5), nrows=1, ncols=2, sharex=True, sharey=False)

for i in xrange(2):
    x = d[i].flatten()
    y = c[i].flatten()
    # make a mask for fill values
    mask = np.logical_or(x<=0, y<=0)
    x = x[~mask]
    y = y[~mask]
    # plot the scatter
    ax[i].plot(x,y, 'k.',markersize=1, alpha=1, markeredgewidth=0.0, rasterized=True)
    ax[i].set_xlim(0, 1e3)
    ax[i].set_ylim(0, 1e3)
    ax[i].plot([0, 1e3], [0, 1e3], 'grey')
    # robust line...
    ransac = linear_model.RANSACRegressor()
    ransac.fit(x.reshape((-1, 1)), y)
    inlier_mask = ransac.inlier_mask_
    outlier_mask = np.logical_not(inlier_mask)
    line_X = np.arange(x.min(), x.max())[:, np.newaxis]
    line_y = ransac.predict(line_X)
    line_y_ransac = ransac.predict(line_X)
    ax[i].plot(line_X, line_y_ransac, 'darkred')
    # add info to the plot
    tmpl = "r$^2$=%.2f y=%.2f x + %.2f" % (ransac.score(x.reshape((-1, 1)), y), ransac.estimator_.coef_, ransac.estimator_.intercept_)
    ax[i].text(0.05, 0.9,tmpl,size=8,
         horizontalalignment='left',
         verticalalignment='center',
         transform = ax[i].transAxes)


# add labels
ax[0].set_xlabel(r"MCD45 $\sigma^2$ [km$^4$] | FireCCI41")
ax[0].set_ylabel(r"MCD45 $\sigma^2$ [km$^4$] | FireCCI50")
# add labels
ax[1].set_xlabel(r"MCD64 $\sigma^2$ [km$^4$] | FireCCI41")
ax[1].set_ylabel(r"MCD64 $\sigma^2$ [km$^4$] | FireCCI50")
# add labels
ax[2].set_xlabel(r"FireCCI41 $\sigma^2$ [km$^4$]")
ax[2].set_ylabel(r"FireCCI50 $\sigma^2$ [km$^4$]")

plt.tight_layout()

plt.savefig("stability_plot2.pdf", dpi=300, bbox_inches='tight')


# Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
nbins=300
k = kde.gaussian_kde([x,y])
xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))

# Make the plot
plt.pcolormesh(xi, yi, z.reshape(xi.shape))
plt.show()

# Change color palette
plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.Greens_r)
plt.show()

import statsmodels.api as sm


dens_u = sm.nonparametric.KDEMultivariate(data=[x,y],  var_type='cc', bw='normal_reference')

z = dens_u.pdf(np.vstack([xi.flatten(), yi.flatten()]))



"""
Now number of obs plot
"""


from mpl_toolkits.basemap import Basemap
import osr, gdal
import matplotlib.pyplot as plt
import numpy as np



def convertXY(xy_source, inproj, outproj):
    # function to convert coordinates

    shape = xy_source[0,:,:].shape
    size = xy_source[0,:,:].size

    # the ct object takes and returns pairs of x,y, not 2d grids
    # so the the grid needs to be reshaped (flattened) and back.
    ct = osr.CoordinateTransformation(inproj, outproj)
    xy_target = np.array(ct.TransformPoints(xy_source.reshape(2, size).T))

    xx = xy_target[:,0].reshape(shape)
    yy = xy_target[:,1].reshape(shape)

    return xx, yy

# Read the data and metadata
ds = gdal.Open(r'world.vrt')
# reproject
wgs84_warp = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs '
ds = gdal.Warp('', ds, format='MEM', dstSRS=wgs84_warp, resampleAlg=gdal.GRA_Average)


data = ds.ReadAsArray()
gt = ds.GetGeoTransform()
proj = ds.GetProjection()
xres = gt[1]
yres = gt[5]
# get the edge coordinates and add half the resolution
# to go to center coordinates
xmin = gt[0] + xres * 0.5
xmax = gt[0] + (xres * ds.RasterXSize) - xres * 0.5
ymin = gt[3] + (yres * ds.RasterYSize) + yres * 0.5
ymax = gt[3] - yres * 0.5
xy_source = np.mgrid[xmin:xmax+xres:xres, ymax+yres:ymin:yres]

data[data==-999]=np.nan
data[data<5]=np.nan


# Create the figure and basemap object
fig = plt.figure(figsize=(12, 6))
m = Basemap(projection='robin', lon_0=0, resolution='c')

# Create the projection objects for the convertion
# original (Albers)
inproj = osr.SpatialReference()
inproj.ImportFromWkt(proj)

# Get the target projection from the basemap object
outproj = osr.SpatialReference()
outproj.ImportFromProj4(m.proj4string)

# Convert from source projection to basemap projection
xx, yy = convertXY(xy_source, inproj, outproj)

# plot the data (first layer)
# norm=mpl.colors.LogNorm(),
im1 = m.pcolormesh(xx, yy, data[3,:,:].T, cmap=plt.cm.cubehelix_r, vmin=0.1, vmax=450, rasterized=True)

# annotate
#m.drawcountries()
m.drawcoastlines(linewidth=.99)

cbar = plt.colorbar(im1, shrink=0.5)
cbar.set_label("Number of valid collocates \n [2002-2015]")

plt.tight_layout()
plt.savefig("nObs.pdf", dpi=300, bbox_inches='tight')