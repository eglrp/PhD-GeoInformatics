# Fill SRTM DEM nodata areas with ASTER DEM

import numpy as np
import pylab
import os
import scipy.stats
import gdal
import scipy.ndimage
import matplotlib.image
from matplotlib import pyplot

if True:
    srtmFileName = "D:/Data/Development/Projects/PhD GeoInformatics/Data/CGA/SRTM/s34_e021_1arc_v3.tif"
    aterFileName = "D:/Data/Development/Projects/PhD GeoInformatics/Data/CGA/ASTER DEM/ASTGTM2_S34E021_dem.tif"
    fillSrtmFileName = "D:/Data/Development/Projects/PhD GeoInformatics/Data/CGA/SRTM/s34_e021_1arc_v3_AsterFill.tif"
else:
    srtmFileName = "D:/Data/Development/Projects/PhD GeoInformatics/Data/CGA/SRTM/s34_e022_1arc_v3.tif"
    aterFileName = "D:/Data/Development/Projects/PhD GeoInformatics/Data/CGA/ASTER DEM/ASTGTM2_S34E022_dem.tif"


srtmDs = gdal.Open(srtmFileName, gdal.GA_ReadOnly)
asterDs = gdal.Open(aterFileName, gdal.GA_ReadOnly)

srtmIm = srtmDs.GetRasterBand(1).ReadAsArray()
asterIm = asterDs.GetRasterBand(1).ReadAsArray()

#these rasters are identical sizes and projections

nodata = -32767
mask = (srtmIm == nodata)
fillIm = srtmIm.copy()
fillIm[mask] = asterIm[mask]

srtmIm[mask] = 0   # for display

pylab.figure()
ax1 = pylab.subplot(131)
pylab.imshow(srtmIm)
pylab.title('SRTM')
pylab.subplot(132, sharex=ax1, sharey=ax1)
pylab.imshow(asterIm)
pylab.title('ASTER')
pylab.subplot(133, sharex=ax1, sharey=ax1)
pylab.imshow(np.abs(fillIm-asterIm))
pylab.title('|ASTER-SRTM|')
pylab.colorbar()

pylab.figure()
ax1 = pylab.subplot(131)
pylab.imshow(srtmIm)
pylab.title('SRTM')
pylab.subplot(132, sharex=ax1, sharey=ax1)
pylab.imshow(asterIm)
pylab.title('ASTER')
pylab.subplot(133, sharex=ax1, sharey=ax1)
pylab.imshow(fillIm)
pylab.title('SRTM[nodata] = ASTER')


# write out filled raster
driver = gdal.GetDriverByName('GTiff')
outDs = driver.CreateCopy(fillSrtmFileName, srtmDs, strict=0)
outDs.GetRasterBand(1).WriteArray(fillIm)
outDs.GetRasterBand(1).FlushCache()
outDs.FlushCache()
outDs = None

# mosaic dem's
# gdalwarp s34_e021_1arc_v3_AsterFill.tif s34_e022_1arc_v3.tif s34_e021_e022_1arc_v3_Mosaic.tif

Re = 6371.
h = 832.
incidence = -11.873132
view = np.arcsin(np.sin(incidence*np.pi/180.)*Re/(Re+h))*180/np.pi