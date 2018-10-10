import gdal
import ogr
import numpy as np
import osr
import pylab
from scipy import stats as stats
from matplotlib import patches
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from scipy import ndimage as ndimage
from scipy.stats import gaussian_kde

import skimage.morphology as morphology

# Python Imaging Library imports
from PIL import Image
from PIL import ImageDraw

import cv2

demFile = r"V:\Data\NGI\GEF DEM\3323d_2015_1001_GEF_DEM_SGM3_clip.tif"

demDs = gdal.OpenEx(demFile, gdal.OF_RASTER)
if demDs is None:
    print "Open failed./n"

print 'Driver: ', demDs.GetDriver().ShortName, '/', \
    demDs.GetDriver().LongName
print 'Size is ', demDs.RasterXSize, 'x', demDs.RasterYSize, \
    'x', demDs.RasterCount
print 'Projection is ', demDs.GetProjection()
geotransform = demDs.GetGeoTransform()
if not geotransform is None:
    print 'Origin = (', geotransform[0], ',', geotransform[3], ')'
    print 'Pixel Size = (', geotransform[1], ',', geotransform[5], ')'
    pixelSize = geotransform[1]

dem = demDs.GetRasterBand(1).ReadAsArray()
demDs = None

pylab.close('all')

se = morphology.disk(15)
# filtDem = morphology.opening(dem, se)           # opencv will be faster
filtDem = cv2.morphologyEx(dem, cv2.MORPH_OPEN, se)

se = morphology.disk(10)
# filtDem = morphology.opening(dem, se)           # opencv will be faster
filtDem = cv2.morphologyEx(filtDem, cv2.MORPH_OPEN, se)

se = morphology.disk(5)
# filtDem = morphology.opening(dem, se)           # opencv will be faster
filtDem = cv2.morphologyEx(filtDem, cv2.MORPH_OPEN, se)

pylab.figure()
ax1 = pylab.subplot(211)
pylab.imshow(dem)
pylab.subplot(212, sharex=ax1, sharey=ax1)
pylab.imshow(filtDem)


ls = LightSource(azdeg=315, altdeg=45)
pylab.figure()
ax1 = pylab.subplot(311)
pylab.imshow(ls.hillshade(dem, vert_exag=1., dx=.5, dy=.5), cmap='gray')
pylab.subplot(312, sharex=ax1, sharey=ax1)
pylab.imshow(ls.hillshade(filtDem, vert_exag=1., dx=.5, dy=.5), cmap='gray')
pylab.subplot(313, sharex=ax1, sharey=ax1)
pylab.imshow(ls.hillshade(dem-filtDem, vert_exag=1., dx=.5, dy=.5), cmap='gray')

# filtDem15to5 = filtDem
# filtDem5to15 = filtDem
# filtDem15 = filtDem

ls = LightSource(azdeg=315, altdeg=45)
pylab.figure()
ax1 = pylab.subplot(221)
pylab.imshow(ls.hillshade(dem, vert_exag=1., dx=.5, dy=.5), cmap='gray')
pylab.subplot(222, sharex=ax1, sharey=ax1)
pylab.imshow(ls.hillshade(filtDem15, vert_exag=1., dx=.5, dy=.5), cmap='gray')
pylab.subplot(223, sharex=ax1, sharey=ax1)
pylab.imshow(ls.hillshade(filtDem15to5, vert_exag=1., dx=.5, dy=.5), cmap='gray')
pylab.subplot(224, sharex=ax1, sharey=ax1)
pylab.imshow(ls.hillshade(filtDem5to15, vert_exag=1., dx=.5, dy=.5), cmap='gray')
