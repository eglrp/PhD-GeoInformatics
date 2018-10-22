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
demFile = r"V:/Data/NGI/GEF DEM/3323d_2015_1001_GEF_DEM_Photoscan_clip.tif"
filtDemFile = r"V:\Data\NGI\GEF DEM\3323d_2015_1001_GEF_DEM_Photoscan_clip_filt.tif"
plantHeightFile = r"V:\Data\NGI\GEF DEM\3323d_2015_1001_GEF_DEM_Photoscan_clip_hgt.tif"

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
# demDs = None

pylab.close('all')

se = morphology.disk(5)
# filtDem = morphology.opening(dem, se)           # opencv will be faster
filtDem = cv2.morphologyEx(dem, cv2.MORPH_ERODE, se)
se = morphology.disk(5)
filtDem = cv2.morphologyEx(filtDem, cv2.MORPH_ERODE, se)


# filtDem = cv2.medianBlur(dem, 5)  # 5 is the biggest size for use with float32
# filtDem = cv2.bilateralFilter(dem, 13, 150, 150)

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


###########################################################################################
#  write filtered DEM & plant height files

filtDemDs = gdal.GetDriverByName('GTiff').CreateCopy(filtDemFile, demDs, options=["TILED=YES", "COMPRESS=DEFLATE"])
filtDemDs.GetRasterBand(1).WriteArray(filtDem)   # Writes my array to the raster
filtDemDs.FlushCache()
filtDemDs = None

plantHeightDs = gdal.GetDriverByName('GTiff').CreateCopy(plantHeightFile, demDs, options=["TILED=YES", "COMPRESS=DEFLATE"])
plantHeightDs.GetRasterBand(1).WriteArray(dem-filtDem)   # Writes my array to the raster
plantHeightDs.FlushCache()
plantHeightDs = None

(dem - filtDem < 0).sum()

##############################################################################################


blurDem = cv2.medianBlur(dem, 5)  # 5 is the biggest size for use with float32
for i in range(0,20):
    blurDem = cv2.medianBlur(blurDem, 5)  # 5 is the biggest size for use with float32
# blurDem = cv2.medianBlur(blurDem, 5)  # 5 is the biggest size for use with float32
# blurDem = cv2.medianBlur(blurDem, 5)  # 5 is the biggest size for use with float32
# blurDem = cv2.medianBlur(blurDem, 5)  # 5 is the biggest size for use with float32
# blurDem = cv2.medianBlur(blurDem, 5)  # 5 is the biggest size for use with float32

ls = LightSource(azdeg=315, altdeg=45)
pylab.figure()
ax1 = pylab.subplot(311)
pylab.imshow(ls.hillshade(dem, vert_exag=1., dx=.5, dy=.5), cmap='gray')
pylab.subplot(312, sharex=ax1, sharey=ax1)
pylab.imshow(ls.hillshade(blurDem, vert_exag=1., dx=.5, dy=.5), cmap='gray')
pylab.subplot(313, sharex=ax1, sharey=ax1)
pylab.imshow(ls.hillshade(dem-blurDem, vert_exag=1., dx=.5, dy=.5), cmap='gray')


#############################################################################################
# experiment with srtm

srtmFile = r"D:\Data\Development\Projects\PhD GeoInformatics\Data\CGA\SRTM\s34_e023-e024_1arc_v3_Mosaic_NGI_ProjGrid.tif"
plantHeightFile = r"V:\Data\NGI\GEF DEM\3323d_2015_1001_GEF_DEM_Photoscan_clip_hgt2.tif"

srtmDs = gdal.OpenEx(srtmFile, gdal.OF_RASTER)
if srtmDs is None:
    print "Open failed./n"

print 'Driver: ', srtmDs.GetDriver().ShortName, '/', \
    srtmDs.GetDriver().LongName
print 'Size is ', srtmDs.RasterXSize, 'x', srtmDs.RasterYSize, \
    'x', srtmDs.RasterCount
print 'Projection is ', srtmDs.GetProjection()
geotransform = srtmDs.GetGeoTransform()
if not geotransform is None:
    print 'Origin = (', geotransform[0], ',', geotransform[3], ')'
    print 'Pixel Size = (', geotransform[1], ',', geotransform[5], ')'
    pixelSize = geotransform[1]

srtm = srtmDs.GetRasterBand(1).ReadAsArray()

pylab.figure()
ax1 = pylab.subplot(311)
pylab.imshow(dem)
pylab.subplot(312, sharex=ax1, sharey=ax1)
pylab.imshow(srtm)
pylab.subplot(313, sharex=ax1, sharey=ax1)
pylab.imshow(dem-srtm)
# pylab.colorbar()


se = morphology.disk(15)
# filtDem = morphology.opening(dem, se)           # opencv will be faster
filtDem2 = cv2.morphologyEx(dem-srtm, cv2.MORPH_OPEN, se)
filtDem2 += srtm

ls = LightSource(azdeg=315, altdeg=45)
pylab.figure()
ax1 = pylab.subplot(311)
pylab.imshow(ls.hillshade(dem, vert_exag=1., dx=.5, dy=.5), cmap='gray')
pylab.subplot(312, sharex=ax1, sharey=ax1)
pylab.imshow(ls.hillshade(filtDem2, vert_exag=1., dx=.5, dy=.5), cmap='gray')
pylab.subplot(313, sharex=ax1, sharey=ax1)
pylab.imshow(ls.hillshade(dem-filtDem2, vert_exag=1., dx=.5, dy=.5), cmap='gray')


plantHeightDs = gdal.GetDriverByName('GTiff').CreateCopy(plantHeightFile, demDs, options=["TILED=YES", "COMPRESS=DEFLATE"])
plantHeightDs.GetRasterBand(1).WriteArray(dem-filtDem2)   # Writes my array to the raster
plantHeightDs.FlushCache()
plantHeightDs = None

(dem - filtDem < 0).sum()


