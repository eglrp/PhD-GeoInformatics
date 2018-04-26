import gdal
import numpy as np
import glob
import os
import pylab


# This finds AOD stats for the little karoo study area for NGI survey days
# it is crude as it does not take into account the sub areas surveyed on each day

sourceDir = 'D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MOD04'

def world2Pixel(geoMatrix, x, y):
    """
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate
    """
    ulX = geoMatrix[0]
    ulY = geoMatrix[3]
    xDist = geoMatrix[1]
    yDist = geoMatrix[5]
    rtnX = geoMatrix[2]
    rtnY = geoMatrix[4]
    pixel = ((x - ulX) / xDist)
    line = ((y - ulY) / yDist)
    return (pixel, line)


# The HDF files are in some weird unspecified projection, QGIS somehow knows how to figure it out.
# So I created these geotiffs with QGIS from the HDF's
resList = []
for file in glob.glob(os.path.join(sourceDir, 'MOD04_L2*.tif')):
    print file
    # break

    ds = gdal.OpenEx(file, gdal.OF_RASTER)
    gt = ds.GetGeoTransform()
    im = ds.ReadAsArray()
    origin = world2Pixel(gt, 21.5, -33.)
    brCnr = world2Pixel(gt, 22.5, -34.)
    sides = np.round(brCnr) - np.round(origin)
    origin = np.int32(origin)
    sides = np.int32(sides)
    brCnr = np.int32(np.round(brCnr))
    imRoi = im[origin[1]:origin[1]+sides[1], origin[0]:origin[0]+sides[0]]
    res = {}
    res['imRoi'] = imRoi
    res['file'] = file
    res['max'] = imRoi.max()
    validMask = imRoi>-1.
    res['min'] = imRoi[validMask].min()
    res['mean'] = imRoi[validMask].mean()
    res['median'] = np.median(imRoi[validMask])
    resList.append(res)

    ds = None

mx = [res['max'] for res in resList]
mn = [res['min'] for res in resList]
mean = [res['mean'] for res in resList]
median = [res['median'] for res in resList]
shapes = [res['imRoi'].shape for res in resList]

pylab.figure()
pylab.plot(mx, 'x-', label='max')
pylab.plot(mn, 'x-', label='min')
pylab.plot(mean, 'x-', label='mean')
pylab.plot(median, 'x-', label='median')
pylab.legend()

print 'Mean AOD over surveys: %.3f' % (np.mean(mean))
print 'Median AOD over surveys: %.3f' % (np.median(mean))


# subDsList = ds.GetSubDatasets()
# idx = [i for (i, subDsItem) in enumerate(subDsList) if subDsItem[0].endswith('AOD_550_Dark_Target_Deep_Blue_Combined')]
# subDs = gdal.OpenEx(subDsList[idx[0]][0], gdal.OF_RASTER)
# im = subDs.ReadAsArray()
# nodataMask = im < -1   # ?
