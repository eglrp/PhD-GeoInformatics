import pylab
import os
import subprocess
import glob
import gdal
import ogr
# import rasterio
import numpy as np
from rasterio.transform import Affine

# NB setup env to run Anaconda GDAL
os.environ['PATH'] += "C:\ProgramData\Anaconda3\envs\py27\Library\\bin"
os.environ['GDAL_DATA'] = "C:\ProgramData\Anaconda3\envs\py27\Library\share\gdal"
os.environ['GDAL_DRIVER_PATH']="C:\ProgramData\Anaconda3\envs\py27\Library\\bin\gdalplugins"

gdaltindexExe = "C:\ProgramData\Anaconda3\envs\py27\Library\\bin\gdaltindex.exe"  # Use anaconda version as osgeo4w wont run under Anaconda
rootDir = "D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\XCalib Experiments\Calibrated2"
doRaw = False
ovwLevel = 0
if doRaw:
    tileIndexWildCard = '*_CMP.tif'
    rootDir = "D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\XCalib Experiments"
else:
    tileIndexWildCard = '*_CMP_XCALIB.tif'
    rootDir = "D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\XCalib Experiments\Calibrated2"

# read a ROI of a raster with windowCnrs in projected co-ords, optionally read the same ROI but from the overviews
# NB NOTE - overviews are not tapped/grid aligned and so will not spatially align properly between adjacent images
# we should recode overviews as our own downsampled ims - should be quite easy with numpy (if not slow)
def ReadRasterWindow(rasterDs, roiCnrs = None, ovwLevel = None):
    geoTransform = np.array(rasterDs.GetGeoTransform())
    if ovwLevel is not None:  # adjust the geotransform to be for larger pixels
        geoTransform[1] = geoTransform[1] * (2. ** (ovwLevel+1))
        geoTransform[5] = geoTransform[5] * (2. ** (ovwLevel+1))
    affine = Affine.from_gdal(*geoTransform)

    if roiCnrs is None:
        origin = np.int32([0, 0])
        counts = np.int32([rasterDs.XSize, rasterDs.YSize])
    else:
        roiPlCnrs = np.array(~affine * roiCnrs.transpose()).transpose()
        origin = np.min(roiPlCnrs, axis=0)
        counts = np.max(roiPlCnrs, axis=0) - origin
        # print 'origin: ' + str(origin)
        # print 'counts: ' + str(counts)
        origin = np.int32(np.round(origin))  # only round after counts calc
        counts = np.int32(np.round(counts))
        if ovwLevel is not None:
            band = rasterDs.GetRasterBand(1).GetOverview(ovwLevel)
        else:
            band = rasterDs.GetRasterBand(1)
        rasterSize = np.array([band.XSize, band.YSize])
        # Note A - it is possible that the ROI extends past the image boundary due to rounding - fix this below
        mask = (origin + counts) > rasterSize
        counts[mask] = (rasterSize - origin)[mask]
        # print 'origin: ' + str(origin)
        # print 'counts: ' + str(counts)

    rasterRoi = np.zeros([counts[1], counts[0], rasterDs.RasterCount], dtype=np.float32)
    for b in range(0, rasterDs.RasterCount):
        if ovwLevel is not None:
            band = rasterDs.GetRasterBand(b + 1).GetOverview(ovwLevel)
        else:
            band = rasterDs.GetRasterBand(b + 1)
        bandArray = band.ReadAsArray(origin[0], origin[1], counts[0], counts[1], None, None, gdal.GDT_Float32)
        rasterRoi[:, :, b] = bandArray

    # use last band to get mask - assume mask same for all bands
    rasterMask = band.GetMaskBand().ReadAsArray(origin[0], origin[1], counts[0], counts[1])
    return rasterRoi, rasterMask  #, origin, counts

if doRaw:
    subDirs = ['Source2']
else:
    subDirs = os.listdir(rootDir)
    results = {}

for dirItem in subDirs:
    if not os.path.isdir(os.path.join(rootDir, dirItem)):
        continue
    sourceDir = os.path.join(rootDir, dirItem)
    print 'Processing {0:%s}' % (sourceDir)
    print 'Generating tile index'
    tileIndexFileName = os.path.join(sourceDir, os.path.basename(sourceDir) + '_TileIndex.shp')

    if False:
        if os.path.exists(tileIndexFileName):
            for file in glob.glob(os.path.splitext(tileIndexFileName)[0] + '.*'):
                print 'Removing ' + file
                os.remove(file)

    if not os.path.exists(tileIndexFileName):
        for file in glob.glob(os.path.join(sourceDir, tileIndexWildCard)):
            subprocess.call('{0} "{1}" "{2}"'.format(gdaltindexExe, tileIndexFileName, file), shell=True, env=os.environ)
            print 'Adding ' + file

    # open tile index field and read into list of dicts
    tileIndexDs = gdal.OpenEx(tileIndexFileName, gdal.OF_VECTOR)
    if tileIndexDs is None:
        raise "Open %s failed."%(tileIndexFileName)

    lyr = tileIndexDs.GetLayerByIndex(0)
    lyr.ResetReading()

    tileList = []
    for (i, feat) in enumerate(lyr):
        print '.',
        feat_defn = lyr.GetLayerDefn()
        f = {}
        for i in range(feat_defn.GetFieldCount()):
            field_defn = feat_defn.GetFieldDefn(i)
            f[field_defn.GetName()] = feat.GetField(i)

        geom = feat.GetGeometryRef()
        if geom is not None and (geom.GetGeometryType() == ogr.wkbPolygon):
            f['geom'] = geom.Clone()
            f['points'] = geom.GetBoundary().GetPoints()
        else:
            f['geom'] = geom.Clone()
            print "no polygon geometry"
            # break
        tileList.append(f)
        # gcpList.append(f)

    tileIndexDs = None


    # nested loop to find and measure all overlaps
    absDiffAccum = np.array([0.,0.,0.,0.])
    numPixelAccum = 0.
    reflScale = np.int64(5000)
    maxAbsDiffAccum = np.array([0.,0.,0.,0.])
    maxAccum = np.array([0.,0.,0.,0.])
    normDiffAccum = []
    breakFlag = False
    for (outerI, outerTile) in enumerate(tileList):
        # print 'outerI: ' + str(outerI)
        for (innerI, innerTile) in enumerate(tileList[outerI + 1:]):
            # print 'innerI: ' + str(innerI)
            if outerTile['geom'].Intersects(innerTile['geom']):
                print os.path.basename(outerTile['location']) + ' intersects ' + os.path.basename(innerTile['location'])
                # if '3321D_319_02_0060_RGBN_CMP_XCALIB.tif' == os.path.basename(outerTile['location']) and '3321D_319_02_0062_RGBN_CMP_XCALIB.tif' == os.path.basename(innerTile['location']):
                #     break

                # break
                intersectionPoly = outerTile['geom'].Intersection(innerTile['geom'])
                print 'Intersection area (ha): %3f'%(intersectionPoly.GetArea()/1.e4)
                intersectionCnrs = np.array(intersectionPoly.GetBoundary().GetPoints())[0:-1]
                outerDs = gdal.OpenEx(outerTile['location'], gdal.OF_RASTER)
                innerDs = gdal.OpenEx(innerTile['location'], gdal.OF_RASTER)
                try:
                    outerRaster, outerMask = ReadRasterWindow(outerDs, intersectionCnrs, ovwLevel=ovwLevel)
                    innerRaster, innerMask = ReadRasterWindow(innerDs, intersectionCnrs, ovwLevel=ovwLevel)
                finally:
                    outerDs = None
                    innerDs = None
                # due to Note A, it is possible the above may be different sizes - fix this
                sizeFix = np.array([outerMask.shape, innerMask.shape]).min(axis=0)
                outerRaster = outerRaster[:sizeFix[0], :sizeFix[1]]
                outerMask = outerMask[:sizeFix[0], :sizeFix[1]]
                innerRaster = innerRaster[:sizeFix[0], :sizeFix[1]]
                innerMask = innerMask[:sizeFix[0], :sizeFix[1]]

                intersectionMask = outerMask & innerMask
                # NOTE - some sort of r2 may be good to measure overlapping areas  - it will easier allow comparison to raw uncalib im
                # ya but... we actually want to measure in absolute terms as the overlap areas should be matched in absolute terms
                diffRaster = np.abs(outerRaster - innerRaster)
                diffRaster[np.logical_not(intersectionMask)] = 0.  # broadcast across dims
                if False:
                    pylab.figure()
                    ax = pylab.subplot(2, 2, 1)
                    pylab.imshow(outerRaster[:, :, [3, 0, 1]] / np.max(outerRaster[:, :, [3, 0, 1]], axis=(0, 1)))
                    pylab.subplot(2, 2, 2, sharex=ax, sharey=ax)
                    pylab.imshow(innerRaster[:, :, [3, 0, 1]] / np.max(innerRaster[:, :, [3, 0, 1]], axis=(0, 1)))
                    pylab.subplot(2, 2, 3, sharex=ax, sharey=ax)
                    pylab.imshow(diffRaster[:, :, [3, 0, 1]] / np.max(diffRaster[:, :, [3, 0, 1]], axis=(0, 1)))
                    pylab.subplot(2, 2, 4, sharex=ax, sharey=ax)
                    pylab.imshow(intersectionMask)

                nPixel = np.int64((intersectionMask > 0)).sum()
                sumAbsDiff = np.sum(np.int64(np.abs(diffRaster[intersectionMask > 0])), axis=0)  #mean for each band
                if nPixel > 0:
                    maxAbsDiff = np.max(np.abs(diffRaster[intersectionMask > 0]), axis=0)  #max for each band
                    maxAbsDiffAccum[maxAbsDiff > maxAbsDiffAccum] = maxAbsDiff[maxAbsDiff > maxAbsDiffAccum]
                    mx = outerRaster[intersectionMask > 0].max(axis=0)
                    maxAccum[mx > maxAccum] = mx[mx > maxAccum]
                    normDiffAccum.append((np.float64(sumAbsDiff)/(nPixel*np.std(outerRaster[intersectionMask > 0], axis=0))).tolist())

                if nPixel == 0:
                    print 'No valid intersection area'
                else:
                    print 'Valid intersection area (ha): %.3f'%(nPixel*(.5**2)/1.e4)
                print 'Mean Abs Error: {0}'.format(100. * np.float64(sumAbsDiff) / (nPixel*reflScale))
                numPixelAccum += nPixel
                absDiffAccum += sumAbsDiff
        #         if np.any(100. * np.float32(sumAbsDiff) / (nPixel*reflScale)<0):
        #             breakFlag = True
        #             break
        #     if breakFlag:
        #         break
        # if breakFlag:
        #     break
    result = {}
    result['name'] = os.path.basename(sourceDir)
    if doRaw:
        result['winSize'] = -1
        result['winArea'] = -1
        result['model'] = -1
    else:
        result['winSize'] = (np.int(result['name'][1]), np.int(result['name'][2]))
        result['winArea'] = np.int(result['name'][1]) * np.int(result['name'][2])
        result['model'] = np.int(result['name'][4])

    result['mae'] = 100.*np.float64(absDiffAccum)/(reflScale*numPixelAccum)
    result['mean(mae)'] = np.mean(result['mae'])
    result['maxAe'] = 100.*np.float64(maxAbsDiffAccum)/reflScale
    result['maxVal'] = maxAccum
    result['area'] = 100. * np.float64(numPixelAccum) * (.5 ** 2) / 1.e4
    result['homog.'] = 100. * np.float64(absDiffAccum)/(maxAccum*numPixelAccum)
    result['mean(homog.)'] = np.mean(result['homog.'])
    result['normMAE'] = np.array(normDiffAccum).mean(axis=0)
    result['mean(normMAE)'] = result['normMAE'].mean()

    print 'Mean Abs Error: {0}'.format(result['mae'])
    print 'Mean Abs Error: {0}'.format(result['mean(mae)'])
    print 'Max Abs Error: {0}'.format(result['maxAe'])
    print 'Homogeneity: {0}'.format(result['homog.'])
    print 'Norm Abs Error: {0}'.format(result['normMAE'])
    print 'Mean Homogeneity: {0}'.format(result['mean(homog.)'])
    print 'Ttl intersected area (ha): {0}'.format(result['area'])
    results[result['name']] = result

#
#
# rasterDs = outerDs
# roiCnrs = intersectionCnrs
# ovwLevel = 2
# b = 0

#####################################################################################################
# visualise

models = np.array([r['model'] for r in results.values()])
winSizes = np.array([r['winSize'] for r in results.values()])
winAreas = np.array([r['winArea'] for r in results.values()])
maes = np.array([r['mean(mae)'] for r in results.values()])
homog = np.array([r['mean(homog.)'] for r in results.values()])
normAE = np.array([r['mean(normMAE)'] for r in results.values()])

fontSize = 12.
# mpl.rcParams.update({'font.size': fontSize})

pylab.figure()
for model in np.unique(models):
    modelIdx = models == model
    modelWinAreas = winAreas[modelIdx]
    sortIdx = np.argsort(modelWinAreas)
    pylab.subplot(3,1,1)
    pylab.plot(modelWinAreas[sortIdx], (maes[modelIdx])[sortIdx], 'x-')
    pylab.subplot(3,1,2)
    pylab.plot(modelWinAreas[sortIdx], (homog[modelIdx])[sortIdx], 'x-')
    pylab.subplot(3,1,3)
    pylab.plot(modelWinAreas[sortIdx], (normAE[modelIdx])[sortIdx], 'x-')

pylab.subplot(3, 1, 1)
pylab.xlabel('Win. area (pixels)')
pylab.ylabel('Mean Abs. error (%)')
pylab.grid()
pylab.legend(np.unique(models))
pylab.subplot(3, 1, 2)
pylab.xlabel('Win. area (pixels)')
pylab.ylabel('Homog. (%)')
pylab.grid()
pylab.legend(np.unique(models))
pylab.subplot(3, 1, 3)
pylab.xlabel('Win. area (pixels)')
pylab.ylabel('Norm AE (%)')
pylab.grid()
pylab.legend(np.unique(models))
