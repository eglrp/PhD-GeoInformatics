import pylab
import os
import subprocess
import glob

# NB setup env to run Anaconda GDAL
os.environ['PATH'] += "C:\ProgramData\Anaconda3\envs\py27\Library\\bin"
os.environ['GDAL_DATA'] = "C:\ProgramData\Anaconda3\envs\py27\Library\share\gdal"
os.environ['GDAL_DRIVER_PATH']="C:\ProgramData\Anaconda3\envs\py27\Library\\bin\gdalplugins"

sourceDir = "D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\XCalib Experiments\Calibrated2\w11p1"
gdaltindexExe = "C:\ProgramData\Anaconda3\envs\py27\Library\\bin\gdaltindex.exe"  # Use anaconda version as osgeo4w wont run under Anaconda

print 'Processing {0:%s}' % (sourceDir)
tileIndexFileName = os.path.join(sourceDir, os.path.basename(sourceDir) + '_TileIndex.shp')
# tileIndexFileName = 'TileIndex.shp'
if False:
    if os.path.exists(tileIndexFileName):
        for file in glob.glob(os.path.splitext(tileIndexFileName)[0] + '.*'):
            print 'Removing ' + file
            os.remove(file)

    for file in glob.glob(os.path.join(sourceDir, '*_CMP_XCALIB.tif')):
        subprocess.call('{0} "{1}" "{2}"'.format(gdaltindexExe, tileIndexFileName, file), shell=True, env=os.environ)
        print 'Adding ' + file

###############################################
# import fiona, shapely
#
# tileIndexFile = fiona.open(tileIndexFileName)
# print tileIndexFile.meta
# for f in tileIndexFile:
#     print f
import gdal
import ogr
# import rasterio
import numpy as np
from rasterio.transform import Affine

tileIndexDs = gdal.OpenEx(tileIndexFileName, gdal.OF_VECTOR)
if tileIndexDs is None:
    raise "Open %s failed."%(tileIndexFileName)

lyr = tileIndexDs.GetLayerByIndex(0)
lyr.ResetReading()

# if lyr.GetSpatialRef() is not None and lyr.GetSpatialRef() is not 0:
#     dgpsSrs = lyr.GetSpatialRef()
# else:
#     dgpsSrs.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

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


# read a ROI of a raster with windowCnrs in projected co-ords, optionally read the same ROI but from the overviews
# mmm - this overview thing is not going to work well - the overviews of different ims are not tapped / aligned to the same grid
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
        print 'origin: ' + str(origin)
        print 'counts: ' + str(counts)
        origin = np.int32(np.round(origin))
        counts = np.int32(np.round(counts))

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


# read a ROI of a raster with windowCnrs in projected co-ords, optionally read the same ROI but from the overviews
# mmm - this overview thing is not going to work well - the overviews of different ims are not tapped / aligned to the same grid
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
        print 'origin: ' + str(origin)
        print 'counts: ' + str(counts)
        origin = np.int32(np.round(origin))
        counts = np.int32(np.round(counts))

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


# read a ROI of a raster with windowCnrs in projected co-ords, optionally read the same ROI but from the overviews
# mmm - this overview thing is not going to work well - the overviews of different ims are not tapped / aligned to the same grid
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
        print 'origin: ' + str(origin)
        print 'counts: ' + str(counts)
        origin = np.int32(np.round(origin))  # only round after counts calc
        counts = np.int32(np.round(counts))

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


for (outerI, outerTile) in enumerate(tileList):
    # print 'outerI: ' + str(outerI)
    for (innerI, innerTile) in enumerate(tileList[outerI + 1:]):
        # print 'innerI: ' + str(innerI)
        if outerTile['geom'].Intersects(innerTile['geom']):
            print os.path.basename(outerTile['location']) + ' intersects ' + os.path.basename(innerTile['location'])
            intersectionPoly = outerTile['geom'].Intersection(innerTile['geom'])
            intersectionCnrs = np.array(intersectionPoly.GetBoundary().GetPoints())[0:-1]
            outerDs = gdal.OpenEx(outerTile['location'], gdal.OF_RASTER)
            innerDs = gdal.OpenEx(innerTile['location'], gdal.OF_RASTER)
            outerRaster, outerMask = ReadRasterWindow(outerDs, intersectionCnrs, ovwLevel=0)
            innerRaster, innerMask = ReadRasterWindow(innerDs, intersectionCnrs, ovwLevel=0)
            intersectionMask = outerMask & innerMask
            diffRaster = np.abs(outerRaster - innerRaster)
            diffRaster[np.logical_not(intersectionMask)] = 0.  # broadcast across dims

            # outerAffine = Affine.from_gdal(*outerRaster.GetGeoTransform())
            # innerAffine = Affine.from_gdal(*innerRaster.GetGeoTransform())
            # outerCnrs = np.array(~outerAffine * intersectionCnrs.transpose()).transpose()
            # innerCnrs = np.array(~innerAffine * intersectionCnrs.transpose()).transpose()
            # outerUl = np.min(outerCnrs, axis=0)
            # outerSize = np.max(outerCnrs, axis=0) - outerUl
            # innerUl = np.min(innerCnrs, axis=0)
            # innerSize = np.max(innerCnrs, axis=0) - innerUl


            # tmp = outerRaster.read(None, window=((0,0),(10,10)))
            break
    break


rasterDs = innerDs
roiCnrs = intersectionCnrs
ovwLevel = 0
b = 0

pylab.figure()
ax = pylab.subplot(2,2,1)
pylab.imshow(outerRaster[:,:,[3,0,1]]/np.max(outerRaster[:,:,[3,0,1]], axis=(0,1)))
pylab.subplot(2,2,2, sharex=ax, sharey=ax)
pylab.imshow(innerRaster[:,:,[3,0,1]]/np.max(innerRaster[:,:,[3,0,1]], axis=(0,1)))
pylab.subplot(2,2,3, sharex=ax, sharey=ax)
pylab.imshow(diffRaster[:,:,[3,0,1]]/np.max(diffRaster[:,:,[3,0,1]], axis=(0,1)))
pylab.subplot(2,2,4, sharex=ax, sharey=ax)
pylab.imshow(intersectionMask)

