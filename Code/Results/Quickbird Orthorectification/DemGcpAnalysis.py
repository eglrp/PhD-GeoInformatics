# Analyse differences between GCP heights and DEM - is there something systematic?
import gdal
import ogr
import osr
import numpy as np
import pylab
from matplotlib import patches
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

demFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/CGA/SUDEM/x3324c_15_15_L2a_crop.tif"
gcpGeoLocFile = "C:/Data/Development/Projects/PhD GeoInformatics/Docs/Misc/Baviaanskloof/BaviiaansPeCorrectedGcpMay2017Combined.shp"


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


ds = gdal.OpenEx(demFile, gdal.OF_RASTER)
if ds is None:
    print "Open failed./n"

print 'Driver: ', ds.GetDriver().ShortName,'/', \
      ds.GetDriver().LongName
print 'Size is ',ds.RasterXSize,'x',ds.RasterYSize, \
      'x',ds.RasterCount
print 'Projection is ',ds.GetProjection()
demGeoTransform = ds.GetGeoTransform()
if not demGeoTransform is None:
    print 'Origin = (',demGeoTransform[0], ',',demGeoTransform[3],')'
    print 'Pixel Size = (',demGeoTransform[1], ',',demGeoTransform[5],')'

demSpatialRef = osr.SpatialReference(ds.GetProjection())

dem = ds.ReadAsArray()
ds = None

pylab.figure()
pylab.imshow(dem)

# read in Geo CoOrds of GCP's (and height!)
ds = gdal.OpenEx(gcpGeoLocFile, gdal.OF_VECTOR)
if ds is None:
    print "Open failed./n"

lyr = ds.GetLayerByIndex(0)
lyr.ResetReading()
gcpSpatialRef = lyr.GetSpatialRef()

gcpToDemTransform = osr.CreateCoordinateTransformation(gcpSpatialRef, demSpatialRef)

# gcpList = []
gcpDict = {}
for (i, feat) in enumerate(lyr):
    if i > 190:
        break
    print '.',
    feat_defn = lyr.GetLayerDefn()
    f = {}
    for i in range(feat_defn.GetFieldCount()):
        field_defn = feat_defn.GetFieldDefn(i)
        f[field_defn.GetName()] = feat.GetField(i)
    geom = feat.GetGeometryRef()
    if geom is not None and (geom.GetGeometryType() == ogr.wkbPoint or geom.GetGeometryType() == ogr.wkbPoint25D):
        print "%s %.6f, %.6f" % (f['Comment'], geom.GetX(), geom.GetY())
        f['geom'] = geom
        f['X'] = geom.GetX()
        f['Y'] = geom.GetY()
        f['Z'] = f['GNSS_Heigh']   #? - should be able to get this from geom

        demGeom = geom.Clone()
        demGeom.Transform(gcpToDemTransform)
        demPixel, demLine = world2Pixel(demGeoTransform, demGeom.GetX(), demGeom.GetY())
        f['demX'] = demGeom.GetX()
        f['demY'] = demGeom.GetY()
        f['demPixel'] = demPixel
        f['demLine'] = demLine
        demP, demL = np.int(np.round(demPixel)), np.int(np.round(demLine))
        if demP >= 0 and demP < dem.shape[1] and demL >= 0 and demL < dem.shape[0]:
            f['demZ'] = dem[demL, demP]  #rather interpolate !
        else:
            f['demZ'] = -1.  #rather interpolate !

        # geom.Ta
    else:
        print "no point geometry/n"
    # gcpList.append(f)
    gcpDict[f['Comment']] = f
print ' '
ds = None

demPixel = np.array([gcp['demPixel'] for gcp in gcpDict.values()])
demLine = np.array([gcp['demLine'] for gcp in gcpDict.values()])
demZ = np.array([gcp['demZ'] for gcp in gcpDict.values()])
gcpZ = np.array([gcp['Z'] for gcp in gcpDict.values()])

pylab.hold('on')
pylab.plot(demPixel, demLine, 'kx')

pylab.figure()
pylab.plot(demZ, gcpZ, 'kx')
pylab.xlabel('DEM Z')
pylab.ylabel('GCP Z')