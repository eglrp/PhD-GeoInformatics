# %gui qt
# %matplotlib qt

# Analyse differences between GCP heights and DEM - is there something systematic?
import os
import gdal
import ogr
import osr
import numpy as np
import pylab
from matplotlib import patches
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import interpolate

os.environ['GDAL_DATA'] = "C:\ProgramData\Anaconda3\envs\py27\etc\conda\activate.d\\..\..\..\Library\share\gdal"
demFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/CGA/SUDEM/x3324c_15_15_L2a_crop.tif"
gcpGeoLocFile = "C:/Data/Development/Projects/PhD GeoInformatics/Docs/Misc/Baviaanskloof/BaviiaansPeCorrectedGcpMay2017Combined.shp"
gcpGeoidGeoLocFile = "C:/Data/Development/Projects/PhD GeoInformatics/Docs/Misc/Baviaanskloof/BaviiaansPeCorrectedGcpMay2017Combined_SaGeoid2010.shp"


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
    print "Open failed"

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
demInterp = interpolate.RectBivariateSpline(range(0, dem.shape[0]), range(0, dem.shape[1]), dem)
ds = None


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
            f['demZ'] = demInterp(demLine, demPixel)[0][0]  #  dem[demL, demP]  #rather interpolate !
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

pylab.figure()
pylab.imshow(dem)
pylab.hold('on')
pylab.plot(demPixel, demLine, 'kx')

pylab.figure()
pylab.plot(demZ, gcpZ, 'kx')
pylab.plot([0,demZ.max()], [0,demZ.max()], 'k--')
pylab.xlabel('DEM Z')
pylab.ylabel('GCP Z')


# Repeat for Geoid ref GCP's
ds = gdal.OpenEx(gcpGeoidGeoLocFile, gdal.OF_VECTOR)
if ds is None:
    print "Open failed./n"

lyr = ds.GetLayerByIndex(0)
lyr.ResetReading()
gcpSpatialRef = lyr.GetSpatialRef()

gcpToDemTransform = osr.CreateCoordinateTransformation(gcpSpatialRef, demSpatialRef)

# gcpList = []
gcpGeoidDict = {}
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
        f['Zgeoid'] = geom.GetZ() #f['GNSS_Heigh']   #? - should be able to get this from geom
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
            f['demZ'] = demInterp(demLine, demPixel)[0][0]  #  dem[demL, demP]  #rather interpolate !
        else:
            f['demZ'] = -1.  #ERROR

        # geom.Ta
    else:
        print "no point geometry/n"
    # gcpList.append(f)
    gcpGeoidDict[f['Comment']] = f
print ' '
ds = None

demPixel = np.array([gcp['demPixel'] for gcp in gcpGeoidDict.values()])
demLine = np.array([gcp['demLine'] for gcp in gcpGeoidDict.values()])
demZ = np.array([gcp['demZ'] for gcp in gcpGeoidDict.values()])
gcpZ = np.array([gcp['Z'] for gcp in gcpGeoidDict.values()])
gcpZgeoid = np.array([gcp['Zgeoid'] for gcp in gcpGeoidDict.values()])
gcpX = np.array([gcp['X'] for gcp in gcpGeoidDict.values()])
gcpY = np.array([gcp['Y'] for gcp in gcpGeoidDict.values()])

print 'WGS84 RMS: %f' % (np.sqrt(((demZ - gcpZ)**2).mean()))
print 'GEOID RMS: %f' % (np.sqrt(((demZ - gcpZgeoid)**2).mean()))

pylab.figure()
pylab.imshow(dem)
pylab.hold('on')
pylab.plot(demPixel, demLine, 'kx')

pylab.figure()
pylab.plot(demZ, gcpZ, 'bx')
pylab.plot(demZ, gcpZgeoid, 'rx')
pylab.legend(['WGS84', 'GEOID'])
pylab.plot([0, demZ.max()], [0, demZ.max()], 'k--')
pylab.xlabel('DEM Z')
pylab.ylabel('GCP Z')

#plot the errors from right to left and top down to look for patterns
pylab.figure()
pylab.subplot(2,1,1)
pylab.plot(demLine, demZ-gcpZ, 'bx')
pylab.plot(demLine, demZ-gcpZgeoid, 'rx')
pylab.ylabel('Z error')
pylab.xlabel('Y')
pylab.grid()
pylab.legend(['WGS84', 'GEOID'])
pylab.subplot(2,1,2)
pylab.plot(demPixel, demZ-gcpZ, 'bx')
pylab.plot(demPixel, demZ-gcpZgeoid, 'rx')
pylab.ylabel('Z error')
pylab.xlabel('X')
pylab.grid()
pylab.legend(['WGS84', 'GEOID'])

#plot error bars where they occur in 3d
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
for (x,y,z,zg) in zip(gcpX, gcpY, demZ, gcpZ):
    ax.plot3D([x,x], [y,y], [z, zg], 'b')
    ax.plot3D([x], [y], [z], 'go')
    ax.plot3D([x], [y], [zg], 'ro')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
for (x,y,z,zg) in zip(gcpX, gcpY, demZ, gcpZgeoid):
    ax.plot3D([x,x], [y,y], [z, zg], 'b')
    ax.plot3D([x], [y], [z], 'go')
    ax.plot3D([x], [y], [zg], 'ro')

xyScale = 1. # 46000./10000.
xgrid, ygrid = np.meshgrid(xyScale*np.arange(0, dem.shape[1]), xyScale*np.arange(0, dem.shape[0]))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.plot_wireframe(xgrid, ygrid, dem)
ax.scatter(demPixel*xyScale, demLine*xyScale, gcpZ, 'z', 100, 'r')
# ax.set_aspect('equal')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_xlim([0, xgrid.max()])
ax.set_ylim([0, ygrid.max()])
ax.set_zlim([dem.min(), dem.max()])
ax.set_alpha(0.5)

# Fit a surface through the errors to see if it is systematic
from numpy.polynomial import polynomial
import numpy as np


def polyfit2d(x, y, f, deg):
    x = np.asarray(x)
    y = np.asarray(y)
    f = np.asarray(f)
    deg = np.asarray(deg)
    vander = polynomial.polyvander2d(x, y, deg)
    vander = vander.reshape((-1,vander.shape[-1]))
    f = f.reshape((vander.shape[0],))
    c = np.linalg.lstsq(vander, f)[0]
    return c.reshape(deg+1)

off = gcpZgeoid - demZ
p = polyfit2d(gcpX, gcpY, off, [2, 2])
offHat = polynomial.polyval2d(gcpX, gcpY, p)

print np.abs(off).mean()
print np.abs(offHat-off).mean()
print np.std(off)
print np.std(offHat-off)

xgrid, ygrid = np.meshgrid(np.arange(gcpX.min(), gcpX.max(), .00005), np.arange(gcpY.min(), gcpY.max(), .00005))

polyZ = polynomial.polyval2d(xgrid, ygrid, p)

pylab.figure()
pylab.imshow(polyZ, extent=(gcpX.min(), gcpX.max(), gcpY.max(), gcpY.min()))
# plt.imshow(zz, extent=(x.min(), y.max(), x.max(), y.min()))
pylab.scatter(gcpX, gcpY, c=off)
pylab.plot(gcpX, gcpY, 'kx')


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_wireframe(xgrid, ygrid, polyZ)
ax.scatter(gcpX, gcpY, off, 'z', 10, 'r')
# ax.set_aspect('equal')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_alpha(0.5)
