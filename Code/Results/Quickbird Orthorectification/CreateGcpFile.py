
import gdal
import ogr
import numpy as np
import osr
import pylab
from scipy import stats as stats
from collections import OrderedDict


# This makes a text file for using in PCI that has the Geo co-ords and image locs of my GCP's
# Note that I (stupidly) first collected the image locs in a projected co-ord system, so this needs
# to be converted back to image pixel locs

# file containing actual GCP locs in WGS84
gcpGeoLocFile = "C:/Data/Development/Projects/PhD GeoInformatics/Docs/Misc/Baviaanskloof/BaviiaansPeCorrectedGcpMay2017Combined.shp"

# file containing image locations of GCP locs in UTM 35S
gcpImLocFile = "C:/Data/Development/Projects/PhD GeoInformatics/Docs/Misc/Baviaanskloof/QuickbirdGcpUtmZ35SZ.shp"

# imFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Digital Globe/056549293010_01/PanSharpen/03NOV18082012-P2AS_R1C12-056549293010_01_P001_GdalPanSharp.tif"
imFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Digital Globe/056549293010_01/056549293010_01_P001_PAN/03NOV18082012-P2AS_R1C1-056549293010_01_P001.TIF"
# read in Geo CoOrds of GCP's (and height!)
ds = gdal.OpenEx(gcpGeoLocFile, gdal.OF_VECTOR)
if ds is None:
    print "Open failed./n"

lyr = ds.GetLayerByIndex(0)
lyr.ResetReading()
spatialRef = lyr.GetSpatialRef()


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
    else:
        print "no point geometry/n"
    # gcpList.append(f)
    gcpDict[f['Comment']] = f
print ' '
ds = None


# read in UTM CoOrds of GCP's image locs
ds = gdal.OpenEx(gcpImLocFile, gdal.OF_VECTOR)
if ds is None:
    print "Open failed./n"

lyr = ds.GetLayerByIndex(0)
lyr.ResetReading()
spatialRefIm = lyr.GetSpatialRef()

imGcpList = []
imGcpDict = {}

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
        f['Xi'] = 0   # im locs
        f['Yi'] = 0
    else:
        print "no point geometry/n"
    # imGcpList.append(f)
    imGcpDict[f['Comment']] = f
print ' '
ds = None

# open image to get extents in proj co-ords
# read in Geo CoOrds of GCP's (and height!)
ds = gdal.OpenEx(imFile, gdal.OF_RASTER)
if ds is None:
    print "Open failed./n"

print 'Driver: ', ds.GetDriver().ShortName,'/', \
      ds.GetDriver().LongName
print 'Size is ',ds.RasterXSize,'x',ds.RasterYSize, \
      'x',ds.RasterCount
print 'Projection is ',ds.GetProjection()
geotransform = ds.GetGeoTransform()
if not geotransform is None:
    print 'Origin = (',geotransform[0], ',',geotransform[3],')'
    print 'Pixel Size = (',geotransform[1], ',',geotransform[5],')'

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

for f in imGcpDict.values():
    (xi, yi) = world2Pixel(geotransform, f['X'], f['Y'])
    f['Xi'] = xi
    f['Yi'] = yi
    print "%s %i, %i" % (f['Comment'], xi, yi)
    if xi < 0 or xi > ds.RasterXSize:
        print '-------------------------------xi out of bounds'
    if yi < 0 or yi > ds.RasterYSize:
        print '--------------------------------yi out of bounds'


# Create text file

# Character descriptions:
# I: the GCP's identification number.
#
# X: the geocoded X coordinate.
#
# Y: the geocoded Y coordinate.
#
# P: the pixel location of the GCP on the uncorrected image file.
#
# L: the line location of the GCP on the uncorrected image file.
#
# E: the elevation of the GCP.
#
# D: data to be ignored.
#
# For example, the string IXYE represents the layout where I is the GCP's identification, X and Y are the GCP's geocoded x and y coordinates, and E is the GCP's elevation.
#
# The line in the text file and the Format must match. Lines that do not match the Format are ignored.


#pci format
gcpFile = open("C:/Data/Development/Projects/PhD GeoInformatics/Docs/Misc/Baviaanskloof/QuickbirdGcp3.txt", 'w')
gcpFile.write("I\tP\tL\tX\tY\tE\n")
id = 1
for key in imGcpDict:
    # key = f['Comment']
    fgcp = gcpDict[key]
    fim = imGcpDict[key]
    #print "%s %.6f, %.6f" % (f['Comment'], geom.GetX(), geom.GetY())
    gcpFile.write("%s\t%.3f\t%.3f\t%.9f\t%.9f\t%.2f\n" % (key, fim['Xi'], fim['Yi'], fgcp['X'], fgcp['Y'], fgcp['Z']))
    id = id + 1
    print "%s" % (key)
gcpFile.close()


#arcmap format
transform = osr.CoordinateTransformation(spatialRef, spatialRefIm)
gcpFile = open("C:/Data/Development/Projects/PhD GeoInformatics/Docs/Misc/Baviaanskloof/QuickbirdGcp_Arcmap.txt", 'w')
for key in imGcpDict:
    # key = f['Comment']
    fgcp = gcpDict[key]
    fim = imGcpDict[key]
    # print "%s %.6f, %.6f" % (f['Comment'], geom.GetX(), geom.GetY())
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(fgcp['X'], fgcp['Y'])
    point.Transform(transform)
    gcpFile.write("%.3f %.3f %.3f %.3f\n" % (fim['X'], fim['Y'], point.GetX(), point.GetY()))
    id = id + 1
    print "%s" % (key)
gcpFile.close()
