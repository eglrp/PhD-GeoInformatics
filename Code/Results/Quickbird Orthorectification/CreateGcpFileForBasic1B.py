
import gdal
import ogr
import osr


# For basic level 1b image
# This makes a text file for using in PCI that has the Geo co-ords and image locs of my GCP's
# Note that I (stupidly) first collected the image locs in a projected co-ord system, so this needs
# to be converted back to image pixel locs

# file containing actual GCP locs in WGS84
doEllipsoidalHeight = True
if doEllipsoidalHeight:  # ellipsoidal heights
    gcpGeoLocFile = "C:/Data/Development/Projects/PhD GeoInformatics/Docs/Misc/Baviaanskloof/BaviiaansPeCorrectedGcpMay2017Combined.shp"
else:
    gcpGeoLocFile = "C:/Data/Development/Projects/PhD GeoInformatics/Docs/Misc/Baviaanskloof/BaviiaansPeCorrectedGcpMay2017Combined_SaGeoid2010.shp"
# file containing image locations of GCP locs in UTM 35S
gcpImLocFile = "C:/Data/Development/Projects/PhD GeoInformatics/Docs/Misc/Baviaanskloof/Quickbird-056844553010_01-Gcp.shp"

# imFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Digital Globe/056549293010_01/PanSharpen/03NOV18082012-P2AS_R1C12-056549293010_01_P001_GdalPanSharp.tif"
# imFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Digital Globe/056549293010_01/056549293010_01_P001_PAN/03NOV18082012-P2AS_R1C1-056549293010_01_P001.TIF"
# imFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Digital Globe/056549293010_01/056549293010_01_P001_PAN/Assembled/03NOV18082012-P2AS_R1C12-056549293010_01_P001.TIF"
# nb the im is in wgs84
imFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Digital Globe/056844553010_01/056844553010_01_P001_PAN/03NOV18082012-P1BS-056844553010_01_P001.TIF"


# read in Geo CoOrds of GCP's (and height!)
ds = ogr.Open(gcpGeoLocFile)
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
        if doEllipsoidalHeight:
            f['Z'] = f['GNSS_Heigh']   #? - should be able to get this from geom
        else:
            f['Z'] = geom.GetZ()        # this has been xformed from elippsoidal to Sa Geoid 2010

    else:
        print "no point geometry/n"
    # gcpList.append(f)
    gcpDict[f['Comment']] = f
print ' '
ds = None


# read in UTM CoOrds of GCP's image locs
ds = ogr.Open(gcpImLocFile)
if ds is None:
    print "Open failed./n"

lyr = ds.GetLayerByIndex(0)
lyr.ResetReading()
spatialRefImLoc = lyr.GetSpatialRef()

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
ds = gdal.Open(imFile)
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

geotransform = [0., 0., 0., 0., 0., 0.]
if True:  # hack for projectionless image
    geotransform[0] = 24.30712463
    geotransform[3] = -33.59966141
    geotransform[1] = (24.50332750 - 24.30712463) / ds.RasterXSize
    geotransform[5] = (-33.59966141 - -33.74584942) / ds.RasterYSize
# 24.3071246305866104,-33.7458494236889948 : 24.5033274972737871,-33.5996614073278437
# extents from qgis
# u'24.30712463 -33.74584942, 24.30712463 -33.59966141, 24.50332750 -33.59966141, 24.50332750 -33.74584942, 24.30712463 -33.74584942'

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
  line = abs((y - ulY) / yDist)
  return (pixel, line)

# xform = osr.CoordinateTransformation(spatialRefImLoc, osr.SpatialReference(ds.GetProjectionRef()))

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
if doEllipsoidalHeight:
    gcpFile = open("C:/Data/Development/Projects/PhD GeoInformatics/Docs/Misc/Baviaanskloof/Quickbird-056844553010-Gcp.txt", 'w')
else:
    gcpFile = open("C:/Data/Development/Projects/PhD GeoInformatics/Docs/Misc/Baviaanskloof/Quickbird-056844553010-Gcp-SaGeoid2010.txt", 'w')
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



if False:
    #arcmap format
    transform = osr.CoordinateTransformation(spatialRef, spatialRefIm)
    if doEllipsoidalHeight:
        gcpFile = open("C:/Data/Development/Projects/PhD GeoInformatics/Docs/Misc/Baviaanskloof/Quickbird-056844553010-Gcp-Arcmap.txt", 'w')
    else:
        gcpFile = open("C:/Data/Development/Projects/PhD GeoInformatics/Docs/Misc/Baviaanskloof/Quickbird-056844553010-Gcp-Arcmap-SaGeoid2010.txt", 'w')

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
