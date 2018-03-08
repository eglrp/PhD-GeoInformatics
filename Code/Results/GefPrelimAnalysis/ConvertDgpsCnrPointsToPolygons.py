import gdal
import ogr
import numpy as np
import osr
import pylab
from scipy import stats as stats
from matplotlib import patches
import matplotlib.pyplot as plt
from scipy import ndimage as ndimage

# Python Imaging Library imports
from PIL import Image
from PIL import ImageDraw

import os
os.environ["GDAL_DATA"] = 'C:\OSGeo4W64\share\gdal'


# |layerid=0|subset="Comment" LIKE 'H%'
correctedShapeFileNames = [
    'C:\Data\Development\Projects\PhD GeoInformatics\Data\GEF Field Trial\DGPS Sept 2017\Corrected\Point_ge.shp',
    'C:\Data\Development\Projects\PhD GeoInformatics\Data\GEF Sampling\DGPS Dec 2017\Corrected\Point_ge.shp']

outShapeFileName = 'C:\Data\Development\Projects\PhD GeoInformatics\Data\GEF Sampling\GefFieldSamplingGt.shp'

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

## read dgps points into dictionary
dgpsSrs = osr.SpatialReference()
dgpsDict = {}
for correctedShapeFileName in correctedShapeFileNames:
    ds = gdal.OpenEx(correctedShapeFileName, gdal.OF_VECTOR)
    if ds is None:
        print "Open failed./n"

    lyr = ds.GetLayerByIndex(0)
    lyr.ResetReading()
    if lyr.GetSpatialRef() is not None and lyr.GetSpatialRef() is not 0:
        dgpsSrs = lyr.GetSpatialRef()
    else:
        dgpsSrs.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

    for (i, feat) in enumerate(lyr):
        print '.',
        feat_defn = lyr.GetLayerDefn()
        f = {}
        for i in range(feat_defn.GetFieldCount()):
            field_defn = feat_defn.GetFieldDefn(i)
            f[field_defn.GetName()] = feat.GetField(i)
        key = "%s-%s" % (f['Datafile'][:-4], f['Comment'])

        geom = feat.GetGeometryRef()
        if geom is not None and (geom.GetGeometryType() == ogr.wkbPoint or geom.GetGeometryType() == ogr.wkbPointZM):
            # concat(left(Datafile,strpos(Datafile,'.cor')-1), '-',  Comment)
            print "%s: %.6f, %.6f" % (key, geom.GetX(), geom.GetY())
            f['geom'] = geom.Clone()
            f['X'] = geom.GetX()
            f['Y'] = geom.GetY()
        else:
            print "%s has no point geometry" % (key)
            break
        # gcpList.append(f)
        f['PlotName'] = f['Datafile'][:-4]
        dgpsDict[key] = f
    print ' '

ds = None

## Create output shapefile
# set up the shapefile driver

driver = ogr.GetDriverByName("ESRI Shapefile")
ds = driver.CreateDataSource(outShapeFileName)
# osr.SpatialReference().Impo(4326)
layer = ds.CreateLayer(outShapeFileName[:-4], dgpsSrs, ogr.wkbMultiPolygon)
# Add the fields we're interested in
# field_name = ogr.FieldDefn("Name", ogr.OFTString)
# field_name.SetWidth(64)
layer.CreateField(ogr.FieldDefn("Name", ogr.OFTString))
layer.CreateField(ogr.FieldDefn("WoodyCS", ogr.OFTReal))


plotNames = np.array([f['PlotName'] for f in dgpsDict.values()])

for plotName in np.unique(plotNames):
    idx = plotNames == plotName
    plotPoints = np.array(dgpsDict.values())[idx]
    comments = [plotPoint['Comment'] for plotPoint in plotPoints]
    cnrIdx = [comment.startswith('H') for comment in comments]
    plotCnrs = plotPoints[cnrIdx]
    print plotName,
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetField("Name", plotName)
    feature.SetField("WoodyCS", 0.)

    # plotLinRing = ogr.Geometry(ogr.wkbLinearRing)
    plotGeomColl = ogr.Geometry(ogr.wkbGeometryCollection)
    for plotCnr in plotCnrs:
        plotGeomColl.AddGeometry(plotCnr['geom'])
        # plotLinRing.AddPoint(plotCnr['X'], plotCnr['Y'], plotCnr['geom'].GetZ())
        print '.',
    print
    plotPoly = plotGeomColl.ConvexHull()
        # ogr.Geometry(ogr.wkbPolygon)
    # plotPoly.AddGeometry(plotLinRing)
    feature.SetGeometry(plotPoly)
    layer.CreateFeature(feature)
    # Dereference the feature
    feature = None

ds = None



if False:
    # Process the text file and add the attributes and features to the shapefile
    for row in reader:
        # create the feature
        feature = ogr.Feature(layer.GetLayerDefn())
        # Set the attributes using the values from the delimited text file
        feature.SetField("Name", row['Name'])
        feature.SetField("Region", row['Region'])
        feature.SetField("Latitude", row['Latitude'])
        feature.SetField("Longitude", row['Longitude'])
        feature.SetField("Elevation", row['Elev'])

        # create the WKT for the feature using Python string formatting
        wkt = "POINT(%f %f)" % (float(row['Longitude']), float(row['Latitude']))

        # Create the point from the Well Known Txt
        point = ogr.CreateGeometryFromWkt(wkt)

        # Set the feature geometry using the point
        feature.SetGeometry(point)
        # Create the feature in the layer (shapefile)
        layer.CreateFeature(feature)
        # Dereference the feature
        feature = None