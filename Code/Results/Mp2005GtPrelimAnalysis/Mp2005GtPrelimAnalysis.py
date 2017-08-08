import gdal
import ogr
import numpy as np
import osr
import pylab
from scipy import stats as stats
from collections import OrderedDict

# take a rough first look for feature correlations with cs

# cs gt file - NB the locs in this file are rounded to 5 decimal places which is only accurate to something like 2m

# we should combine with gps file which is more accurate
csGtFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Misc/BMR Carbon Stocks/abf_agc_191_plots.shp"
# the file below contains only the plot locations but to greater accuracy than the above file
csGtGpsFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Misc/BMR Carbon Stocks/gps_coords_191plots.shp"

# file containing image locations of GCP locs in UTM 35S
# imFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Digital Globe/056549293010_01/Ortho/R1C12-GdalPanSharp-ArcGcpWarp.tif"
# imFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Digital Globe/056844553010_01/PCI Output/TOA and Haze/TOACorrected_056844553010_01_P001_OrthoPanSharpen_05644015.tif"
imFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Digital Globe/056844553010_01/PCI Output\ATCOR1/ATCORCorrected_056844553010_01_P001_OrthoPanSharpen_05644032.tif"
# read in cs gt
ds = gdal.OpenEx(csGtFile, gdal.OF_VECTOR)
if ds is None:
    print "Open failed./n"

lyr = ds.GetLayerByIndex(0)
lyr.ResetReading()
csGtSpatialRef = lyr.GetSpatialRef()

# gcpList = []
csGtDict = {}
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
        print "%s %.6f, %.6f" % (f['PLOT'], geom.GetX(), geom.GetY())
        f['geom'] = geom
        f['X'] = geom.GetX()
        f['Y'] = geom.GetY()
    else:
        print "no point geometry/n"
    # gcpList.append(f)
    csGtDict[f['PLOT']] = f
print ' '

ds = None

# read in cs gt locations
ds = gdal.OpenEx(csGtGpsFile, gdal.OF_VECTOR)
if ds is None:
    print "Open failed./n"

lyr = ds.GetLayerByIndex(0)
lyr.ResetReading()
csGtSpatialRef = lyr.GetSpatialRef()

csGtGpsDict = {}
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
        print "%s %.6f, %.6f" % (f['CODE'], geom.GetX(), geom.GetY())
        f['geom'] = geom
        f['X'] = geom.GetX()
        f['Y'] = geom.GetY()
    else:
        print "no point geometry/n"
    # gcpList.append(f)
    csGtGpsDict[f['CODE']] = f
print ' '

ds = None

# copy across gps locs to csGtDict
# [s for s in csGtGpsDict.keys() if "GH" in str(s)]
# [s for s in csGtDict.keys() if "GH" in str(s)]
# [s for s in csGtGpsDict.keys() if "RHOL" in str(s)]
#
gpsInd = (np.zeros((1, csGtGpsDict.__len__())))

for k in csGtDict.keys():
    # fixing some quirks in the keys
    gpsK = k
    gpsK = str(gpsK).replace('KADSS', 'KASS')
    gpsK = str(gpsK).replace('GHOL1_', 'RHOL1_')
    if gpsK == 'GHOL3_5':
        gpsK = 'GHOL3_'
    # gpsK = str(gpsK).replace('GHOL3_', 'GHOL3_5')
    if csGtGpsDict.has_key(gpsK):
        csGtDict[k]['X'] = csGtGpsDict[gpsK]['X']
        csGtDict[k]['Y'] = csGtGpsDict[gpsK]['Y']
        # print k
        gpsInd[0, csGtGpsDict.keys().index(gpsK)] = True
    else:
        print k + " not found in gps locs - " + gpsK

np.array(csGtGpsDict.keys())[np.logical_not(gpsInd[0])]

# Read in the image

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


def extract_features(imbuf):
    imbuf = np.float64(imbuf)
    s = np.sum(imbuf, 2)
    cn = imbuf/np.tile(s[:,:,None], (1, 1, imbuf.shape[2]))
    ndvi = (imbuf[:,:,3] - imbuf[:,:,0])/(imbuf[:,:,3] + imbuf[:,:,0])
    ir_rat = imbuf[:,:,3]/imbuf[:,:,0]
    feat = {}
    feat['r_n'] = cn[:,:,0].mean()
    feat['g_n'] = cn[:,:,1].mean()
    feat['b_n'] = cn[:,:,2].mean()
    feat['ir_n'] = cn[:,:,3].mean()
    feat['NDVI'] = ndvi.mean()
    feat['ir_rat'] = ir_rat.mean()
    feat['i'] = (s/np.prod(s.shape[0:2])).mean()
    feat['i_std'] = (s/np.prod(s.shape[0:2])).std()
    return feat

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

#

transform = osr.CoordinateTransformation(csGtSpatialRef, osr.SpatialReference(ds.GetProjection()))
i = 0
winSize = (5, 5)
plotDict = {}
plotTagcDict = {}
for plot in csGtDict.values():
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(plot['X'], plot['Y'])
    point.Transform(transform)  # xform into im projection
    (pixel, line) = world2Pixel(geotransform, point.GetX(), point.GetY())
    # not all the point fall inside the image
    if pixel >= 0 and line >=0 and pixel < ds.RasterXSize and line < ds.RasterYSize:
        imbuf = np.zeros((winSize[0], winSize[1], 4), dtype=float)
        for b in range(1,5):
            imbuf[:, :, b-1] = ds.GetRasterBand(b).ReadAsArray(np.int(np.round(pixel))-(winSize[0]-1)/2,
                                                               np.int(np.round(line))-(winSize[1]-1)/2,
                                                               winSize[0], winSize[1])
            # imbuf[:, :, b-1] = ds.GetRasterBand(b).ReadAsArray(np.int(np.round(pixel))-(winSize[0]-1),
            #                                                    np.int(np.round(line))-(winSize[1]-1),
            #                                                    winSize[0], winSize[1])
            # imbuf[:, :, b-1] = ds.GetRasterBand(b).ReadAsArray(np.int(np.round(pixel)),
            #                                                    np.int(np.round(line)),
            #                                                    winSize[0], winSize[1])
            if np.all(imbuf==0):
                print "imbuf zero"
            feat = extract_features(imbuf)
            feat['TAGC'] = csGtDict[plot['PLOT']]['TAGC']
            plotDict[plot['PLOT']] = feat
            # plotTagcDict[plot['PLOT']] = csGtDict[plot['PLOT']]['TAGC']
        print plot['PLOT']
        i = i +1
    else:
        print "x-" + plot['PLOT']
print i


ndvi = np.array([plot['NDVI'] for plot in plotDict.values()])
gn = np.array([plot['g_n'] for plot in plotDict.values()])
std = np.array([plot['i_std'] for plot in plotDict.values()])
ir_rat = np.array([plot['ir_rat'] for plot in plotDict.values()])
tagc = np.array([plot['TAGC'] for plot in plotDict.values()])
plotNames = plotDict.keys()

pylab.figure()
pylab.subplot(2,2,1)
pylab.plot(ndvi, tagc, 'kx')
pylab.xlabel('NDVI')
pylab.ylabel('TAGC')
pylab.subplot(2,2,2)
pylab.plot(gn, tagc, 'kx')
pylab.xlabel('gn')
pylab.ylabel('TAGC')
pylab.subplot(2,2,3)
pylab.plot(std, tagc, 'kx')
pylab.xlabel('i_std')
pylab.ylabel('TAGC')
pylab.subplot(2,2,4)
pylab.plot(ir_rat, tagc, 'kx')
pylab.xlabel('ir_rat')
pylab.ylabel('TAGC')


idx = np.array(['ST' in str(s) and not 'DST' in str(s) for s in plotNames])

pylab.figure()
pylab.subplot(2,2,1)
pylab.plot(ndvi[idx], tagc[idx], 'kx')
pylab.xlabel('NDVI')
pylab.ylabel('TAGC')
pylab.subplot(2,2,2)
pylab.plot(gn[idx], tagc[idx], 'kx')
pylab.xlabel('gn')
pylab.ylabel('TAGC')
pylab.subplot(2,2,3)
pylab.plot(std[idx], tagc[idx], 'kx')
pylab.xlabel('i_std')
pylab.ylabel('TAGC')
pylab.subplot(2,2,4)
pylab.plot(ir_rat[idx], tagc[idx], 'kx')
pylab.xlabel('ir_rat')
pylab.ylabel('TAGC')


ds = None
