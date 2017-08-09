import gdal
import ogr
import numpy as np
import osr
import pylab
from scipy import stats as stats
from matplotlib import patches

from collections import OrderedDict

# take a rough first look for feature correlations with cs

# cs gt file - NB the locs in this file are rounded to 5 decimal places which is only accurate to something like 2m

# we should combine with gps file which is more accurate
csGtFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Misc/BMR Carbon Stocks/abf_agc_191_plots.shp"
# the file below contains only the plot locations but to greater accuracy than the above file
csGtGpsFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Misc/BMR Carbon Stocks/gps_coords_191plots.shp"

# file containing image locations of GCP locs in UTM 35S
# imFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Digital Globe/056549293010_01/Ortho/R1C12-GdalPanSharp-ArcGcpWarp.tif"
# imFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Digital Globe/056844553010_01/PCI Output/ATCOR1/ATCORCorrected_056844553010_01_P001_OrthoPanSharpen_05644032.tif"
imFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Digital Globe/056844553010_01/PCI Output/TOA and Haze/TOACorrected_056844553010_01_P001_OrthoPanSharpen_05644015.tif"
# imFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Digital Globe/056844553010_01/PCI Output/Separate Pan and MS/TOA/PanSharpToaOrtho.tif"
# imFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Digital Globe/056844553010_01/PCI Output/Separate Pan and MS/ATCOR1/PansharpAtcorOrtho.tif"
# imFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Digital Globe/056844553010_01/PCI Output/056844553010_01_P001_OrthoPanSharpen.tif"

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
    imbuf = np.float64(imbuf)/100.   # values are in percent?
    s = np.sum(imbuf, 2)
    cn = imbuf/np.tile(s[:,:,None], (1, 1, imbuf.shape[2]))
    b_i = 0
    g_i = 1
    r_i = 2
    ir_i = 3
    ndvi = (imbuf[:,:,ir_i] - imbuf[:,:,r_i])/(imbuf[:,:,ir_i] + imbuf[:,:,r_i])
    ir_rat = imbuf[:,:,ir_i]/imbuf[:,:,r_i]
    L = 0.02
    savi = (1 + L)*(imbuf[:,:,ir_i] - imbuf[:,:,r_i])/(L + imbuf[:,:,ir_i] + imbuf[:,:,r_i])
    feat = {}
    feat['r_n'] = cn[:,:,r_i].mean()
    feat['g_n'] = cn[:,:,g_i].mean()
    feat['b_n'] = cn[:,:,b_i].mean()
    feat['ir_n'] = cn[:,:,ir_i].mean()
    feat['NDVI'] = ndvi.mean()
    feat['SAVI'] = savi.mean()
    feat['ir_rat'] = ir_rat.mean()
    feat['i'] = (s/np.prod(s.shape[0:2])).mean()
    feat['i_std'] = (s/np.prod(s.shape[0:2])).std()
    feat['NDVI_std'] = ndvi.std()

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
winSize = (8, 8)
plotDict = {}
plotTagcDict = {}
class_labels = ['OL','DST','ST']
for plot in csGtDict.values():
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(plot['X'], plot['Y'])
    point.Transform(transform)  # xform into im projection
    (pixel, line) = world2Pixel(geotransform, point.GetX(), point.GetY())
    # not all the point fall inside the image
    if pixel >= 0 and line >=0 and pixel < ds.RasterXSize and line < ds.RasterYSize:
        imbuf = np.zeros((winSize[0], winSize[1], 4), dtype=float)
        for b in range(1,5):
            # imbuf[:, :, b-1] = ds.GetRasterBand(b).ReadAsArray(np.int(np.round(pixel))-(winSize[0]-1)/2,
            #                                                    np.int(np.round(line))-(winSize[1]-1)/2,
            #                                                    winSize[0], winSize[1])
            # imbuf[:, :, b-1] = ds.GetRasterBand(b).ReadAsArray(np.int(np.round(pixel))-(winSize[0]-1),
            #                                                    np.int(np.round(line))-(winSize[1]-1),
            #                                                    winSize[0], winSize[1])
            # imbuf[:, :, b-1] = ds.GetRasterBand(b).ReadAsArray(np.int(np.round(pixel)),
            #                                                    np.int(np.round(line)),
            #                                                    winSize[0], winSize[1])
            # imbuf[:, :, b-1] = ds.GetRasterBand(b).ReadAsArray(np.int(np.round(pixel))-(winSize[0]-1),
            #                                                    np.int(np.round(line)),
            #                                                    winSize[0], winSize[1])
            # this option performs best
            imbuf[:, :, b-1] = ds.GetRasterBand(b).ReadAsArray(np.int(np.round(pixel)),
                                                               np.int(np.round(line))-(winSize[0]-1),
                                                               winSize[0], winSize[1])
            if np.all(imbuf==0):
                print "imbuf zero"
            feat = extract_features(imbuf)
            fields = ['TAGC', 'Z_P_AFRA', 'ALLOMETRY', 'HERB', 'LITTER']
            for f in fields:
                feat[f] = csGtDict[plot['PLOT']][f]
            if 'DST' in plot['PLOT']:
                ci = 1
            elif 'ST' in plot['PLOT']:
                ci = 2
            elif 'OL' in plot['PLOT']:
                ci = 0
            else:
                ci = 2
            feat['classi'] = ci
            feat['class'] = class_labels[ci]
            feat['thumbnail'] = imbuf
            plotDict[plot['PLOT']] = feat
            # plotTagcDict[plot['PLOT']] = csGtDict[plot['PLOT']]['TAGC']
        print plot['PLOT']
        i = i +1
    else:
        print "x-" + plot['PLOT']
print i

#
# ndvi = np.array([plot['NDVI'] for plot in plotDict.values()])
# gn = np.array([plot['g_n'] for plot in plotDict.values()])
# std = np.array([plot['i_std'] for plot in plotDict.values()])
# ir_rat = np.array([plot['ir_rat'] for plot in plotDict.values()])
# tagc = np.array([plot['TAGC'] for plot in plotDict.values()])
plotNames = plotDict.keys()


yfields = ['TAGC', 'Z_P_AFRA', 'ALLOMETRY', 'HERB', 'LITTER']
xfields = ['NDVI', 'g_n', 'b_n', 'ir_n', 'ir_rat', 'NDVI_std']

pylab.figure()
x = np.array([plot['NDVI'] for plot in plotDict.values()])
for yi, yf in enumerate(yfields):
    pylab.subplot(2, 3, yi+1)
    y = np.array([plot[yf] for plot in plotDict.values()])
    (slope, intercept, r, p, stde) = stats.linregress(x, y)
    pylab.plot(x, y, 'kx')
    # pylab.plot(mPixels[::step], dRawPixels[::step], color='k', marker='.', linestyle='', markersize=.5)
    xl = pylab.gca().get_xlim()
    yl = pylab.gca().get_ylim()
    pylab.text((xl[0] + np.diff(xl)*0.05)[0], (yl[0] + np.diff(yl)*0.8)[0], str.format('$R^2$ = {0:.2f}',
                                                                                       np.round(r**2, 2)))
    pylab.xlabel('NDVI')
    pylab.ylabel(yf)


pylab.figure()
y = np.log10([plot['TAGC'] for plot in plotDict.values()])
for xi, xf in enumerate(xfields):
    pylab.subplot(2, 3, xi+1)
    x = np.array([plot[xf] for plot in plotDict.values()])
    (slope, intercept, r, p, stde) = stats.linregress(x, y)
    pylab.plot(x, y, 'kx')
    # pylab.plot(mPixels[::step], dRawPixels[::step], color='k', marker='.', linestyle='', markersize=.5)
    xl = pylab.gca().get_xlim()
    yl = pylab.gca().get_ylim()
    pylab.text((xl[0] + np.diff(xl)*0.05)[0], (yl[0] + np.diff(yl)*0.8)[0], str.format('$R^2$ = {0:.2f}',
                                                                                       np.round(r**2, 2)))
    pylab.xlabel(xf)
    pylab.ylabel('TAGC')


pylab.figure()
x = np.array([plot['NDVI'] for plot in plotDict.values()])
classi = np.array([plot['classi'] for plot in plotDict.values()])
class_lab = np.array([plot['class'] for plot in plotDict.values()])
colours = ['r','m','b']
plot_names = plotDict.keys()
for yi, yf in enumerate(yfields):
    pylab.subplot(2, 3, yi+1)
    y = np.array([plot[yf] for plot in plotDict.values()])
    (slope, intercept, r, p, stde) = stats.linregress(x, y)
    for idx,col in enumerate(colours):
        class_idx = classi == idx
        pylab.plot(x[class_idx], y[class_idx], col + 'x')

        for xx,yy,ll in zip(x[class_idx], y[class_idx], np.array(plot_names)[class_idx]):
            pylab.text(xx+.001, yy+.001, ll, fontdict={'size':6, 'color': col})
    # pylab.plot(mPixels[::step], dRawPixels[::step], color='k', marker='.', linestyle='', markersize=.5)

    xl = pylab.gca().get_xlim()
    yl = pylab.gca().get_ylim()
    pylab.text((xl[0] + np.diff(xl)*0.05)[0], (yl[0] + np.diff(yl)*0.8)[0], str.format('$R^2$ = {0:.2f}',
                                                                                       np.round(r**2, 2)))
    pylab.xlabel('NDVI')
    pylab.ylabel(yf)
    pylab.legend(class_labels)


#with images
pylab.figure()
x = np.array([plot['NDVI'] for plot in plotDict.values()])
classi = np.array([plot['classi'] for plot in plotDict.values()])
class_lab = np.array([plot['class'] for plot in plotDict.values()])
colours = ['r','m','b']
plot_names = plotDict.keys()
for yi, yf in enumerate(['TAGC']):
    ax = pylab.subplot(1, 1, yi+1)
    y = np.log10(np.array([plot[yf] for plot in plotDict.values()]))
    (slope, intercept, r, p, stde) = stats.linregress(x, y)
    ylim = [np.min(y), np.max(y)]
    xlim = [np.min(x), np.max(x)]
    xd = np.diff(xlim)[0]
    yd = np.diff(ylim)[0]
    pylab.axis([np.min(x), np.max(x), np.min(y), np.max(y)])
    pylab.hold('on')
    #pylab.plot(x, y, 'k.')

    for idx,col in enumerate(colours):
        class_idx = classi == idx
        # pylab.plot(x[class_idx], y[class_idx], col + 'x')
        for xx,yy,ll in zip(x[class_idx], y[class_idx], np.array(plot_names)[class_idx]):
            imbuf = plotDict[ll]['thumbnail'].copy()
            imbuf[:,:,3] = imbuf[:,:,3]/1.5

            # pylab.text(xx+.001, yy+.001, ll, fontdict={'size':6, 'color': col})
            ims = 20.
            extent = [xx-xd/(2*ims), xx+xd/(2*ims), yy-yd/(2*ims), yy + yd/(2*ims)]
            pylab.imshow(imbuf[:,:,[2,1,0]]/30., extent=extent, aspect='auto') #zorder=-1,
            ax.add_patch(patches.Rectangle((xx-xd/(2*ims), yy-yd/(2*ims)), xd/ims, yd/ims, fill=False, edgecolor=col, linewidth=2.))
    # pylab.plot(mPixels[::step], dRawPixels[::step], color='k', marker='.', linestyle='', markersize=.5)

    xl = pylab.gca().get_xlim
    yl = pylab.gca().get_ylim()
    pylab.text((xl[0] + np.diff(xl)*0.05)[0], (yl[0] + np.diff(yl)*0.8)[0], str.format('$R^2$ = {0:.2f}',
                                                                                       np.round(r**2, 2)))
    pylab.xlabel('NDVI')
    pylab.ylabel(yf)
    pylab.legend(class_labels)

#todo x check if atcor refl vals are in % - yes they are
#todo x define ol, dst, st classes and visualise
#todo research "regression analysis" in general and in python
#todo plot with labels and or images to try and figure out what are problem plots and patterns
#todo show pattern of correlation with spatial accuracy and window size - somehow simulate different accuracy and plot size conditions
#todo make a separate feature of area around window - should give some insight into eg soil condition
#todo brownness or soil index
#todo make features that use classification within window eg perhaps ttl vegetation pixels, or mean soil colour and mean veg colour
#todo investigate visualisation with qgis/arc python libraries
#kko, kq are pristine

#
# import sys
# sys.path.append("C:/OSGeo4W64/apps/Python27/Lib/site-packages")
# sys.path.append("C:/OSGeo4W64/apps/qgis/python/qgis")
# sys.path.append("C:/OSGeo4W64/apps/qgis/python/")
# sys.path.append("C:/OSGeo4W64/bin/")

# start pycharm with batch file
if False:
    %gui qt
    import qgis.core
    import qgis.gui

    app = qgis.core.QgsApplication()
    canvas = qgis.gui.QgsMapCanvas()
