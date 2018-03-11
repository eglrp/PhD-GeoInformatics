# look at the relationship between NDVI and actual yc from GEF sampling data to date


import gdal
import ogr
import numpy as np
import osr
import pylab
from scipy import stats as stats
from matplotlib import patches
import matplotlib.pyplot as plt
from scipy import ndimage as ndimage
from scipy.stats import gaussian_kde

# Python Imaging Library imports
from PIL import Image
from PIL import ImageDraw

tchnuganuImageFile = "V:/Data/NGI/Rectified/3323D_2015_1001/RGBN/XCALIB/o3323d_2015_1001_02_0081_RGBN_XCALIB.tif"  # "V:/Data/NGI/Rectified/3324C_2015_1004/RGBN/o3324c_2015_1004_02_0044_RGBN_XCALIB.tif"
vdwImageFile = "V:/Data/NGI/Rectified/3323D_2015_1001/RGBN/XCALIB/o3323d_2015_1001_02_0078_RGBN_XCALIB.tif"  # ""V:/Data/NGI/Rectified/3323D_2015_1001/RGBN/o3323d_2015_1001_02_0077_Lo25Wgs84_RGBN_XCALIB.tif"
# samplingPtFile = "C:/Data/Development/Projects/PhD GeoInformatics/Data/GEF Sampling/GEF Sampling Points.shp"
samplingPlotGtFile = "C:/Data/Development/Projects/PhD GeoInformatics/Data/GEF Sampling/GEF Plot Polygons with Yc.shp"


pylab.close('all')

def World2Pixel(geoMatrix, x, y):
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



def ExtractPatchFeatures(imbuf, mask):
    mask = np.bool8(mask)
    imbuf_mask = np.ndarray(shape=(np.int32(mask.sum()), imbuf.shape[2]), dtype=np.float64)
    for i in range(0, imbuf.shape[2]):
        band = imbuf[:, :, i]
        imbuf_mask[:, i] = np.float64(band[mask]) / 100.
    # imbuf_mask[:, 3] = imbuf_mask[:,  3]/2.
    s = np.sum(imbuf_mask[:,:3], 1)   # NNB only sum r,g,b as ir confuses things in g_n
    cn = imbuf_mask / np.tile(s[:, None], (1, imbuf_mask.shape[1]))
    b_i = 2
    g_i = 1
    r_i = 0
    ir_i = 3
    ndvi = (imbuf_mask[:, ir_i] - imbuf_mask[:, r_i]) / (imbuf_mask[:, ir_i] + imbuf_mask[:, r_i])
    ir_rat = imbuf_mask[:, ir_i] / imbuf_mask[:, r_i]
    L = 0.05
    savi = (1 + L) * (imbuf_mask[:, ir_i] - imbuf_mask[:, r_i]) / (L + imbuf_mask[:, ir_i] + imbuf_mask[:, r_i])
    feat = {}
    feat['r_n'] = cn[:, r_i].mean()
    feat['g_n'] = cn[:, g_i].mean()
    feat['b_n'] = cn[:, b_i].mean()
    feat['ir_n'] = cn[:, ir_i].mean()
    feat['NDVI'] = ndvi.mean()
    feat['SAVI'] = savi.mean()
    feat['ir_rat'] = ir_rat.mean()
    feat['i'] = (s / imbuf_mask.shape[1]).mean()
    feat['i_std'] = (s / imbuf_mask.shape[1]).std()
    feat['NDVI_std'] = ndvi.std()

    return feat

def ExtractAllFeatures(ds, csGtSpatialRef, csGtDict, plotFigures=False):  # , axis_signs=[1, 1]):
    # Note         A way around all this confusion may to be to layout plot co-ords in world co-ords and then convert to
    #              to pixel co-ords using world2pixel.  this would avoid all the fiddling and confusion in pixel space

    geotransform = ds.GetGeoTransform()
    transform = osr.CoordinateTransformation(csGtSpatialRef, osr.SpatialReference(ds.GetProjection()))
    i = 0
    plot_dict = {}
    # plotTagcDict = {}
    # class_labels = ['Pristine', 'Moderate', 'Severe']
    max_im_vals = np.zeros((4))
    for plot in csGtDict.values():
        # transform plot corners into ds pixel space
        plotCnrsWorld = plot['points']
        plotCnrsPixel = []
        for cnr in plotCnrsWorld:
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(cnr[0], cnr[1])
            point.Transform(transform)  # xform into im projection
            (pixel, line) = World2Pixel(geotransform, point.GetX(), point.GetY())
            plotCnrsPixel.append((pixel,line))
        plotCnrsPixel = np.array(plotCnrsPixel)

        # not all the point fall inside the image
        if np.all(plotCnrsPixel >=0) and  np.all(plotCnrsPixel[:,0] < ds.RasterXSize) \
                and np.all(plotCnrsPixel[:,1] < ds.RasterYSize):

            # get winddow extents
            ulCnr = np.floor(np.min(plotCnrsPixel, 0))
            lrCnr = np.ceil(np.max(plotCnrsPixel, 0))
            plotSizePixel = np.int32(lrCnr - ulCnr)

            # make a mask for this plot
            img = Image.fromarray(np.zeros((plotSizePixel[1], plotSizePixel[0])))

            # Draw a rotated rectangle on the image.
            draw = ImageDraw.Draw(img)
            # rect = get_rect(x=120, y=80, width=100, height=40, angle=30.0)
            draw.polygon([tuple((p - ulCnr)) for p in plotCnrsPixel], fill=1)
            # Convert the Image data to a numpy array.
            plotMask = np.asarray(img)

            # extract image patch with mask
            imbuf = np.zeros((plotSizePixel[1], plotSizePixel[0], 4), dtype=float)
            for b in range(1, 5):
                imbuf[:, :, b - 1] = ds.GetRasterBand(b).ReadAsArray(ulCnr[0], ulCnr[1], plotSizePixel[0],
                                                                     plotSizePixel[1])

            # imbuf[:, :, 3] = imbuf[:, :, 3] / 2  # hack for NGI XCALIB
            if np.all(imbuf == 0):
                print plot['ID'] + ": imbuf zero, assume NODATA ommitting"
                break
            # for b in range(0, 4):
            #     imbuf[:, :, b] = imbuf[:, :, b] / max_im_vals_[b]
            if not plotMask.shape == imbuf.shape[0:2]:
                print "error - mask and buf different sizes"
                raise Exception("error - mask and buf different sizes")
            feat = ExtractPatchFeatures(imbuf.copy(), plotMask)

            fields = plot.keys()
            for f in fields:
                feat[f] = csGtDict[plot['ID']][f]
            feat['thumbnail'] = np.float32(imbuf.copy())
            plot_dict[plot['ID']] = feat
            # plotTagcDict[plot['PLOT']] = csGtDict[plot['PLOT']]['TAGC']
            tmp = np.reshape(feat['thumbnail'], (np.prod(plotSizePixel), 4))
            # max_tmp = tmp.max(axis=0)
            max_tmp = np.percentile(tmp, 98., axis=0)
            max_im_vals[max_tmp > max_im_vals] = max_tmp[max_tmp > max_im_vals]
            # print plot['PLOT']
            i = i + 1
        else:
            print "x-" + plot['ID']

    print i
    for k, v in plot_dict.iteritems():
        thumb = v['thumbnail']
        max_im_vals[1] = max_im_vals[1]
        for b in range(0, 4):
            thumb[:, :, b] = thumb[:, :, b] / max_im_vals[b]
            thumb[:, :, b][thumb[:, :, b] > 1.] = 1.
        # thumb[:, :, 0] = thumb[:, :, 0] / 1.5
        # thumb[:, :, 1] = thumb[:, :, 1] * 1.2
        # thumb[:, :, 1][thumb[:, :, 1] > 1.] = 1.
        plot_dict[k]['thumbnail'] = thumb

    return plot_dict


def scatterd(x, y, labels=None, class_labels=None, thumbnails=None, regress=True, xlabel=None, ylabel=None):
    if class_labels is None:
        class_labels = np.zeros(x.__len__())
    classes = np.unique(class_labels)
    colours = ['r', 'm', 'b', 'g', 'y', 'k', 'o']
    ylim = [np.min(y), np.max(y)]
    xlim = [np.min(x), np.max(x)]
    xd = np.diff(xlim)[0]
    yd = np.diff(ylim)[0]
    pylab.axis([np.min(x), np.max(x), np.min(y), np.max(y)])
    pylab.hold('on')
    ax = pylab.gca()
    handles = [0, 0, 0]

    for ci, (class_label, colour) in enumerate(zip(classes, colours[:classes.__len__()])):
        class_idx = class_labels == class_label
        if thumbnails is None:
            pylab.plot(x[class_idx], y[class_idx], colour + 'x')

        for xyi, (xx, yy) in enumerate(zip(x[class_idx], y[class_idx])):  # , np.array(plot_names)[class_idx]):
            if labels is not None:
                pylab.text(xx - .0015, yy - .0015, np.array(labels)[class_idx][xyi],
                           fontdict={'size': 8, 'color': colour})

            if thumbnails is not None:
                imbuf = np.array(thumbnails)[class_idx][xyi]
                ims = 20.
                extent = [xx - xd / (2 * ims), xx + xd / (2 * ims), yy - yd / (2 * ims), yy + yd / (2 * ims)]
                pylab.imshow(imbuf[:, :, :3], extent=extent, aspect='auto')  # zorder=-1,
                handles[ci] = ax.add_patch(
                    patches.Rectangle((xx - xd / (2 * ims), yy - yd / (2 * ims)), xd / ims, yd / ims, fill=False,
                                      edgecolor=colour, linewidth=2.))
                # pylab.plot(mPixels[::step], dRawPixels[::step], color='k', marker='.', linestyle='', markersize=.5)
        if regress and classes.__len__() > 1:  # and False:
            (slope, intercept, r, p, stde) = stats.linregress(x[class_idx], y[class_idx])
            pylab.text(xlim[0] + xd * 0.7, ylim[0] + yd * 0.05 * (ci + 2),
                       str.format('{1}: $R^2$ = {0:.2f}', np.round(r ** 2, 2), classes[ci]),
                       fontdict={'size': 10, 'color': colour})

    if regress:
        (slope, intercept, r, p, stde) = stats.linregress(x, y)
        pylab.text((xlim[0] + xd * 0.4), (ylim[0] + yd * 0.05), str.format('All: $R^2$ = {0:.2f}', np.round(r ** 2, 2)),
                   fontdict={'size': 12})

    if xlabel is not None:
        pylab.xlabel(xlabel, fontdict={'size': 12})
    if ylabel is not None:
        pylab.ylabel(ylabel, fontdict={'size': 12})
    # pylab.ylabel(yf)
    if classes.__len__() > 1:
        pylab.legend(handles, classes, fontsize=12)


###########################################################################################################

# read in sampling pts
ds = gdal.OpenEx(samplingPlotGtFile, gdal.OF_VECTOR)
if ds is None:
    print "Open failed./n"

lyr = ds.GetLayerByIndex(0)
lyr.ResetReading()
samplingPlotSpatialRef = lyr.GetSpatialRef()

# gcpList = []
samplingPlotGtDict = {}
for (i, feat) in enumerate(lyr):
    print '.',
    feat_defn = lyr.GetLayerDefn()
    f = {}
    for i in range(feat_defn.GetFieldCount()):
        field_defn = feat_defn.GetFieldDefn(i)
        f[field_defn.GetName()] = feat.GetField(i)
    geom = feat.GetGeometryRef()
    if geom is not None and (geom.GetGeometryType() == ogr.wkbPolygon):
        print "%s N Pts: %d" % (f['ID'], geom.GetGeometryRef(0).GetPointCount())
        f['geom'] = geom.Clone()
        f['points'] = geom.GetGeometryRef(0).GetPoints()[:-1]
        # pixCnr = []
        # for point in f['points']:
        #     pixCnr.append(World2Pixel(geotransform, point[0], point[1]))

    else:
        print "no polygon geometry/n"
    # gcpList.append(f)
    samplingPlotGtDict[f['ID']] = f
print ' '

ds = None


# Read in the images
tchnuganuDs = gdal.OpenEx(tchnuganuImageFile, gdal.OF_RASTER)
if tchnuganuDs is None:
    print "Open failed./n"

print 'Driver: ', tchnuganuDs.GetDriver().ShortName, '/', \
    tchnuganuDs.GetDriver().LongName
print 'Size is ', tchnuganuDs.RasterXSize, 'x', tchnuganuDs.RasterYSize, \
    'x', tchnuganuDs.RasterCount
print 'Projection is ', tchnuganuDs.GetProjection()
geotransform = tchnuganuDs.GetGeoTransform()
if not geotransform is None:
    print 'Origin = (', geotransform[0], ',', geotransform[3], ')'
    print 'Pixel Size = (', geotransform[1], ',', geotransform[5], ')'


tchnuganuPlotDict = ExtractAllFeatures(tchnuganuDs, samplingPlotSpatialRef, samplingPlotGtDict, plotFigures=True)
tchnuganuDs = None

vdwDs = gdal.OpenEx(vdwImageFile, gdal.OF_RASTER)
if vdwDs is None:
    print "Open failed./n"
geotransform = vdwDs.GetGeoTransform()
if not geotransform is None:
    print 'Origin = (', geotransform[0], ',', geotransform[3], ')'
    print 'Pixel Size = (', geotransform[1], ',', geotransform[5], ')'

vdwPlotDict = ExtractAllFeatures(vdwDs, samplingPlotSpatialRef, samplingPlotGtDict, plotFigures=True)
vdwDs = None

plotDict = tchnuganuPlotDict.copy()
plotDict.update(vdwPlotDict.copy())




featureName = 'NDVI'

featureVal = np.array([plot[featureName] for plot in plotDict.values()])
# featureVal = np.log10(np.array([plot[featureName] for plot in plotDict.values()]))

featureGrid = np.linspace(featureVal.min(), featureVal.max(), 100)

pylab.figure()

classes = np.array([p['DegrClass'] for p in plotDict.values()])
class_labels = ['Severe', 'Moderate', 'Pristine'] #np.unique(classes)
class_num = [35, 25, 30]
featureValSub = np.array([])
for i, cl in enumerate(class_labels):
    idx = classes == cl
    #idx = idx[:class_num[i]]
    #allIdx.append(idx)
    classFeatureVal = featureVal[idx][:class_num[i]]
    featureValSub = np.concatenate((featureValSub, classFeatureVal))
    kde = gaussian_kde(classFeatureVal)  # , bw_method=bandwidth / height.std(ddof=1))
    # ndviGrid = np.linspace(ndvi.min(), ndvi.max(), 100)
    featureKde = kde.evaluate(featureGrid)

    pylab.subplot(2, 2, i+1)
    pylab.plot(featureGrid, featureKde)
    pylab.xlabel(featureName)
    pylab.ylabel('Density')
    pylab.title(cl)
    pylab.grid()
    pylab.text(0.9*featureGrid.max(), 0.9*featureKde.max(),'N=%d'%(class_num[i]))

kde = gaussian_kde(np.array(featureValSub))  # , bw_method=bandwidth / height.std(ddof=1))
featureKde = kde.evaluate(featureGrid)

pylab.subplot(2, 2, 4)
pylab.plot(featureGrid, featureKde)
pylab.xlabel(featureName)
pylab.ylabel('Density')
pylab.title('All')
pylab.grid()


gn = np.array([plot['g_n'] for plot in plotDict.values()])
ndvi = np.array([plot['NDVI'] for plot in plotDict.values()])
id = np.array([plot['ID'] for plot in plotDict.values()])
thumbnails = [plot['thumbnail'] for plot in plotDict.values()]

pylab.figure()
scatterd(gn, ndvi, labels=id, class_labels=classes, thumbnails=thumbnails, regress=False, xlabel='gn', ylabel='NDVI')


# plotDict.values()[0].keys()


