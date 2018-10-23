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

imageFile = r"D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\058217622010_01\PCI Output\ATCOR\Basic\ATCORCorrected_o17OCT01084657_R1C12-058217622010_01_P001_14368025_PanSharp.tif"
imageFile = r"D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\058217622010_01\PCI Output\ATCOR\Basic\ATCORCorrected_o17OCT01084657-M2AS_R1C12-058217622010_01_P001_14368025.pix"
imageFile = r"D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\058217622010_01\PCI Output\ATCOR\SRTM\ATCORCorrected_o17OCT01084657-M2AS_R1C12-058217622010_01_P001_14368011.pix"
imageFile = r"D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\058217622010_01\PCI Output\ATCOR\SRTM+BRDF1\ATCORCorrected_o17OCT01084657-M2AS_R1C12-058217622010_01_P001_14344011.pix"
imageFile = r"D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\058217622010_01\PCI Output\Ortho\o17OCT01084657-M2AS_R1C12-058217622010_01_P001.TIF"
imageFile = r"D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\058217622010_01\PCI Output\ATCOR\SRTM+AdjCorr\ATCORCorrected_o17OCT01084657-M2AS_R1C12-058217622010_01_P001_14368043.pix"

demFile = r"V:\Data\NGI\GEF DEM\3323d_2015_1001_GEF_DEM_SGM3.tif"
slopeFile = r"V:\Data\NGI\GEF DEM\3323d_2015_1001_GEF_DEM_SGM3_slope.pix"
heightFile = r"V:\Data\NGI\GEF DEM\3323d_2015_1001_GEF_DEM_Photoscan_clip_hgt2.tif"
# imageFile = r"D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\058217622010_01\PCI Output\ATCOR\SRTM+AdjCorr\ATCORCorrected_o17OCT01084657_R1C12-058217622010_01_P001_14368043_PanSharpen.pix"
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
        imbuf_mask[:, i] = np.float64(band[mask]) / 5000.  # 5000 is scale for MODIS / XCALIB
    # imbuf_mask[:, 3] = imbuf_mask[:,  3]/2.
    # wv3 bands
    if np.any(imbuf_mask<0):
        print 'imbuf_mask < 0'
        # print np.where(imbuf_mask<0)
    s = np.sum(imbuf_mask[:,[1,2,4,6]], 1)   # NNB only sum r,g,b as ir confuses things in g_n
    #  s = np.sum(imbuf_mask[:,:4], 1)   # ??? check this
    cn = imbuf_mask / np.tile(s[:, None], (1, imbuf_mask.shape[1]))
    # wv 3 channels
    b_i = 1
    g_i = 2
    r_i = 4
    ir_i = 6
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


def ExtractDemPatchFeatures(imbuf, mask):
    mask = np.bool8(mask)
    imbuf_mask = np.ndarray(shape=(np.int32(mask.sum()), imbuf.shape[2]), dtype=np.float64)
    for i in range(0, imbuf.shape[2]):
        band = imbuf[:, :, i]
        imbuf_mask[:, i] = np.float64(band[mask])   # / 5000.  # 5000 is scale for MODIS / XCALIB
    # imbuf_mask[:, 3] = imbuf_mask[:,  3]/2.
    # wv3 bands
    if np.any(imbuf_mask<0):
        print 'imbuf_mask < 0'
        # print np.where(imbuf_mask<0)
    feat = {}

    feat['mean'] = imbuf_mask.mean()
    feat['std'] = imbuf_mask.std()
    feat['mean-min'] = (imbuf_mask-imbuf_mask.min()).mean()     #assuming flat which is wrong

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
    max_im_vals = np.zeros((ds.RasterCount))
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
                and np.all(plotCnrsPixel[:,1] < ds.RasterYSize) and plot.has_key('Yc') and plot['Yc'] > 0.:

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

            # adjust yc if it exists
            if plot.has_key('Yc') and plot['Yc']>0:
                plot['YcPp'] = plot['Yc'] / plotMask.sum()  # the average per pixel in the mask
                plot['YcPm2'] = plot['Yc'] / (plot['Size']**2)  # the average per m2 in the theoretical plot size
            else:
                print '%s - no yc' % (plot['ID'])
                # plot['YcPp'] = 0.
                # plot['YcPm2'] = 0.

            # extract image patch with mask
            imbuf = np.zeros((plotSizePixel[1], plotSizePixel[0], ds.RasterCount), dtype=float)
            for b in range(1, ds.RasterCount+1):
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
            tmp = np.reshape(feat['thumbnail'], (np.prod(plotSizePixel), ds.RasterCount))
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
        for b in range(0, ds.RasterCount):
            thumb[:, :, b] = thumb[:, :, b] / max_im_vals[b]
            thumb[:, :, b][thumb[:, :, b] > 1.] = 1.
        # thumb[:, :, 0] = thumb[:, :, 0] / 1.5
        # thumb[:, :, 1] = thumb[:, :, 1] * 1.2
        # thumb[:, :, 1][thumb[:, :, 1] > 1.] = 1.
        plot_dict[k]['thumbnail'] = thumb[:,:,[4,2,1]]

    return plot_dict


def ExtractAllDemFeatures(ds, csGtSpatialRef, csGtDict, plotFigures=False):  # , axis_signs=[1, 1]):
    # Note         A way around all this confusion may to be to layout plot co-ords in world co-ords and then convert to
    #              to pixel co-ords using world2pixel.  this would avoid all the fiddling and confusion in pixel space

    geotransform = ds.GetGeoTransform()
    transform = osr.CoordinateTransformation(csGtSpatialRef, osr.SpatialReference(ds.GetProjection()))
    i = 0
    plot_dict = {}
    # plotTagcDict = {}
    # class_labels = ['Pristine', 'Moderate', 'Severe']
    max_im_vals = np.zeros((ds.RasterCount))
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
                and np.all(plotCnrsPixel[:,1] < ds.RasterYSize) and plot.has_key('Yc') and plot['Yc'] > 0.:

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

            # adjust yc if it exists
            if plot.has_key('Yc') and plot['Yc']>0:
                plot['YcPp'] = plot['Yc'] / plotMask.sum()  # the average per pixel in the mask
                plot['YcPm2'] = plot['Yc'] / (plot['Size']**2)  # the average per m2 in the theoretical plot size
            else:
                print '%s - no yc' % (plot['ID'])
                # plot['YcPp'] = 0.
                # plot['YcPm2'] = 0.

            # extract image patch with mask
            imbuf = np.zeros((plotSizePixel[1], plotSizePixel[0], ds.RasterCount), dtype=float)
            for b in range(1, ds.RasterCount+1):
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
            feat = ExtractDemPatchFeatures(imbuf.copy(), plotMask)

            fields = plot.keys()
            for f in fields:
                feat[f] = csGtDict[plot['ID']][f]
            feat['thumbnail'] = np.float32(imbuf.copy())
            plot_dict[plot['ID']] = feat
            # plotTagcDict[plot['PLOT']] = csGtDict[plot['PLOT']]['TAGC']
            tmp = np.reshape(feat['thumbnail'], (np.prod(plotSizePixel), ds.RasterCount))
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
        # max_im_vals[1] = max_im_vals[1]
        for b in range(0, ds.RasterCount):
            thumb[:, :, b] = thumb[:, :, b] / max_im_vals[b]
            thumb[:, :, b][thumb[:, :, b] > 1.] = 1.
        # thumb[:, :, 0] = thumb[:, :, 0] / 1.5
        # thumb[:, :, 1] = thumb[:, :, 1] * 1.2
        # thumb[:, :, 1][thumb[:, :, 1] > 1.] = 1.
        plot_dict[k]['thumbnail'] = thumb.squeeze()

    return plot_dict

def ScatterD(x, y, labels=None, class_labels=None, thumbnails=None, regress=True, xlabel=None, ylabel=None):
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
                           fontdict={'size': 9, 'color': colour, 'weight': 'bold'})

            if thumbnails is not None:
                imbuf = np.array(thumbnails)[class_idx][xyi]
                ims = 20.
                extent = [xx - xd / (2 * ims), xx + xd / (2 * ims), yy - yd / (2 * ims), yy + yd / (2 * ims)]
                #pylab.imshow(imbuf[:, :, :3], extent=extent, aspect='auto')  # zorder=-1,
                pylab.imshow(imbuf, extent=extent, aspect='auto')  # zorder=-1,
                handles[ci] = ax.add_patch(
                    patches.Rectangle((xx - xd / (2 * ims), yy - yd / (2 * ims)), xd / ims, yd / ims, fill=False,
                                      edgecolor=colour, linewidth=2.))
                # pylab.plot(mPixels[::step], dRawPixels[::step], color='k', marker='.', linestyle='', markersize=.5)
        if regress and classes.__len__() > 1 and False:
            (slope, intercept, r, p, stde) = stats.linregress(x[class_idx], y[class_idx])
            pylab.text(xlim[0] + xd * 0.7, ylim[0] + yd * 0.05 * (ci + 2),
                       str.format('{1}: $R^2$ = {0:.2f}', np.round(r ** 2, 2), classes[ci]),
                       fontdict={'size': 10, 'color': colour})

    if regress:
        (slope, intercept, r, p, stde) = stats.linregress(x, y)
        pylab.text((xlim[0] + xd * 0.7), (ylim[0] + yd * 0.05), str.format('$R^2$ = {0:.2f}', np.round(r ** 2, 2)),
                   fontdict={'size': 12})
        print str.format('$R^2$ = {0:.2f}', np.round(r ** 2, 2))
        print str.format('P (slope=0) = {0:.2f}', np.round(p, 3))
        print str.format('Slope = {0:.2f}', np.round(slope, 3))
        print str.format('Std error of slope = {0:.2f}', np.round(stde, 3))
        yhat = x*slope + intercept
        print str.format('RMS error = {0:.2f}', np.round(np.sqrt(np.mean((y-yhat)**2)), 3))


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

    id = f['ID']
    if id[0] == 'S' or id[:3] == 'TCH':
        f['DegrClass'] = 'Severe'
    elif id[0] == 'M':
        f['DegrClass'] = 'Moderate'
    elif id[0] == 'P' or id[:3] == 'INT':
        f['DegrClass'] = 'Pristine'
    else:
        f['DegrClass'] = '?'

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
imageDs = gdal.OpenEx(imageFile, gdal.OF_RASTER)
if imageDs is None:
    print "Open failed./n"

print 'Driver: ', imageDs.GetDriver().ShortName, '/', \
    imageDs.GetDriver().LongName
print 'Size is ', imageDs.RasterXSize, 'x', imageDs.RasterYSize, \
    'x', imageDs.RasterCount
print 'Projection is ', imageDs.GetProjection()
geotransform = imageDs.GetGeoTransform()
if not geotransform is None:
    print 'Origin = (', geotransform[0], ',', geotransform[3], ')'
    print 'Pixel Size = (', geotransform[1], ',', geotransform[5], ')'
    pixelSize = geotransform[1]


plotDict = ExtractAllFeatures(imageDs, samplingPlotSpatialRef, samplingPlotGtDict, plotFigures=True)
imageDs = None


# plotDict = plotDict.copy()
# plotDict.update(vdwPlotDict.copy())


gn = np.array([plot['g_n'] for plot in plotDict.values()])
ndvi = np.array([plot['NDVI'] for plot in plotDict.values()])
rn = np.array([plot['r_n'] for plot in plotDict.values()])
bn = np.array([plot['b_n'] for plot in plotDict.values()])
ir_rat = np.array([plot['ir_rat'] for plot in plotDict.values()])
ir_n = np.array([plot['ir_n'] for plot in plotDict.values()])
id = np.array([plot['ID'] for plot in plotDict.values()])
ycpp = np.array([plot['YcPp'] for plot in plotDict.values()])*(100.**2)/(pixelSize**2)
yc = np.array([plot['YcHa'] for plot in plotDict.values()])
classes = np.array([plot['DegrClass'] for plot in plotDict.values()])

abg = np.array([plot['AgbHa'] for plot in plotDict.values()])
litter = np.array([plot['LitterHa'] for plot in plotDict.values()])

thumbnails = [plot['thumbnail'] for plot in plotDict.values()]
# classes = np.array([p['DegrClass'] for p in plotDict.values()])
# class_labels = ['Severe', 'Moderate', 'Pristine'] #np.unique(classes)

fig = pylab.figure()
ScatterD(np.log10(ndvi), yc/1000., class_labels=classes, labels=None, thumbnails=thumbnails, regress=True, xlabel='log(NDVI)', ylabel='AGB (t/ha)')
pylab.grid()
# fig.tight_layout()

fig = pylab.figure()
ScatterD(np.log10(gn), yc/1000., class_labels=classes, labels=None, thumbnails=thumbnails, regress=True, xlabel='log(gN)', ylabel='AGB (t/ha)')

fig = pylab.figure()
ScatterD(np.log10(rn), yc/1000., class_labels=classes, labels=id, thumbnails=thumbnails, regress=True, xlabel='log(rN)', ylabel='AGB (t/ha)')
pylab.grid()

fig = pylab.figure()
ScatterD(np.log10(rn), yc/1000., class_labels=classes, labels=id, thumbnails=thumbnails, regress=True, xlabel='log(rN)', ylabel='AGB (t/ha)')
pylab.grid()

fig = pylab.figure()
ScatterD(np.log10(rn), abg/1000., class_labels=classes, labels=id, thumbnails=thumbnails, regress=True, xlabel='log(rN)', ylabel='AGB (t/ha)')
pylab.grid()

fig = pylab.figure()
ScatterD(np.sqrt(ndvi), abg/1000., class_labels=classes, labels=id, thumbnails=thumbnails, regress=True, xlabel='log(NDVI)', ylabel='AGB (t/ha)')
pylab.grid()

fig = pylab.figure()
ScatterD(rn**2, abg/1000., class_labels=classes, labels=None, thumbnails=thumbnails, regress=True, xlabel='log(rN)', ylabel='AGB (t/ha)')
pylab.grid()


####################################################################################################################
# Read in the slope/dem
if False:
    heightDs = gdal.OpenEx(heightFile, gdal.OF_RASTER)
    if heightDs is None:
        print "Open failed./n"

    print 'Driver: ', heightDs.GetDriver().ShortName, '/', \
        heightDs.GetDriver().LongName
    print 'Size is ', heightDs.RasterXSize, 'x', heightDs.RasterYSize, \
        'x', heightDs.RasterCount
    print 'Projection is ', heightDs.GetProjection()
    geotransform = heightDs.GetGeoTransform()
    if not geotransform is None:
        print 'Origin = (', geotransform[0], ',', geotransform[3], ')'
        print 'Pixel Size = (', geotransform[1], ',', geotransform[5], ')'
        pixelSize = geotransform[1]

    heightPlotDict = ExtractAllDemFeatures(heightDs, samplingPlotSpatialRef, samplingPlotGtDict, plotFigures=True)
    heightDs = None


    mean = np.array([plot['mean'] for plot in heightPlotDict.values()])
    std = np.array([plot['std'] for plot in heightPlotDict.values()])
    hgt = np.array([plot['mean-min'] for plot in heightPlotDict.values()])
    ycpp = np.array([plot['YcPp'] for plot in heightPlotDict.values()])*(100.**2)/(pixelSize**2)
    yc = np.array([plot['YcHa'] for plot in heightPlotDict.values()])
    classes = np.array([plot['DegrClass'] for plot in heightPlotDict.values()])
    abg = np.array([plot['AgbHa'] for plot in heightPlotDict.values()])
    litter = np.array([plot['LitterHa'] for plot in heightPlotDict.values()])
    thumbnails = np.array([plot['thumbnail'] for plot in heightPlotDict.values()])

    pylab.close('all')
    fig = pylab.figure()
    ScatterD(hgt, yc/1000., class_labels=classes, labels=None, thumbnails=thumbnails, regress=True, xlabel='log(mean(height))', ylabel='AGB (t/ha)')
    pylab.grid()

    # try only severe/mod class - this area is easier to find ground because it is flatter and has less vegetation
    fig = pylab.figure()
    subIdx = (classes == 'Severe') | (classes == 'Moderate')
    ScatterD(mean[subIdx], yc[subIdx]/1000., class_labels=classes[subIdx], labels=None, thumbnails=thumbnails[subIdx], regress=True, xlabel='log(mean(height))', ylabel='AGB (t/ha)')
    pylab.grid()


####################################################################################################################

# lock at ycha KDS density
from scipy.stats import gaussian_kde
class_labels = np.unique(classes)
ycGrid = np.linspace(yc.min(), yc.max(), 50)

# class_num = [30, 30, 30]
# featureValSub = np.array([])
fontSize=12
pylab.figure()
for i, cl in enumerate(class_labels):
    idx = classes == cl
    #idx = idx[:class_num[i]]
    #allIdx.append(idx)
    classYc = yc[idx]
    kde = gaussian_kde(classYc)  # , bw_method=bandwidth / height.std(ddof=1))
    # ndviGrid = np.linspace(ndvi.min(), ndvi.max(), 100)
    ycKde = kde.evaluate(ycGrid)

    pylab.subplot(2, 2, i+1)
    pylab.plot(ycGrid, ycKde)
    pylab.xlabel('Yc (t/ha)', fontdict={'size':fontSize})
    pylab.ylabel('Prob. Density', fontdict={'size':fontSize})
    pylab.title(cl, fontdict={'size':fontSize})
    pylab.grid()
    pylab.axis('tight')
    pylab.text(0.9*ycGrid.max(), 0.9*ycGrid.max(),'N=%d'%(idx.sum()), fontdict={'size':fontSize})
    print '%s N: %d' % (cl, idx.sum())

kde = gaussian_kde(np.array(yc))  # , bw_method=bandwidth / height.std(ddof=1))
ycKde = kde.evaluate(ycGrid)

pylab.subplot(2, 2, 4)
pylab.plot(ycGrid, ycKde)
pylab.xlabel('Yc (t/ha)', fontdict={'size': fontSize})
pylab.ylabel('Prob. Density', fontdict={'size': fontSize})
pylab.title('All', fontdict={'size': fontSize})
pylab.grid()
pylab.axis('tight')
pylab.text(0.9 * ycGrid.max(), 0.9 * ycGrid.max(), 'N=%d' % (idx.sum()), fontdict={'size': fontSize})

pylab.tight_layout()

#--------------------------------------------------------------------------------------------------------------
# simulate more plots to see spread of yc

ycBs = np.array([])
nClasses = [20, 30, 30]    #['Moderate', 'Pristine', 'Severe']
pylab.figure()
for i, cl in enumerate(class_labels):
    idx = classes == cl
    #idx = idx[:class_num[i]]
    #allIdx.append(idx)
    classYc = yc[idx]
    classYcBs = np.random.choice(classYc, nClasses[i], replace=True)   #generate a bootstrap sample N[i]
    ycBs = np.concatenate((ycBs, classYcBs))

    kde = gaussian_kde(classYcBs)  # , bw_method=bandwidth / height.std(ddof=1))
    # ndviGrid = np.linspace(ndvi.min(), ndvi.max(), 100)
    ycKde = kde.evaluate(ycGrid)

    pylab.subplot(2, 2, i+1)
    pylab.plot(ycGrid, ycKde)
    pylab.xlabel('Yc (t/ha)', fontdict={'size':fontSize})
    pylab.ylabel('Prob. Density', fontdict={'size':fontSize})
    pylab.title(cl, fontdict={'size':fontSize})
    pylab.grid()
    pylab.axis('tight')
    pylab.text(0.9*ycGrid.max(), 0.9*ycGrid.max(),'N=%d'%(idx.sum()), fontdict={'size':fontSize})
    print '%s N: %d' % (cl, idx.sum())

kde = gaussian_kde(np.array(ycBs))  # , bw_method=bandwidth / height.std(ddof=1))
ycKde = kde.evaluate(ycGrid)

pylab.subplot(2, 2, 4)
pylab.plot(ycGrid, ycKde)
pylab.xlabel('Yc (t/ha)', fontdict={'size': fontSize})
pylab.ylabel('Prob. Density', fontdict={'size': fontSize})
pylab.title('All', fontdict={'size': fontSize})
pylab.grid()
pylab.axis('tight')
pylab.text(0.9 * ycGrid.max(), 0.9 * ycGrid.max(), 'N=%d' % (idx.sum()), fontdict={'size': fontSize})




##############################################################################################################
# play with regression analysis

import statsmodels.api as sm

x = sm.add_constant(np.log10(ndvi))
model = sm.OLS(yc,x).fit()
predictions = model.predict(x)
rms = np.sqrt(((yc-predictions)**2).mean())
print 'rms=%.3f'%(rms)
model.summary()

x = sm.add_constant(np.array(np.log10([ndvi, gn])).transpose())
model = sm.OLS(yc,x).fit_regularized(method='elastic_net', alpha=5, L1_wt=1.0)
predictions = model.predict(x)
rms = np.sqrt(((yc-predictions)**2).mean())
print 'rms=%.3f'%(rms)
model.summary()
print (model.params)


x = np.log10(ndvi)
model = sm.OLS(yc, x).fit()
predictions = model.predict(x)
rms = np.sqrt(((yc-predictions)**2).mean())
print 'rms=%.3f'%(rms)
# Print out the statistics
model.summary()


(slope, intercept, r, p, stde) = stats.linregress(ndvi, yc)
##############################################################################################################
# play with sklearn

from sklearn import datasets
from sklearn.model_selection import cross_val_predict, cross_validate, GridSearchCV
from sklearn import linear_model
from sklearn import metrics
from sklearn.feature_selection import f_regression, mutual_info_regression, RFECV, RFE
from sklearn.svm import SVR, NuSVR
import matplotlib.pyplot as plt

# make an np.array of features
featKeysOrig = ['i', 'r_n', 'g_n', 'b_n', 'ir_n', 'ir_rat', 'SAVI', 'NDVI', 'i_std', 'NDVI_std']
# featKeysOrig = ['i', 'r_n', 'g_n', 'b_n', 'ir_n', 'SAVI', 'NDVI', 'i_std', 'NDVI_std']   # hack for wv3
X = []
featKeys = []
for featKey in featKeysOrig:
    f = np.array([plot[featKey] for plot in plotDict.values()])
    X.append(f)
    featKeys.append(featKey)
for featKey in featKeysOrig:
    f = np.array([plot[featKey] for plot in plotDict.values()])
    X.append(np.log10(f))
    featKeys.append(str.format('log10({0})', featKey))
for featKey in featKeysOrig:
    f = np.array([plot[featKey] for plot in plotDict.values()])
    X.append(f**2)
    featKeys.append(str.format('{0}^2', featKey))
for featKey in featKeysOrig:
    f = np.array([plot[featKey] for plot in plotDict.values()])
    X.append(np.sqrt(f))
    featKeys.append(str.format('{0}^.5', featKey))

id = np.array([plot['ID'] for plot in plotDict.values()])
featKeys = np.array(featKeys)
X = np.array(X).transpose()
np.isnan(X).sum()
X[np.isnan(X)] = -1000  # NB HACK wv3 weirdness for now

# plot['Yc'] is the mean Yc per pixel in the plot
if False:
    y = np.array([plot['YcPp'] for plot in plotDict.values()])*(100.**2)/(pixelSize**2)
    y = np.array([plot['YcHa'] for plot in plotDict.values()])
else:
    y = np.array([plot['AgbHa'] for plot in plotDict.values()])
    litterHa = np.array([plot['LitterHa'] for plot in plotDict.values()])
    agbIdx = litterHa>=0.
    # exclude plots with no litter data
    # y = y[agbIdx]
    # X = X[agbIdx, :]

f_test, _ = f_regression(X, y)
f_test /= np.max(f_test)
print 'Best F test feature: ' + featKeys[np.argmax(f_test)]
idx = np.argsort(-f_test)
print 'Features ranked by F test: %s'%(featKeys[idx])

mi = mutual_info_regression(X, y)
mi /= np.max(mi)
print 'Best MI feature: ' + featKeys[np.argmax(mi)]
idx = np.argsort(-mi)
print 'Features ranked by MI: %s' % (featKeys[idx])

rfe = RFE(linear_model.LinearRegression(), 2, step=1)
rfe.fit(X, y)
print 'Features ranked by RFE: %s' % (featKeys[rfe.ranking_.argsort()])
print 'Model selected by RFE: %s' % (featKeys[rfe.support_])

lasso = linear_model.LassoCV()
lasso.fit(X,y/1000.)
print 'Best Lasso feature: ' + featKeys[np.argmax(np.abs(lasso.coef_))]
print 'Main Lasso features: %s' % (featKeys[np.abs(lasso.coef_)>5])
print 'Features ranked by Lasso: %s' % (featKeys[np.argsort(-np.abs(lasso.coef_))])

lasso.coef_[np.argsort(-np.abs(lasso.coef_))]
np.argsort(-np.abs(lasso.coef_))

from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2

# SelectKBest
skb = SelectKBest(mutual_info_regression, k=2)   #this is just a ranking
skb.fit(X, y)
predicted = cross_val_predict(linear_model.LinearRegression(), X[:, skb.get_support()], y, cv=y.__len__()-1)
print 'Main Lasso SelectKBest: %s' % (featKeys[skb.get_support()])
print 'Features ranked by SelectKBest: %s' % (featKeys[np.argsort(-np.abs(skb.scores_))])


# KNN
from sklearn import neighbors
from sklearn import preprocessing
knn = neighbors.KNeighborsRegressor(3)
predicted = cross_val_predict(knn, preprocessing.scale(X[:, [11, 17]]), y, cv=y.__len__()-1)


# SVR
svr = GridSearchCV(SVR(kernel='rbf', gamma=0.1), cv=5,
                   param_grid={"C": [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8],
                               "gamma": np.logspace(-2, 4, 5)})
svr.fit(X[:, [11, 17]], y)
svr.best_params_

# PCA
from sklearn.decomposition import PCA
pca = PCA(n_components=.9, whiten=False)
pca.fit(X, y)
Xpca = pca.transform(X)
predicted = cross_val_predict(linear_model.LinearRegression(), Xpca, y, cv=y.__len__()-1)
rms = np.sqrt(((y-predicted)**2).mean())
r2 = metrics.r2_score(yc, predicted)
print 'R2 = %.3f'%r2
print 'rms = %.3f'%rms


predicted = cross_val_predict(linear_model.LinearRegression(), X[:, [19, 25]], y, cv=y.__len__()-1)
predicted = cross_val_predict(linear_model.LinearRegression(), X[:, [24, 19, 10, 25]], y, cv=y.__len__()-1)
predicted = cross_val_predict(linear_model.LinearRegression(), X[:, [1,11]], y, cv=y.__len__()-1)
predicted = cross_val_predict(linear_model.LassoCV(), X, y, cv=y.__len__() - 1)
predicted = cross_val_predict(svr.best_estimator_, X[:, [11, 17]], y, cv=y.__len__() - 1)
rms = np.sqrt(((y-predicted)**2).mean())
r2 = metrics.r2_score(yc, predicted)
print 'R2 = %.3f'%r2
print 'rms = %.3f'%rms


fig, ax = plt.subplots()
ax.scatter(y, predicted, edgecolors=(0, 0, 0))
ax.plot([y.min(), y.max()], [y.min(), y.max()], 'k--', lw=4)
ax.set_xlabel('Measured')
ax.set_ylabel('Predicted')
plt.show()



############################################################
# results for report
# cross validation for score variance
# feats = [19, 10, 25]
feats = [24, 19, 10, 25]
feats = [10, 24, 18, 25]  # for june data from Lasso

feats = [17, 11]  # wv ms (atcor basic) abg
feats = [7, 11, 17, 24]  # wv ms (atcor srtm+adjcorr) y=abg
# Features: ['NDVI' 'log10(r_n)' 'log10(NDVI)' 'ir_n^2']
# R2: 0.8298
# RMSE: 7.89
# Method 1: RMSE: 7.956, 5-95% CI: 0.297 - 16.288
# Method 2: RMSE: 6.129, 5-95% CI: 0.296 - 16.288


feats = [11, 17,  5,  7]  # wv ms nir=b6 (atcor srtm+adjcorr) y=abg
# R2: 0.8333
# RMSE: 7.81
# Method 1: RMSE: 7.870, 5-95% CI: 0.633 - 15.172
# Method 2: RMSE: 6.262, 5-95% CI: 0.625 - 15.172


feats = [11,  5, 17, 25]  # wv ms nir=b6, s=sum(r,g,b,nir) (atcor srtm+adjcorr) y=abg
# Features: ['log10(r_n)' 'ir_rat' 'log10(NDVI)' 'ir_rat^2']
# R2: 0.8419
# RMSE: 7.61
# Method 1: RMSE: 7.666, 5-95% CI: 0.994 - 15.103
# Method 2: RMSE: 6.175, 5-95% CI: 0.993 - 15.103

feats = [11,  5, 17, 25]  # wv ms nir=b6, s=sum(r,g,b,nir) (ortho) y=abg
# Features: ['log10(r_n)' 'ir_rat' 'log10(NDVI)' 'ir_rat^2']
# R2: 0.8413
# RMSE: 7.62
# Method 1: RMSE: 7.662, 5-95% CI: 0.746 - 14.831
# Method 2: RMSE: 6.177, 5-95% CI: 0.743 - 14.829

feats = [11,  9, 17,  5]  # wv ms nir=b6, s=sum(r,g,b,nir) (atcor srtm+adjcorr) y = ycpp
# Features: ['log10(r_n)' 'ir_rat' 'log10(NDVI)' 'ir_rat^2']
# R2: 0.7566
# RMSE: 7.45
# Method 1: RMSE: 7.512, 5-95% CI: 0.511 - 14.360
# Method 2: RMSE: 6.084, 5-95% CI: 0.508 - 14.360

#[11, 17,  5,  9, 10, 18,
feats = [11, 17,  5]  # wv ms nir=b6, s=sum(r,g,b,nir) (atcor srtm+adjcorr) y = ycha
# Features: ['log10(r_n)' 'log10(NDVI)' 'ir_rat']
# R2: 0.7941
# RMSE: 6.98
# Method 1: RMSE: 7.035, 5-95% CI: 0.515 - 13.139
# Method 2: RMSE: 5.724, 5-95% CI: 0.515 - 13.139


feats = [11, 17,  9,  5] # wv ms nir=b6, s=sum(r,g,b,nir) (atcor srtm+adjcorr) y = agb, incl trial litter data
# Features: ['log10(r_n)' 'log10(NDVI)' 'NDVI_std' 'ir_rat']
# R2: 0.8469
# RMSE: 7.22
# Method 1: RMSE: 7.273, 5-95% CI: 0.809 - 13.621
# Method 2: RMSE: 6.034, 5-95% CI: 0.809 - 13.613


# feats = [11, 17, 25, 9]

# 10, 24, 18, 25,

# feats = [11, 17]
# feats = [11]
# feats = np.arange(0, X.shape[1])
yt = y/1000.
# est = linear_model.RidgeCV()
# est = SVR(kernel='rbf', C=1e6, gamma=0.01)
est = linear_model.LinearRegression()
scores = cross_validate(est, X[:, feats], yt, scoring=('r2', 'neg_mean_squared_error'), cv=yt.__len__()-1)
predicted = cross_val_predict(est, X[:, feats], yt, cv=yt.__len__()-1)
r2 = metrics.r2_score(yt, predicted)
rmse = np.sqrt(metrics.mean_squared_error(yt, predicted))
print 'Features: %s' % (featKeys[feats])
print 'R2: {0:.4f}'.format(r2)
print 'RMSE: {0:.2f}'.format(rmse)
mse = (-scores['test_neg_mean_squared_error'])
print 'Method 1: RMSE: {0:.3f}, 5-95% CI: {1:.3f} - {2:.3f}'.format(np.sqrt(mse.mean()), np.sqrt(np.percentile(mse,5)), np.sqrt(np.percentile(mse,95)))
rms = np.sqrt(-scores['test_neg_mean_squared_error'])
print 'Method 2: RMSE: {0:.3f}, 5-95% CI: {1:.3f} - {2:.3f}'.format(rms.mean(), np.percentile(rms,5), np.percentile(rms,95))


fig, ax = plt.subplots()
ax.scatter(yt, predicted, edgecolors=(0, 0, 0))
ax.plot([yt.min(), yt.max()], [yt.min(), yt.max()], 'k--', lw=4)
ax.set_xlabel('Measured AGB (t/ha)')
ax.set_ylabel('Predicted AGB (t/ha)')
plt.xlim(0, yt.max())
plt.ylim(0, yt.max())
plt.grid()
plt.text(42, 15, str.format('$R^2$ = {0:.2f}', r2),
           fontdict={'size': 11})
plt.text(42, 11, str.format('RMSE = {0:.2f} t/ha', rmse),
           fontdict={'size': 11})
plt.show()


fig = pylab.figure()
ScatterD(rn, yc/1000., class_labels=classes, labels=None, thumbnails=thumbnails, regress=False, xlabel='$R/(R+G+B)$', ylabel='Woody C (t/ha)')
pylab.grid()


#####################################################
lr = linear_model.LinearRegression()
y = yc

# cross_val_predict returns an array of the same size as `y` where each entry
# is a prediction obtained by cross validation:
x = sm.add_constant(np.log10(ndvi))
# x = np.column_stack(([ndvi]))
predicted = cross_val_predict(lr, x, yc, cv=10)

fig, ax = plt.subplots()
ax.scatter(y, predicted, edgecolors=(0, 0, 0))
ax.plot([y.min(), y.max()], [y.min(), y.max()], 'k--', lw=4)
ax.set_xlabel('Measured')
ax.set_ylabel('Predicted')

plt.show()

rms = np.sqrt(((y-predicted)**2).mean())
r2 = metrics.r2_score(yc, predicted)
print '$R^2$ = %.3f'%r2
print 'rms = %.3f'%rms
