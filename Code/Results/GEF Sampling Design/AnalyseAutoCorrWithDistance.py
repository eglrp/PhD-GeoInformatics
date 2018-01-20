# Look at how autocorr of a patch changes with distance for each stratum
# to get an idea of how far apart plots in a cluster should be and also just how far apart plots in general should be


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
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

# Python Imaging Library imports
from PIL import Image
from PIL import ImageDraw

tchnuganuImageFile = "V:/Data/NGI/Rectified/3323D_2015_1001/RGBN/XCALIB/o3323d_2015_1001_02_0081_RGBN_XCALIB.tif"  # "V:/Data/NGI/Rectified/3324C_2015_1004/RGBN/o3324c_2015_1004_02_0044_RGBN_XCALIB.tif"
vdwSewefonteinImageFile = "V:/Data/NGI/Rectified/3323D_2015_1001/RGBN/XCALIB/o3323d_2015_1001_02_0078_RGBN_XCALIB.tif"  # ""V:/Data/NGI/Rectified/3323D_2015_1001/RGBN/o3323d_2015_1001_02_0077_Lo25Wgs84_RGBN_XCALIB.tif"
samplingAreaFile = "C:/Data/Development/Projects/PhD GeoInformatics/Data/GEF Sampling/GEF Sampling Areas V2.shp"
samplingPtFile = "C:/Data/Development/Projects/PhD GeoInformatics/Data/GEF Sampling/GEF Sampling Points.shp"
pylab.close('all')

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


# def extract_patch_features(imbuf):
#     imbuf = np.float64(imbuf)  # values are in percent?
#     s = np.sum(imbuf, 2)
#     cn = imbuf / np.tile(s[:, :, None], (1, 1, imbuf.shape[2]))
#     b_i = 2
#     g_i = 1
#     r_i = 0
#     ir_i = 3
#     ndvi = (imbuf[:, :, ir_i] - imbuf[:, :, r_i]) / (imbuf[:, :, ir_i] + imbuf[:, :, r_i])
#     ir_rat = imbuf[:, :, ir_i] / imbuf[:, :, r_i]
#     L = 0.05
#     savi = (1 + L) * (imbuf[:, :, ir_i] - imbuf[:, :, r_i]) / (L + imbuf[:, :, ir_i] + imbuf[:, :, r_i])
#     feat = {}
#     feat['r_n'] = cn[:, :, r_i].mean()
#     feat['g_n'] = cn[:, :, g_i].mean()
#     feat['b_n'] = cn[:, :, b_i].mean()
#     feat['ir_n'] = cn[:, :, ir_i].mean()
#     feat['NDVI'] = ndvi.mean()
#     feat['SAVI'] = savi.mean()
#     feat['ir_rat'] = ir_rat.mean()
#     feat['i'] = (s / np.prod(s.shape[0:2])).mean()
#     feat['i_std'] = (s / np.prod(s.shape[0:2])).std()
#     feat['NDVI_std'] = ndvi.std()
#
#     return feat


def extract_patch_features(imbuf, mask):
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


def extract_all_features(ds, cs_gt_spatial_ref, cs_gt_dict, win_size=np.array((20, 20)),
                         win_origin='TL', win_rotation=27., plotFigures=False):  # , axis_signs=[1, 1]):
    # win_sizes - the OL and other window sizes in pixels
    # win_origin - the corner of the plot where the CS plot gps loc was recorded (BL,TL,TR,BR = SW,NW,NE,SE)
    #              the image (and any numpy matrix) origin is TL and increasing pixels L->R is W->E and T->B is N->S
    #              so eg if the plot origin is BL, one would want to start reading at TL to end at BL
    # win_rotation - the rotation of the plots (magnetic N offset from TN)
    # Note         A way around all this confusion may to be to layout plot co-ords in world co-ords and then convert to
    #              to pixel co-ords using world2pixel.  this would avoid all the fiddling and confusion in pixel space
    win_cnrs = np.array([[0, 0], [0, 1], [1, 1], [1, 0]])  # TL origin by default

    if win_origin == 'BL':
        win_cnrs[:, 1] = win_cnrs[:, 1] - 1
    elif win_origin == 'TL':
        # win_cnrs[:,1] = win_cnrs[:,1]-1
        win_cnrs = win_cnrs
    elif win_origin == 'TR':
        win_cnrs[:, 0] = win_cnrs[:, 0] - 1
        # win_cnrs[:,1] = win_cnrs[:,1]-1
    elif win_origin == 'BR':
        win_cnrs[:, 0] = win_cnrs[:, 0] - 1
        win_cnrs[:, 1] = win_cnrs[:, 1] - 1
    else:
        raise Exception('Unknown win_origin')

    # where is the right place to apply this?
    # win_cnrs[:, 0] = win_cnrs[:, 0] * axis_signs[0]
    # win_cnrs[:, 1] = win_cnrs[:, 1] * axis_signs[1]

    win_coord = win_cnrs.copy()

    theta = -np.radians(win_rotation)  # -ve angle for inverted y axis (valid for rotation<90)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s], [s, c]])
    # win_coords_ = [[], []]
    win_coord = win_cnrs * np.tile(win_size, [4, 1])  # note -1
    win_coord_ = np.dot(win_coord, R.transpose())
        # where is the right place to apply this?  also, it can be incorporated into R more neatly
        # win_coords[i][:, 0] = win_coords[i][:, 0] * axis_signs[0]
        # win_coords[i][:, 1] = win_coords[i][:, 1] * axis_signs[1]
        # win_coords_[i][:, 0] = win_coords_[i][:, 0] * axis_signs[0]
        # win_coords_[i][:, 1] = win_coords_[i][:, 1] * axis_signs[1]

    if plotFigures:
        f = pylab.figure()
        p = pylab.plot(win_coord[:, 0], win_coord[:, 1])
        pylab.plot(win_coord_[:, 0], win_coord_[:, 1], 'r')
        pylab.axis('equal')
        f.gca().invert_yaxis()

    win_mask = np.zeros(win_size)
    win_mask_size = []
    if win_rotation is not None:
        print 'rotating'

        # for win_mask, win_rotation, win_size in izip(win_masks, win_rotations, win_sizes):
        if True:
            mn = (np.min(win_coord_, 0))
            mx = (np.max(win_coord_, 0))
            win_size_r = np.int32(np.ceil(mx - mn))

            img = Image.fromarray(np.zeros(win_size_r))

            # Draw a rotated rectangle on the image.
            draw = ImageDraw.Draw(img)
            # rect = get_rect(x=120, y=80, width=100, height=40, angle=30.0)
            draw.polygon([tuple((p - mn)) for p in win_coord_], fill=1)
            # Convert the Image data to a numpy array.
            win_mask = np.asarray(img)
            # win_masks[i] = np.flipud(win_masks[i])  # ??? why is this necessary
        else:
            win_mask = ndimage.rotate(win_mask, win_rotation)
            # if axis_signs[0] < 0:
            #     win_masks[i] = np.fliplr(win_masks[i])
            # if axis_signs[1] < 0:
            #     win_masks[i] = np.flipud(win_masks[i])

            win_mask_size = win_mask.shape

    if plotFigures:
        pylab.figure()
        pylab.imshow(np.bool8(win_mask))

    geotransform = ds.GetGeoTransform()
    transform = osr.CoordinateTransformation(cs_gt_spatial_ref, osr.SpatialReference(ds.GetProjection()))
    i = 0
    plot_dict = {}
    # plotTagcDict = {}
    # class_labels = ['Pristine', 'Moderate', 'Severe']
    max_im_vals = np.zeros((4))
    for plot in cs_gt_dict.values():
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(plot['X'], plot['Y'])
        point.Transform(transform)  # xform into im projection
        (pixel, line) = world2Pixel(geotransform, point.GetX(), point.GetY())
        # not all the point fall inside the image
        if pixel >= 0 and line >= 0 and pixel < ds.RasterXSize and line < ds.RasterYSize:
            # win_size = win_sizes[0]
            # win_mask = win_masks[0]
            win_coord = win_coord_

            # ci = class_labels.index(plot['DegrClass'])

            # think about window size for co-ord xform
            mn = (np.min(win_coord, 0))
            mx = (np.max(win_coord, 0))
            win_size_r = np.int32(np.ceil(mx - mn))
            imbuf = np.zeros((win_size_r[0], win_size_r[1], 4), dtype=float)

            for b in range(1, 5):
                pixel_start = np.int(np.round((pixel + mn[0])))
                line_start = np.int(np.round((line + mn[1])))
                imbuf[:, :, b - 1] = ds.GetRasterBand(b).ReadAsArray(pixel_start, line_start, win_size_r[0],
                                                                     win_size_r[1])

            # imbuf[:, :, 3] = imbuf[:, :, 3] / 2  # hack for NGI XCALIB
            if np.all(imbuf == 0):
                print plot['ID'] + ": imbuf zero, assume NODATA ommitting"
                break
            # for b in range(0, 4):
            #     imbuf[:, :, b] = imbuf[:, :, b] / max_im_vals_[b]
            if not win_mask.shape == imbuf.shape[0:2]:
                print "error - mask and buf different sizes"
                raise Exception("error - mask and buf different sizes")
            feat = extract_patch_features(imbuf.copy(), win_mask)
            # feat['classi'] = ci
            # feat['class'] = class_labels[ci]
            # fields = ['TAGC', 'Z_P_AFRA', 'ALLOMETRY', 'HERB', 'LITTER']
            fields = plot.keys()
            for f in fields:
                feat[f] = cs_gt_dict[plot['ID']][f]
            feat['Xp'] = point.GetX()  # projected co-ords
            feat['Yp'] = point.GetY()

            feat['thumbnail'] = np.float32(imbuf.copy())
            plot_dict[plot['ID']] = feat
            # plotTagcDict[plot['PLOT']] = cs_gt_dict[plot['PLOT']]['TAGC']
            tmp = np.reshape(feat['thumbnail'], (np.prod(win_size_r), 4))
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
ds = gdal.OpenEx(samplingPtFile, gdal.OF_VECTOR)
if ds is None:
    print "Open failed./n"

lyr = ds.GetLayerByIndex(0)
lyr.ResetReading()
samplingPtSpatialRef = lyr.GetSpatialRef()

# gcpList = []
samplingPtDict = {}
for (i, feat) in enumerate(lyr):
    print '.',
    feat_defn = lyr.GetLayerDefn()
    f = {}
    for i in range(feat_defn.GetFieldCount()):
        field_defn = feat_defn.GetFieldDefn(i)
        f[field_defn.GetName()] = feat.GetField(i)
    geom = feat.GetGeometryRef()
    if geom is not None and (geom.GetGeometryType() == ogr.wkbPoint or geom.GetGeometryType() == ogr.wkbPoint25D):
        print "%s %.6f, %.6f" % (f['ID'], geom.GetX(), geom.GetY())
        f['geom'] = geom
        f['X'] = geom.GetX()
        f['Y'] = geom.GetY()
    else:
        print "no point geometry/n"
    # gcpList.append(f)
    samplingPtDict[f['ID']] = f
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


tchnuganuPlotDict = extract_all_features(tchnuganuDs, samplingPtSpatialRef, samplingPtDict, win_size=np.array([40, 40]),
                                  win_rotation=0., plotFigures=True)
tchnuganuDs = None

vdwSeweDs = gdal.OpenEx(vdwSewefonteinImageFile, gdal.OF_RASTER)
if vdwSeweDs is None:
    print "Open failed./n"
geotransform = vdwSeweDs.GetGeoTransform()
if not geotransform is None:
    print 'Origin = (', geotransform[0], ',', geotransform[3], ')'
    print 'Pixel Size = (', geotransform[1], ',', geotransform[5], ')'

vdwSewePlotDict = extract_all_features(vdwSeweDs, samplingPtSpatialRef, samplingPtDict, win_size=np.array([20, 20]),
                                   win_rotation=0., plotFigures=True)
vdwSeweDs = None

plotDict = tchnuganuPlotDict.copy()
plotDict.update(vdwSewePlotDict.copy())

featureName = 'NDVI'
featureVal = np.array([plot[featureName] for plot in plotDict.values()])
# featureVal = np.log10(np.array([plot[featureName] for plot in plotDict.values()]))
featureVal = np.float_power(10, np.array([plot[featureName] for plot in plotDict.values()]))
xp = np.array([plot['Xp'] for plot in plotDict.values()])
yp = np.array([plot['Yp'] for plot in plotDict.values()])
xyp = np.hstack((np.vstack(xp),np.vstack(yp)))

dxy = pdist(xyp, 'euclidean')
pylab.figure()
pylab.imshow(squareform(dxy))
pylab.colorbar()

# dfeat = pdist(np.vstack(featureVal), 'minkowski', p=1)
dfeat = pdist(np.vstack(featureVal), 'cityblock', p=1)
# dfeat = squareform(pdist(np.vstack(featureVal), 'cityblock'))
pylab.figure()
pylab.imshow(squareform(dfeat))
pylab.colorbar()

dxyGrid = np.linspace(dfeat.min(), dfeat.max(), 100)

pylab.figure()
classes = np.array([p['DegrClass'] for p in plotDict.values()])
class_labels = ['Severe', 'Moderate', 'Pristine'] #np.unique(classes)
class_num = [35, 25, 30]
featureValSub = np.array([])

for i, cl in enumerate(class_labels):
    idx = classes == cl
    #idx = idx[:class_num[i]]
    #allIdx.append(idx)
    dxy = pdist(xyp[idx], 'euclidean')
    dfeat = pdist(np.vstack(featureVal[idx]), 'minkowski', p=1)
    def variogram(distXy, distFeat, nbins=100):
        sortIdx = np.argsort(distXy)
        distXyGrid = np.linspace(distXy.min(), np.min((distXy.max(), 1500)), nbins+1)
        vg = np.zeros((nbins, 1))
        for i in range(0, nbins):
            idx = np.logical_and(distXy[sortIdx] > distXyGrid[i], distXy[sortIdx] <= distXyGrid[i+1])
            vg[i] = np.median(distFeat[sortIdx][idx])
        return vg, distXyGrid[:-1] + np.diff(distXyGrid).mean()
    featVg, dBins = variogram(dxy, dfeat, nbins=30)

    pylab.subplot(2, 2, i+1)
    pylab.plot(dBins, featVg, 'k-')
    pylab.xlabel('Distance (m)')
    pylab.ylabel('NDVI')
    pylab.title(cl)
    pylab.grid()


# plotDict.values()[0].keys()
