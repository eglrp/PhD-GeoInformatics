import gdal
import ogr
import numpy as np
import osr
import pylab
from scipy import stats as stats
from matplotlib import patches
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
imFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Digital Globe/056844553010_01/PCI Output/TOA and Haze/TOACorrected_056844553010_01_P001_OrthoPanSharpen_05644015.tif"
imFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Digital Globe/056844553010_01/056844553010_01_P001_MUL/PCI Ortho/03NOV18082012-M1BS-056844553010_01_P001_PCiOrtho.tif"

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


def extract_patch_features(imbuf):
    imbuf = np.float64(imbuf)/100.   # values are in percent?
    s = np.sum(imbuf, 2)
    cn = imbuf/np.tile(s[:,:,None], (1, 1, imbuf.shape[2]))
    b_i = 0
    g_i = 1
    r_i = 2
    ir_i = 3
    ndvi = (imbuf[:,:,ir_i] - imbuf[:,:,r_i])/(imbuf[:,:,ir_i] + imbuf[:,:,r_i])
    ir_rat = imbuf[:,:,ir_i]/imbuf[:,:,r_i]
    L = 0.05
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


def extract_all_features(ds, cs_gt_spatial_ref, cs_gt_dict, win_sizes=[np.array((8, 8)), np.array((42, 42))],
                     win_offsets=None):

    if win_offsets is None:   #setup defaults as determined to be optimal from experiments
        win_offsets = [np.array((0, 0)), np.array((0, 0))]
        for wi, win_size in enumerate(win_sizes):
            win_offsets[wi][0] = 0   #win_size[0]/2
            win_offsets[wi][1] = -win_size[1]

    transform = osr.CoordinateTransformation(cs_gt_spatial_ref, osr.SpatialReference(ds.GetProjection()))
    i = 0
    plot_dict = {}
    plotTagcDict = {}
    class_labels = ['OL','DST','ST']
    max_im_vals = np.zeros((4))
    for plot in cs_gt_dict.values():
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(plot['X'], plot['Y'])
        point.Transform(transform)  # xform into im projection
        (pixel, line) = world2Pixel(geotransform, point.GetX(), point.GetY())
        # not all the point fall inside the image
        if pixel >= 0 and line >=0 and pixel < ds.RasterXSize and line < ds.RasterYSize:
            win_size = win_sizes[0]
            win_offset = win_offsets[0]
            if 'DST' in plot['PLOT']:
                ci = 1
            elif 'ST' in plot['PLOT']:
                ci = 2
            elif 'OL' in plot['PLOT']:
                ci = 0
                win_size = win_sizes[1]
                win_offset = win_offsets[1]
            else:
                ci = 2
            imbuf = np.zeros((win_size[0], win_size[1], 4), dtype=float)
    
            for b in range(1,5):
                # imbuf[:, :, b-1] = ds.GetRasterBand(b).ReadAsArray(np.int(np.round(pixel))-(win_size[0]-1)/2,
                #                                                    np.int(np.round(line))-(win_size[1]-1)/2,
                #                                                    win_size[0], win_size[1])
                # imbuf[:, :, b-1] = ds.GetRasterBand(b).ReadAsArray(np.int(np.round(pixel))-(win_size[0]-1),
                #                                                    np.int(np.round(line))-(win_size[1]-1),
                #                                                    win_size[0], win_size[1])
                # imbuf[:, :, b-1] = ds.GetRasterBand(b).ReadAsArray(np.int(np.round(pixel)),
                #                                                    np.int(np.round(line)),
                #                                                    win_size[0], win_size[1])
                # imbuf[:, :, b-1] = ds.GetRasterBand(b).ReadAsArray(np.int(np.round(pixel))-(win_size[0]-1),
                #                                                    np.int(np.round(line)),
                #                                                    win_size[0], win_size[1])
                # this option performs best (bottom left / SW cnrs)
                # pixel_start = np.int(np.round((pixel - (win_size[0] - 1)/2) + win_offset[0]))
                # line_start = np.int(np.round((line - (win_size[1] - 1)/2) + win_offset[1]))
                pixel_start = np.int(np.round((pixel + win_offset[0])))
                line_start = np.int(np.round((line + win_offset[1])))
                imbuf[:, :, b-1] = ds.GetRasterBand(b).ReadAsArray(pixel_start, line_start, win_size[0], win_size[1])

            if np.all(imbuf==0):
                print "imbuf zero"
            # for b in range(0, 4):
            #     imbuf[:, :, b] = imbuf[:, :, b] / max_im_vals_[b]
            feat = extract_patch_features(imbuf)
            feat['classi'] = ci
            feat['class'] = class_labels[ci]
            fields = ['TAGC', 'Z_P_AFRA', 'ALLOMETRY', 'HERB', 'LITTER']
            for f in fields:
                feat[f] = cs_gt_dict[plot['PLOT']][f]
            feat['thumbnail'] = imbuf
            plot_dict[plot['PLOT']] = feat
                # plotTagcDict[plot['PLOT']] = cs_gt_dict[plot['PLOT']]['TAGC']
            tmp = np.reshape(imbuf, (np.prod(win_size), 4))
            # max_tmp = tmp.max(axis=0)
            max_tmp = np.percentile(tmp, 98., axis=0)
            max_im_vals[max_tmp > max_im_vals] = max_tmp[max_tmp > max_im_vals]
            # print plot['PLOT']
            i = i +1
        # else:
        #     print "x-" + plot['PLOT']
    
    # print i
    thumbs = np.array([plot['thumbnail'] for plot in plot_dict.values()])
    for thumb in thumbs:
        for b in range(0, 4):
            thumb[:,:,b] = thumb[:,:,b]/max_im_vals[b]
        thumb[:, :, 0] = thumb[:, :, 0] / 1.5
    
    return plot_dict


def scatterd(x, y, labels=None, class_labels=None, thumbnails=None, regress=True, xlabel=None, ylabel=None):

    if class_labels is None:
        class_labels = np.zeros(x.__len__())
    classes = np.unique(class_labels)
    colours = ['r','m','b','g','y','k','o']
    ylim = [np.min(y), np.max(y)]
    xlim = [np.min(x), np.max(x)]
    xd = np.diff(xlim)[0]
    yd = np.diff(ylim)[0]
    pylab.axis([np.min(x), np.max(x), np.min(y), np.max(y)])
    pylab.hold('on')
    ax = pylab.gca()
    handles = [0,0,0]

    for ci, (class_label, colour) in enumerate(zip(classes, colours[:classes.__len__()])):
        class_idx = class_labels == class_label
        if thumbnails is None:
            pylab.plot(x[class_idx], y[class_idx], colour + 'x')

        for xyi, (xx, yy) in enumerate(zip(x[class_idx], y[class_idx])):    #, np.array(plot_names)[class_idx]):
            if labels is not None:
                pylab.text(xx-.0015, yy-.0015, np.array(labels)[class_idx][xyi], fontdict={'size':8, 'color': colour})

            if thumbnails is not None:
                imbuf = np.array(thumbnails)[class_idx][xyi]
                ims = 20.
                extent = [xx - xd / (2 * ims), xx + xd / (2 * ims), yy - yd / (2 * ims), yy + yd / (2 * ims)]
                pylab.imshow(imbuf[:, :, 2::-1], extent=extent, aspect='auto')  # zorder=-1,
                handles[ci] = ax.add_patch(
                    patches.Rectangle((xx - xd / (2 * ims), yy - yd / (2 * ims)), xd / ims, yd / ims, fill=False,
                                      edgecolor=colour, linewidth=2.))
                # pylab.plot(mPixels[::step], dRawPixels[::step], color='k', marker='.', linestyle='', markersize=.5)
        if regress and classes.__len__() > 1:  # and False:
            (slope, intercept, r, p, stde) = stats.linregress(x[class_idx], y[class_idx])
            pylab.text(xlim[0] + xd*0.7, ylim[0] + yd*0.05*(ci + 2),
                       str.format('{1}: $R^2$ = {0:.2f}', np.round(r**2, 2), classes[ci]), fontdict={'size':10, 'color': colour})

    if regress:
        (slope, intercept, r, p, stde) = stats.linregress(x, y)
        pylab.text((xlim[0] + xd*0.4), (ylim[0] + yd*0.05), str.format('All: $R^2$ = {0:.2f}', np.round(r**2, 2)),
                   fontdict={'size':12})

    if xlabel is not None:
        pylab.xlabel(xlabel, fontdict={'size':12})
    if ylabel is not None:
        pylab.ylabel(ylabel, fontdict={'size':12})
    # pylab.ylabel(yf)
    if classes.__len__() > 1:
        pylab.legend(handles, classes, fontsize=12)

###########################################################################################################

# read in cs gt
ds = gdal.OpenEx(csGtFile, gdal.OF_VECTOR)
if ds is None:
    print "Open failed./n"

lyr = ds.GetLayerByIndex(0)
lyr.ResetReading()
cs_gt_spatial_ref = lyr.GetSpatialRef()

# gcpList = []
cs_gt_dict = {}
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
    cs_gt_dict[f['PLOT']] = f
print ' '

ds = None

# read in cs gt locations
ds = gdal.OpenEx(csGtGpsFile, gdal.OF_VECTOR)
if ds is None:
    print "Open failed./n"

lyr = ds.GetLayerByIndex(0)
lyr.ResetReading()
cs_gt_spatial_ref = lyr.GetSpatialRef()

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

# copy across gps locs to cs_gt_dict (they are more accurate than the ones in cs_gt_dict)
# [s for s in csGtGpsDict.keys() if "GH" in str(s)]
# [s for s in cs_gt_dict.keys() if "GH" in str(s)]
# [s for s in csGtGpsDict.keys() if "RHOL" in str(s)]
#
gpsInd = (np.zeros((1, csGtGpsDict.__len__())))

for k in cs_gt_dict.keys():
    # fixing some quirks in the keys
    gpsK = k
    gpsK = str(gpsK).replace('KADSS', 'KASS')
    gpsK = str(gpsK).replace('GHOL1_', 'RHOL1_')
    if gpsK == 'GHOL3_5':
        gpsK = 'GHOL3_'
    # gpsK = str(gpsK).replace('GHOL3_', 'GHOL3_5')
    if csGtGpsDict.has_key(gpsK):
        cs_gt_dict[k]['X'] = csGtGpsDict[gpsK]['X']
        cs_gt_dict[k]['Y'] = csGtGpsDict[gpsK]['Y']
        # print k
        gpsInd[0, csGtGpsDict.keys().index(gpsK)] = True
    else:
        print k + " not found in gps locs - " + gpsK

np.array(csGtGpsDict.keys())[np.logical_not(gpsInd[0])]

# Read in the image


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


if False:  # for MS image and excl OL
    plot_dict = extract_all_features(ds, cs_gt_spatial_ref, cs_gt_dict,
                                         win_sizes=[np.array((3, 3)), np.array((11, 11))])

    class_lab = np.array([plot['class'] for plot in plot_dict.values()])
    class_labi = np.array([plot['classi'] for plot in plot_dict.values()])
    ol_idx = class_lab == 'OL'

    plot_dict2 = dict((k, plot_dict[k]) for k in np.array(plot_dict.keys())[np.logical_not(ol_idx)])
    plot_dict = plot_dict2


# plot_dict = extract_all_features(ds, cs_gt_spatial_ref, cs_gt_dict)
plot_dict = extract_all_features(ds, cs_gt_spatial_ref, cs_gt_dict,
                                     win_sizes=[np.array((32, 32)), np.array((32, 32))])
# these are "best" options from window pos/size experiments below
plot_dict = extract_all_features(ds, cs_gt_spatial_ref, cs_gt_dict, win_sizes=[np.array([18,18]), np.array([32,32])])
plot_dict = extract_all_features(ds, cs_gt_spatial_ref, cs_gt_dict, win_sizes=[np.array([22,22]), np.array([50,50])],
                                 win_offsets=[np.array([0, 0]), np.array([0, 0])])
plot_dict = extract_all_features(ds, cs_gt_spatial_ref, cs_gt_dict, win_sizes=[np.array([18,18]), np.array([50,50])],
                                 win_offsets=[np.array([0, -18]), np.array([0, 0])])


# ndvi = np.array([plot['NDVI'] for plot in plot_dict.values()])
# gn = np.array([plot['g_n'] for plot in plot_dict.values()])
# std = np.array([plot['i_std'] for plot in plot_dict.values()])
# ir_rat = np.array([plot['ir_rat'] for plot in plot_dict.values()])
# tagc = np.array([plot['TAGC'] for plot in plot_dict.values()])
plotNames = plot_dict.keys()


yfields = ['TAGC', 'Z_P_AFRA', 'ALLOMETRY', 'HERB', 'LITTER']
xfields = ['NDVI', 'g_n', 'b_n', 'SAVI', 'ir_rat', 'NDVI_std']

pylab.figure()
x = np.array([plot['NDVI'] for plot in plot_dict.values()])
for yi, yf in enumerate(yfields):
    pylab.subplot(2, 3, yi+1)
    y = np.array([plot[yf] for plot in plot_dict.values()])
    scatterd(x, y, regress=True, xlabel='NDVI', ylabel=yf)


pylab.figure()
y = np.log10([plot['TAGC'] for plot in plot_dict.values()])
for xi, xf in enumerate(xfields):
    pylab.subplot(2, 3, xi+1)
    x = np.array([plot[xf] for plot in plot_dict.values()])
    scatterd(x, y, regress=True, xlabel=xf, ylabel='log10(TAGC)')


pylab.figure()
x = np.array([plot['NDVI'] for plot in plot_dict.values()])
class_lab = np.array([plot['class'] for plot in plot_dict.values()])
plot_names = plot_dict.keys()
for yi, yf in enumerate(yfields):
    pylab.subplot(2, 3, yi+1)
    y = np.array([plot[yf] for plot in plot_dict.values()])
    scatterd(x, y, labels=plot_names, class_labels=class_lab, regress=True, xlabel='NDVI', ylabel=yf)


#with images
pylab.figure()
x = np.array([plot['NDVI'] for plot in plot_dict.values()])
class_lab = np.array([plot['class'] for plot in plot_dict.values()])
plot_names = plot_dict.keys()
thumbs = np.array([plot['thumbnail'] for plot in plot_dict.values()])

for yi, yf in enumerate(['TAGC']):
    ax = pylab.subplot(1, 1, yi+1)
    y = np.log10(np.array([plot[yf] for plot in plot_dict.values()]))
    scatterd(x, y, labels=None, class_labels=class_lab, thumbnails=thumbs, xlabel='NDVI', ylabel='log10(TAGC)')
    pylab.grid(zorder=-1)



################################################################
# check effect of window offset
# win_sizes = [np.array((18, 18)), np.array((30, 30))]
#win_sizes = [np.array((8, 8)), np.array((42, 42))]
win_sizes = [np.array((3, 3)), np.array((11, 11))]

win_offsets = [np.array((0, 0)), np.array((0, 0))]
#xy_offset = np.arange(-64, 68, 8)
xy_offset = np.arange(-1.5, 2., 0.5)
res = np.zeros((xy_offset.__len__(), xy_offset.__len__(), 5))
for xi, xoff in enumerate(xy_offset):
    for yi, yoff in enumerate(xy_offset):
        for wi, win_size in enumerate(win_sizes):
            win_offsets[wi][0] = np.round(win_size[0] * xoff)
            win_offsets[wi][1] = np.round(win_size[1] * yoff)

        print win_offsets
        plot_dict = extract_all_features(ds, cs_gt_spatial_ref, cs_gt_dict, win_sizes=win_sizes, win_offsets=win_offsets)

        x = np.array([plot['NDVI'] for plot in plot_dict.values()])
        y = np.log10([plot['TAGC'] for plot in plot_dict.values()])
        class_lab = np.array([plot['class'] for plot in plot_dict.values()])
        class_labi = np.array([plot['classi'] for plot in plot_dict.values()])
        (slope, intercept, r, p, stde) = stats.linregress(x, y)
        res[yi, xi, 0] = r**2
        for ci,cclass in enumerate(['OL','DST','ST']):
            class_idx = class_lab == cclass
            (slope, intercept, r, p, stde) = stats.linregress(x[class_idx], y[class_idx])
            res[yi, xi, ci+1] = r**2
        class_idx = np.logical_not(class_lab == 'OL')
        (slope, intercept, r, p, stde) = stats.linregress(x[class_idx], y[class_idx])
        res[yi, xi, 4] = r**2

xgrid, ygrid = np.meshgrid(xy_offset, xy_offset)

res[:,:,0].max()
my, mx = np.unravel_index(res[:,:,0].argmax(), res[:,:,0].shape)
print xgrid[my, mx], ygrid[my, mx]


#fig = pylab.figure()
res_lab = ['All','OL','DST','ST','*ST']
fig = plt.figure()
for ci in range(0, res.shape[2]):
    ax = fig.add_subplot(2,3,ci+1, projection='3d')
    ax.plot_surface(xgrid, ygrid, res[:,:,ci])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(res_lab[ci])
# ax = fig.add_subplot(132, projection='3d')
# ax.plot_surface(xgrid, ygrid, res[:,:,1])
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax = fig.add_subplot(133, projection='3d')
# ax.plot_surface(xgrid, ygrid, res[:,:,2])
# ax.set_xlabel('x')
# ax.set_ylabel('y')


pylab.figure()
for i in range(0, res.shape[2]):
    pylab.subplot(2,3,i+1)
    pylab.imshow(res[:,:,i], extent=[-1.5, 1.5, 1.5, -1.5])
    pylab.colorbar()
    pylab.title(res_lab[i])

# check default
plot_dict = extract_all_features(ds, cs_gt_spatial_ref, cs_gt_dict, win_sizes=win_sizes)
x = np.array([plot['NDVI'] for plot in plot_dict.values()])
y = np.log10([plot['TAGC'] for plot in plot_dict.values()])
(slope, intercept, r, p, stde) = stats.linregress(x, y)
print r**2


#############################################################################33
# check effect of changing window size

# for differently placed windows
# win_sizes = [np.array((3, 3)), np.array((17, 17))]
win_sizes = [np.array((3, 3)), np.array((11, 11))]
win_offsets = [np.array((0, 0)), np.array((0, 0))]
incr_array = np.arange(2, 50, 2)
res = np.zeros((incr_array.__len__(), 5))
for i, incr in enumerate(incr_array):
    win_sizes_ = win_sizes + incr
    for wi, win_size in enumerate(win_sizes_):
        win_offsets[wi][0] = 0 #-win_size[0]/2    # CHHANGE THESE
        win_offsets[wi][1] = -win_size[1]
    plot_dict = extract_all_features(ds, cs_gt_spatial_ref, cs_gt_dict, win_sizes=win_sizes_,
                                     win_offsets = win_offsets)

    x = np.array([plot['NDVI'] for plot in plot_dict.values()])
    y = np.log10([plot['TAGC'] for plot in plot_dict.values()])
    class_lab = np.array([plot['class'] for plot in plot_dict.values()])
    (slope, intercept, r, p, stde) = stats.linregress(x, y)
    res[i, 0] = r ** 2
    for ci, cclass in enumerate(['OL', 'DST', 'ST']):
        class_idx = class_lab == cclass
        (slope, intercept, r, p, stde) = stats.linregress(x[class_idx], y[class_idx])
        res[i, ci + 1] = r ** 2
    class_idx = np.logical_not(class_lab == 'OL')
    (slope, intercept, r, p, stde) = stats.linregress(x[class_idx], y[class_idx])
    res[i, 4] = r ** 2

res_lab = ['All','OL','DST','ST','*ST']
pylab.figure()
pylab.plot(np.c_[incr_array+3, incr_array+17, incr_array+3, incr_array+3, incr_array+3], res, 'x-')
pylab.legend(res_lab) #['all', 'OL', 'ST+DST'])
pylab.title('0,-1')


#Note: testing above with different window offsets, the most likely offsets are 0,0 ors 0,-1.  But then 0,0 gives an
# overall r2>5 for windows of ~30 (?)
#NB Note: The OL individual performance is significantly improved for offset 0,0 !!  Overall though 0,-1 gives better results

#check with best params
plot_dict = extract_all_features(ds, cs_gt_spatial_ref, cs_gt_dict, win_sizes=[np.array([18,18]), np.array([32,32])])
plot_dict = extract_all_features(ds, cs_gt_spatial_ref, cs_gt_dict, win_sizes=[np.array([22,22]), np.array([50,50])],
                                 win_offsets=[np.array([0, 0]), np.array([0, 0])])

x = np.array([plot['NDVI'] for plot in plot_dict.values()])
y = np.log10([plot['TAGC'] for plot in plot_dict.values()])
class_lab = np.array([plot['class'] for plot in plot_dict.values()])
(slope, intercept, r, p, stde) = stats.linregress(x, y)
print r ** 2
class_idx = class_lab == 'OL'
(slope, intercept, r, p, stde) = stats.linregress(x[class_idx], y[class_idx])
print r ** 2
(slope, intercept, r, p, stde) = stats.linregress(x[np.logical_not(class_idx)], y[np.logical_not(class_idx)])
print r ** 2



########################################################################################################
# Simulate smaller windows than 5m\

#############################################################################33
# check effect of changing window size

# for differently placed windows
win_sizes = [np.array((9, 9)), np.array((42, 42))]
win_offsets = [np.array((0, 0)), np.array((0, 0))]
incr_array = np.arange(-1, -9, -1)
res = np.zeros((incr_array.__len__(), 5))
for i, incr in enumerate(incr_array):
    win_sizes_ = win_sizes + incr
    for wi, win_size in enumerate(win_sizes_):
        win_offsets[wi][0] = 0 #-win_size[0]/2    # CHHANGE THESE
        win_offsets[wi][1] = -win_size[1]
    plot_dict = extract_all_features(ds, cs_gt_spatial_ref, cs_gt_dict, win_sizes=win_sizes_,
                                     win_offsets = win_offsets)

    x = np.array([plot['NDVI'] for plot in plot_dict.values()])
    y = np.log10([plot['TAGC'] for plot in plot_dict.values()])
    class_lab = np.array([plot['class'] for plot in plot_dict.values()])
    (slope, intercept, r, p, stde) = stats.linregress(x, y)
    res[i, 0] = r ** 2
    for ci, cclass in enumerate(['OL', 'DST', 'ST']):
        class_idx = class_lab == cclass
        (slope, intercept, r, p, stde) = stats.linregress(x[class_idx], y[class_idx])
        res[i, ci + 1] = r ** 2
    class_idx = np.logical_not(class_lab == 'OL')
    (slope, intercept, r, p, stde) = stats.linregress(x[class_idx], y[class_idx])
    res[i, 4] = r ** 2

res_lab = ['All','OL','DST','ST','*ST']
pylab.figure()
pylab.plot(0.6*np.c_[incr_array+9, incr_array+42, incr_array+9, incr_array+9, incr_array+9], res, 'x-')
pylab.legend(res_lab) #['all', 'OL', 'ST+DST'])
pylab.title('0,-1')

#plot for report
res_lab = ['All','OL','DST','ST','*ST']
pylab.figure()
pylab.plot(0.6*np.c_[incr_array+9], res[:, 0], 'x-')
pylab.xlabel('Plot size (m)', fontsize=12)
pylab.ylabel('$R^2$ for NDVI vs log$_{10}$(TAGC)', fontsize=12)
#pylab.legend(res_lab) #['all', 'OL', 'ST+DST'])
pylab.title('Effect of plot size', fontsize=12)
pylab.grid()



if False:

    #for default window placement
    win_sizes = [np.array((4, 4)), np.array((16, 16))]
    incr_array = np.arange(0, 25, 2)
    res = np.zeros((incr_array.__len__(), 3))
    for i, incr in enumerate(incr_array):
        # win_sizes += incr
        plot_dict = extract_all_features(ds, cs_gt_spatial_ref, cs_gt_dict, win_sizes=win_sizes+incr)

        x = np.array([plot['NDVI'] for plot in plot_dict.values()])
        y = np.log10([plot['TAGC'] for plot in plot_dict.values()])
        class_lab = np.array([plot['class'] for plot in plot_dict.values()])
        (slope, intercept, r, p, stde) = stats.linregress(x, y)
        res[i, 0] = r ** 2
        class_idx = class_lab == 'OL'
        (slope, intercept, r, p, stde) = stats.linregress(x[class_idx], y[class_idx])
        res[i, 1] = r ** 2
        (slope, intercept, r, p, stde) = stats.linregress(x[np.logical_not(class_idx)], y[np.logical_not(class_idx)])
        res[i, 2] = r ** 2


    pylab.figure()
    pylab.plot(np.c_[incr_array+4, incr_array+16, incr_array+4], res, 'x-')
    pylab.legend(['all', 'OL', 'ST+DST'])

    #check with best params
    plot_dict = extract_all_features(ds, cs_gt_spatial_ref, cs_gt_dict, win_sizes=[np.array([18,18]), np.array([30,30])])

    x = np.array([plot['NDVI'] for plot in plot_dict.values()])
    y = np.log10([plot['TAGC'] for plot in plot_dict.values()])
    (slope, intercept, r, p, stde) = stats.linregress(x, y)
    print r**2


    #for actual OL window size
    win_size= np.array((3, 3))
    incr_array = np.arange(5, 20, 2)
    res = np.zeros(incr_array.__len__())
    for i, incr in enumerate(incr_array):
        win_size += incr
        plot_dict = extract_all_features(ds, cs_gt_spatial_ref, cs_gt_dict, win_sizes=[win_size, np.array([42, 42])])

        x = np.array([plot['NDVI'] for plot in plot_dict.values()])
        y = np.log10([plot['TAGC'] for plot in plot_dict.values()])
        (slope, intercept, r, p, stde) = stats.linregress(x, y)
        res[i] = r**2

    pylab.figure()
    pylab.plot(3 + incr_array, res, 'kx-')




#todo x check if atcor refl vals are in % - yes they are
#todo x define ol, dst, st classes and visualise
#todo research "regression analysis" in general and in python
#todo x plot with labels and or images to try and figure out what are problem plots and patterns
#todo show pattern of correlation with spatial accuracy and window size - somehow simulate different accuracy and plot size conditions
#todo make a separate feature of area around window - should give some insight into eg soil condition
#todo brownness or soil index
#todo make features that use classification within window eg perhaps ttl vegetation pixels, or mean soil colour and mean veg colour
#todo investigate visualisation with qgis/arc python libraries
#todo x FETCH MY PKG FROM PATTI
#todo x what about xects?  were they used in calculating TAGC? No
#todo consider other no visual variables like altitude and slope, ttl sun hours.
#todo what about finding the ttl amount of sun various areas receive based on DEM?
#todo try doing this with aeral imagery - need NIR though...
#todo x regress on individual classes
#todo refactor code to have functions like scatterd for general scattering
#todo NB go back to orthorectification - why do the errors in my image appear a lot worse than the PCI reported errors?
#todo NB go back to atcor
# Notes on MP's gt
#- OL plots are 25x25m, ST and DST are 5x5m
#- Allometry was done inside these plots - no transects
#- "All C data and tree architecture data was log10 transformed."

# Looking at anomalous plots
#- RHST21 med low NDVI, high TAGC
#  It is on UL of a big tree that is probably missed in feature extraction but not MP GT.  Spatial issues again.
#- RHDST23 mostly bare ground (low NDVI) but high TAGC
#  similar to RHDST6 - there are med dense number of bushes/trees in the vicinity but prob not in the image block.  Again
#  this points to a spatial error and the possibility of improving with a feature that considers surrounding veg.
#- RHDST5-6 med/high TAGC but mostly bare ground / low NDVI
#  RHDST5 has some little low NDVI bossies around - does not make sense that NDVI is so low.  Again this suggests we need
#  other measures than just NDVI
#  RHDST6 has quite a few bushes and trees in the vicinity but perhaps not in its block therefore lowish NDVI.  This could
#  be due to spatial error and could also be improved with a feature that considers the surrounding veg.
#- RHOL15, RHOL21 - high TAGC, low NDVI
#  RHOL15 sits on the UR edge of a big tree which is prob excluded from my analysis but not from MP GT i.e. spatial error
#  RHOL21 sits near big tree but it has lowish NDVI.  There also appear to be some little shrubby things which are prob
#  not picked up in image analysis i.e. we need to look at more than just NDVI and perhaps include texture which would show
#  shrubby things
#- RHST17 - high NDVI, med/low TAGC
#  its not clear from the image why it should have low TAGC, it is covered in veg
#- To summarise: spatial accuracy is NB and part of the cause of our inaccuracies.  Including other features in the
#  regression may help improve things.

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
