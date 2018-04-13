import gdal
import numpy as np
import os
import pylab
import scipy.stats
from matplotlib import pyplot

spotFileName = "D:\Data\Development\Projects\PhD GeoInformatics\Data\SPOT\S131022114824832\Orthorectification\oATCORCorrected_METADATA_00812098_AutoGCP_NgiFormat.tif"  #improved orthorect
modisFileName = 'D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.A2010025.h19v12.005.2010043064233.Lo23.RGBN.tif'
calibRootDir = 'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\XCalib Experiments\Calibrated\\'

mosaicFileName = "D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\XCalib Experiments\Calibrated\w11p1\w11p1_Mosaic10m.tif"
reflScale = 5000.

# for subDir in os.listdir(calibRootDir):
#     if os.path.isdir(os.path.join(calibRootDir, subDir)):
#         print subDir
#         mosaicFileName = subDir + '_Mosaic10m.tif'

spotDs = gdal.Open(spotFileName, gdal.GA_ReadOnly)
ngiCalibDs = gdal.Open(mosaicFileName, gdal.GA_ReadOnly)

spotIm = np.zeros((spotDs.RasterYSize, spotDs.RasterXSize, 3), dtype = np.float32)
ngiCalibIm = np.zeros((ngiCalibDs.RasterYSize, ngiCalibDs.RasterXSize, 4), dtype = np.float32)

# to fix some weirdness in NGI calib mosaic nodata
ngiCalibMask = np.bool8(ngiCalibDs.GetRasterBand(1).GetMaskBand().ReadAsArray())
for b in range(1,4):
    ncd = ngiCalibDs.GetRasterBand(b).ReadAsArray()
    ngiCalibMask = (ngiCalibMask & (ncd > 0))

spotMask = np.bool8(spotDs.GetRasterBand(1).GetMaskBand().ReadAsArray())

if False:
    pylab.figure()
    ax = pylab.subplot(1,2,1)
    pylab.imshow(spotMask)
    pylab.subplot(1,2,2, sharex=ax, sharey=ax)
    pylab.imshow(ngiCalibMask)

for b in range(1,4):
    sb = spotDs.GetRasterBand(b).ReadAsArray()
    sb[~spotMask] = 0
    spotIm[:,:,b-1] = np.float32(sb)/reflScale

for b in range(1, 5):
    ncd = ngiCalibDs.GetRasterBand(b).ReadAsArray()
    ncd[~ngiCalibMask] = 0
    ngiCalibIm[:,:,b-1] = np.float32(ncd)/reflScale

# find ngi extents in spot
spotGeoTransform = spotDs.GetGeoTransform()
ngiCalibGeoTransform = ngiCalibDs.GetGeoTransform()

spotInvGeoTransform = gdal.InvGeoTransform(spotGeoTransform)
(ulx, uly) = gdal.ApplyGeoTransform(ngiCalibGeoTransform, 0, 0)
(brx, bry) = gdal.ApplyGeoTransform(ngiCalibGeoTransform, ngiCalibDs.RasterXSize, ngiCalibDs.RasterYSize)
(refUlx, refUly) = np.int32(gdal.ApplyGeoTransform(spotInvGeoTransform, ulx, uly))
(refBrx, refBry) = np.int32(gdal.ApplyGeoTransform(spotInvGeoTransform, brx, bry))

spotSubIm = spotIm[refUly:refBry, refUlx:refBrx, :]
spotSubMask = spotMask[refUly:refBry, refUlx:refBrx]
mask = ngiCalibMask & spotSubMask

spotDs = None
ngiCalibDs = None

pylab.close('all')
spotBands = (0, 1, 2)
ngiBands = (3, 0, 1)
colors = ( 'orange', 'r', 'g')
legend = ('NIR', 'Red', 'Green')
step = 1000
fn = pylab.gcf().number
pn = 1
#
#
# fontSize = 24.
# mpl.rcParams.update({'font.size': fontSize})

f2 = pylab.figure('Calibrated')
f2.set_size_inches(16., 5.25, forward=True)

## NB
ePixels = []
slopes = []
intercepts = []
for sb, nb in zip(spotBands, ngiBands):
    sPixels = spotSubIm[:, :, sb][mask]
    dCalibPixels = ngiCalibIm[:, :, nb][mask]
    e = sPixels - dCalibPixels
    ePixels.append(e)

    (slope, intercept, r, p, stde) = scipy.stats.linregress(sPixels.flatten(), dCalibPixels.flatten())
    slopes.append(slope)
    intercepts.append(intercept)
    pylab.figure('Calibrated')
    pylab.subplot(1, 3, pn)
    pylab.plot(sPixels[::step], dCalibPixels[::step], color='k', marker='.', linestyle='', markersize=.5)
    pylab.hold(True)
    m = 1. #np.min([np.max([np.percentile(sPixels[::step], 99.99), np.percentile(dCalibPixels[::step], 99.99)]), 1.])
    oneLineH, = pylab.plot([0, m], [0, m], color='k', marker='', linestyle='--', label='1:1')
    pylab.gca().set_xlim([0, m])
    pylab.gca().set_ylim([0, m])
    xl = pylab.gca().get_xlim()
    yl = pylab.gca().get_ylim()
    pylab.text((xl[0] + np.diff(xl)*0.05)[0], (yl[0] + np.diff(yl)*0.8)[0], str.format('$R^2$ = {0:.2f}',
                                                                                       np.round(r**2, 2)))
    # pylab.text((xl[0] + np.diff(xl)*0.05)[0], (yl[0] + np.diff(yl)*0.9)[0], str.format('y = {0:.2f} x + {1:.3f}',
    #                                                                                    slope, intercept))
    pylab.title(legend[sb])
    # pylab.xlabel(r'SPOT5 $\rho_t$')
    # pylab.ylabel(r'DMC $\rho_t$')
    pylab.xlabel(r'SPOT5 surface refl.')
    pylab.ylabel(r'DMC surface refl.')
    pylab.grid('on')
    pylab.legend(handles=[oneLineH], prop={'size':20}, loc=4)
    pylab.tight_layout()
    pyplot.locator_params(axis='y', nbins=5)#to specify number of ticks on both or any single axes
    pyplot.locator_params(axis='x', nbins=5)
    #pylab.hold('on')
    pn += 1
#
# f1.savefig('C:/Data/Development/Projects/PhD GeoInformatics/Docs/My Docs/Thesis/Retrieval of Surface Reflectance '
#            'from Aerial Imagery/Figure 12 - DMC DN and SPOT5 Surface Reflectance Correlation V2.eps', dpi=1200)
# f2.savefig('C:/Data/Development/Projects/PhD GeoInformatics/Docs/My Docs/Thesis/Retrieval of Surface Reflectance '
#            'from Aerial Imagery/Figure 13 - DMC Homogenised and SPOT5 Surface Reflectance Correlation V2.eps', dpi=1200)

print str.format("N (num points plotted): {0}, step: {1}", sPixels[::step].__len__(), step)


# errors
# std dev(abs error), std dev(rms) etc dont really make sense and one measure of error variation is good enough i.e. std dev
# so i.e. use dont use the rhs () bracket vals but rather the third and last rows for error variation
ePixels = np.array(np.float64(ePixels))*100
print str.format("Per band Abs error (%): {0} ({1})", np.abs(ePixels).mean(axis=1).round(decimals=2), np.abs(ePixels).std(axis=1).round(decimals=2))
print str.format("Per band RMS error (%): {0} ({1})", np.sqrt((ePixels**2).mean(axis=1)).round(decimals=2), np.sqrt((ePixels**2).std(axis=1)).round(decimals=2))  # does root-std-sqr make any sense?
print str.format("Per band Std Dev error (%): {0}", ePixels.std(axis=1).round(decimals=2))
print str.format("Mean Abs error (%): {0:.2f} ({1:.2f})", np.abs(ePixels).mean(), np.abs(ePixels).std())
print str.format("Mean RMS error (%): {0:.2f} ({1:.2f})", np.sqrt((ePixels**2).mean()), np.sqrt((ePixels**2).std()))  # does root-std-sqr make any sense?
print str.format("Std Dev error (%): {0:.2f}", ePixels.std())
