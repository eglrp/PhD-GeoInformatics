#show scatter plots of DMC vs MODIS before and after calibration
# %gui qt
# %matplotlib qt

import numpy as np
import pylab
import os
import scipy.stats
import gdal
import scipy.ndimage
# from mayavi import mlab
import matplotlib.image
from matplotlib import pyplot

reflScale = 5000.
if True:
    numBands = 4

    modisFileName = "D:/Data/Development/Projects/PhD GeoInformatics/Data/NGI/Cross Calibration/Mosaics/Modis_CropStudyArea_Float32.tif"
    ngiRawFileName = "D:/Data/Development/Projects/PhD GeoInformatics/Data/NGI/Cross Calibration/Mosaics/StudyAreaUncalibratedMosaic500mCs_CropStudyArea_Float32.tif"
    ngiCalibFileName = "D:/Data/Development/Projects/PhD GeoInformatics/Data/NGI/Cross Calibration/Mosaics/StudyAreaMosaic500mCs_CropStudyArea_Float32.tif"
    #ngiCalibFileName = "F:/MSc GeoInformatics/Data/NGI/Cross Calibration/Mosaics/MasterXCalibMosaic500mCsSeamLineFix_CropStudyArea_Float32.tif"

    modisDs = gdal.Open(modisFileName, gdal.GA_ReadOnly)
    ngiRawDs = gdal.Open(ngiRawFileName, gdal.GA_ReadOnly)
    ngiCalibDs = gdal.Open(ngiCalibFileName, gdal.GA_ReadOnly)

    modisIm = np.zeros((modisDs.RasterYSize, modisDs.RasterXSize, 4), dtype = float)
    ngiRawIm = np.zeros((ngiRawDs.RasterYSize, ngiRawDs.RasterXSize, 4), dtype = float)
    ngiCalibIm = np.zeros((ngiCalibDs.RasterYSize, ngiCalibDs.RasterXSize, 4), dtype = float)

    ngiCalibMask = np.bool8(ngiCalibDs.GetRasterBand(1).GetMaskBand().ReadAsArray())

    for b in range(1,5):
        mb = modisDs.GetRasterBand(b).ReadAsArray()
        modisIm[:,:,b-1] = mb/reflScale
        nrb = ngiRawDs.GetRasterBand(b).ReadAsArray()
        ngiRawIm[:,:,b-1] = nrb #/reflScale
        ncd = ngiCalibDs.GetRasterBand(b).ReadAsArray()
        ngiCalibIm[:,:,b-1] = ncd/reflScale

    modisDs = None
    ngiRawDs = None
    ngiCalibDs = None

        #ngiCalibIm[:,:,b-1] = scipy.ndimage.filters.gaussian_filter(ngiCalibIm[:,:,b-1], 0.15)
else: #exported from acrmap mosaics to same res
    numBands = 3
    modisFileName = "F:/MSc GeoInformatics/Data/NGI/Cross Calibration/Mosaics/ModisStudyArea.png"
    ngiRawFileName = "F:/MSc GeoInformatics/Data/NGI/Cross Calibration/Mosaics/NgRawStudyAreaDs.png"
    ngiCalibFileName = "F:/MSc GeoInformatics/Data/NGI/Cross Calibration/Mosaics/NgiXCalibStudyAreaDs.png"

    modisIm = np.float32(pylab.imread(modisFileName))
    modisIm = modisIm/modisIm.max()
    ngiRawIm = np.float32(pylab.imread(ngiRawFileName))
    ngiRawIm = ngiRawIm/ngiRawIm.max()
    ngiCalibIm = np.float32(pylab.imread(ngiCalibFileName))
    ngiCalibIm = ngiCalibIm/modisIm.max()

errorIm = np.abs(ngiCalibIm - modisIm)

#pylab.close('all')
pylab.figure()
pylab.subplot(2, 2, 1)
pylab.imshow(modisIm[:, :, :3]/modisIm[:, :, :3].max())
pylab.subplot(2, 2, 2)
pylab.imshow(ngiRawIm[:, :, :3]/ngiRawIm[:, :, :3].max())
pylab.subplot(2, 2, 3)
pylab.imshow(ngiCalibIm[:, :, :3]/ngiCalibIm[:, :, :3].max())
pylab.subplot(2, 2, 4)
pylab.imshow(errorIm[:, :, :3]/errorIm[:, :, :3].max())

pylab.close('all')
bands = (3, 0, 1, 2)
colors = ('r', 'g', 'b', 'orange')
legend = ('Red', 'Green', 'Blue', 'NIR')
step = 5
fn = pylab.gcf().number
pn = 1

fontSize = 24.
import matplotlib as mpl
from matplotlib import pyplot
mpl.rcParams.update({'font.size': fontSize})

f1 = pylab.figure('Uncalibrated')
f1.set_size_inches(12, 9, forward=True)
f2 = pylab.figure('Calibrated')
f2.set_size_inches(12, 9, forward=True)

ePixels = []
for b in bands:
    mPixels = modisIm[:, :, b][ngiCalibMask]
    dRawPixels = ngiRawIm[:, :, b][ngiCalibMask]
    dCalibPixels = ngiCalibIm[:, :, b][ngiCalibMask]

    ePixels.append(mPixels - dCalibPixels)

    #dCalibPixels = scipy.ndimage.filters.gaussian_filter(ngiCalibIm[:, :, b], 2)
    #dCalibPixels = dCalibPixels.flatten()

    (slope, intercept, r, p, stde) = scipy.stats.linregress(mPixels, dRawPixels)
    pylab.figure('Uncalibrated')
    pylab.subplot(2, 2, pn)
    pylab.plot(mPixels[::step], dRawPixels[::step], color='k', marker='.', linestyle='', markersize=.5)
    xl = pylab.gca().get_xlim()
    yl = pylab.gca().get_ylim()
    pylab.text((xl[0] + np.diff(xl)*0.05)[0], (yl[0] + np.diff(yl)*0.8)[0], str.format('$R^2$ = {0:.2f}',
                                                                                       np.round(r**2, 2)))
    #pylab.xlabel(r'MODIS $\rho_t$')
    pylab.xlabel(r'MODIS surface refl.')
    pylab.ylabel('DMC DN')
    pylab.title(legend[b])
    pylab.grid('on')
    pylab.tight_layout()
    pyplot.locator_params(axis='y', nbins=5)#to specify number of ticks on both or any single axes
    pyplot.locator_params(axis='x', nbins=6)
    #pylab.hold('on')

    (slope, intercept, r, p, stde) = scipy.stats.linregress(mPixels, dCalibPixels)
    pylab.figure('Calibrated')
    pylab.subplot(2, 2, pn)
    pylab.plot(mPixels[::step], dCalibPixels[::step], color='k', marker='.', linestyle='', markersize=.5)
    pylab.hold(True)
    m = np.max([mPixels.max(), dCalibPixels.max()])
    oneLineH, = pylab.plot([0, m], [0, m], color='k', marker='', linestyle='--', label='1:1')
    pylab.gca().set_xlim([0, m])
    pylab.gca().set_ylim([0, m])
    xl = pylab.gca().get_xlim()
    yl = pylab.gca().get_ylim()
    pylab.text((xl[0] + np.diff(xl)*0.05)[0], (yl[0] + np.diff(yl)*0.8)[0], str.format('$R^2$ = {0:.2f}',
                                                                                       np.round(r**2, 2)))
    pylab.title(legend[b])
    #pylab.xlabel(r'MODIS $\rho_t$')
    pylab.xlabel(r'MODIS surface refl.')
    #pylab.ylabel(r'DMC $\rho_t$')
    pylab.ylabel(r'DMC surface refl.')
    pylab.legend(handles=[oneLineH], prop={'size':20}, loc=4)
    pylab.grid('on')

    pylab.tight_layout()
    pyplot.locator_params(axis='y', nbins=5)#to specify number of ticks on both or any single axes
    pyplot.locator_params(axis='x', nbins=6)

    #pylab.hold('on')
    pn += 1

    bins = 50

    dRawH, dRawB = np.histogram(dRawPixels, bins = bins)
    mH, mB = np.histogram(mPixels, bins = 100)
    dCalibH, dCalibB = np.histogram(dCalibPixels, bins = mB)
    mX = mB[:-1] + np.diff(mB[:2])              #bin centers
    dRawX = dRawB[:-1] + np.diff(dRawB[:2])     #bin centers
    pylab.figure(fn + 3)
    pylab.subplot(2, 2, b+1)
    pylab.plot(mX, mH, linestyle='-', color=colors[b])
    pylab.hold('on')
    pylab.plot(mX, dCalibH, linestyle=':', color=colors[b])
    pylab.plot(dRawX, dRawH, linestyle='--', color=colors[b])
    pylab.tight_layout()
    #title(str.format('Calib r = {}', r))

f1.savefig('C:/Data/Development/Projects/MSc GeoInformatics/Docs/My Docs/Thesis/Retrieval of Surface Reflectance'
           ' from Aerial Imagery/Figure 8 - DMC DN and MODIS Surface Reflectance Correlation.eps', dpi=1200)
f2.savefig('C:/Data/Development/Projects/MSc GeoInformatics/Docs/My Docs/Thesis/Retrieval of Surface Reflectance'
           ' from Aerial Imagery/FFigure 9 - DMC Homogenised and MODIS Surface Reflectance Correlation.eps', dpi=1200)

print str.format("N (num points plotted): {0}, step: {1}", mPixels[::step].__len__(), step)

# errors
# somehow the rms and sdd vals look identical when rounded to 2 decimal places.  How the hell can this be?
ePixels = np.array(ePixels)*100
print str.format("Per band Abs error (%): {0} ({1})", np.abs(ePixels).mean(axis=1), np.abs(ePixels).std(axis=1))
print str.format("Per band RMS error (%): {0} ({1})", np.sqrt((ePixels**2).mean(axis=1)), np.sqrt((ePixels**2).std(axis=1)))  # does root-std-sqr make any sense?
print str.format("Per band Std Dev error (%): {0}", ePixels.std(axis=1))
print str.format("Mean Abs error (%): {0} ({1})", np.abs(ePixels).mean(), np.abs(ePixels).std())
print str.format("Mean RMS error (%): {0} ({1})", np.sqrt((ePixels**2).mean()), np.sqrt((ePixels**2).std()))  # does root-std-sqr make any sense?
print str.format("Std Dev error (%): {0}", ePixels.std())
# std dev = rms when mean = 0 which is more or less the case above


##
#make a scatter plot of SPOT vs NGI

spotFileName = "D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\spotInt16_2.tif"  #old version
spotFileNameNew = "D:\Data\Development\Projects\PhD GeoInformatics\Data\SPOT\S131022114824832\Orthorectification\oATCORCorrected_METADATA_00812098_AutoGCP_NgiFormat.tif"  #improved orthorect
ngiRawFileName = "D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\Cross Calibration\Mosaics\StudyAreaUncalibratedMosaicSpotMask.tif"
ngiCalibFileName = "D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\StudyAreaXCalibMosaicSpotMask.tif"

spotDs = gdal.Open(spotFileName, gdal.GA_ReadOnly)
spotNewDs = gdal.Open(spotFileNameNew, gdal.GA_ReadOnly)
ngiRawDs = gdal.Open(ngiRawFileName, gdal.GA_ReadOnly)
ngiCalibDs = gdal.Open(ngiCalibFileName, gdal.GA_ReadOnly)

spotIm = np.zeros((spotDs.RasterYSize, spotDs.RasterXSize, 3), dtype = np.float32)
spotNewIm = np.zeros((spotNewDs.RasterYSize, spotNewDs.RasterXSize, 3), dtype = np.float32)
ngiRawIm = np.zeros((ngiRawDs.RasterYSize, ngiRawDs.RasterXSize, 3), dtype = np.float32)
ngiCalibIm = np.zeros((ngiCalibDs.RasterYSize, ngiCalibDs.RasterXSize, 3), dtype = np.float32)
diffIm = np.zeros((ngiCalibDs.RasterYSize, ngiCalibDs.RasterXSize, 3), dtype = np.float32)
diffNewIm = np.zeros((ngiCalibDs.RasterYSize, ngiCalibDs.RasterXSize, 3), dtype = np.float32)

ngiCalibMask = np.bool8(ngiCalibDs.GetRasterBand(1).GetMaskBand().ReadAsArray())

mask = np.bool8(spotDs.GetRasterBand(1).GetMaskBand().ReadAsArray())
mask2 = np.bool8(ngiCalibDs.GetRasterBand(1).GetMaskBand().ReadAsArray())
mask3 = np.bool8(spotNewDs.GetRasterBand(1).GetMaskBand().ReadAsArray())
mask = mask & mask2 & mask3

for b in range(1,4):
    sb = spotDs.GetRasterBand(b).ReadAsArray()
    sb[~mask] = 0
    spotIm[:,:,b-1] = np.float32(sb)/reflScale
    sbNew = spotNewDs.GetRasterBand(b).ReadAsArray()
    sbNew[~mask] = 0
    spotNewIm[:,:,b-1] = np.float32(sbNew)/reflScale
    nrb = ngiRawDs.GetRasterBand(b).ReadAsArray()
    nrb[~mask] = 0
    ngiRawIm[:,:,b-1] = np.float32(nrb)
    ncd = ngiCalibDs.GetRasterBand(b).ReadAsArray()
    ncd[~mask] = 0
    ngiCalibIm[:,:,b-1] = np.float32(ncd)/reflScale
    diffIm[:, :, b - 1] = (np.float32(sb)-np.float32(ncd)) / reflScale
    diffNewIm[:, :, b - 1] = (np.float32(sbNew) - np.float32(ncd)) / reflScale

# get ngi ranges to justify c=0
# pylab.figure()
# for b in range(1,4):
#     nrb = ngiRawDs.GetRasterBand(b).ReadAsArray()
#     nrb[~mask] = 0
#     ngiRawIm[:, :, b-1] = np.float32(nrb)
#     # find min vals for justifying C=0
#     print str.format("Band {0}: min {1} max {2} 0.01% {3} 99.99% {4} min% {5}", b, nrb[mask].min(), nrb[mask].max(),
#                      np.percentile(nrb[mask], 0.001), np.percentile(nrb[mask], 99.99),
#                      100.*np.percentile(nrb[mask], 0.001)/np.percentile(nrb[mask], 99.99))
#     pylab.subplot(2,2,b)
#     n, bins, patches = pylab.hist(nrb[mask], 200, range=(0,4096))
#
# pylab.figure()
# for b in range(1,4):
#     nrb = ngiCalibDs.GetRasterBand(b).ReadAsArray()
#     nrb[~mask] = 0
#     ngiRawIm[:, :, b-1] = np.float32(nrb)
#     # find min vals for justifying C=0
#     print str.format("Band {0}: min {1} max {2} 0.01% {3} 99.99% {4} min% {5}", b, nrb[mask].min(), nrb[mask].max(),
#                      np.percentile(nrb[mask], 0.001), np.percentile(nrb[mask], 99.99),
#                      100.*np.percentile(nrb[mask], 0.001)/np.percentile(nrb[mask], 99.99))
#     pylab.subplot(2,2,b)
#     n, bins, patches = pylab.hist(nrb[mask], 200, range=(0,4096))

spotDs = None
spotNewDs = None
ngiRawDs = None
ngiCalibDs = None

if False:
    pylab.figure()
    pylab.subplot(2, 2, 1)
    pylab.imshow(spotIm/spotIm.max())
    pylab.subplot(2, 2, 2)
    pylab.imshow(ngiRawIm/ngiRawIm.max())
    pylab.subplot(2, 2, 3)
    pylab.imshow(ngiCalibIm/ngiCalibIm.max())
    pylab.subplot(2, 2, 4)
    pylab.imshow(mask)

    pylab.figure()
    for i in range(0,3):
        if i == 0:
            ax = pylab.subplot(2, 3, 1)
        else:
            pylab.subplot(2, 3, 1 + i, sharex=ax, sharey=ax)
        pylab.imshow(diffIm[:,:,i], vmin=-0.15, vmax=.15)
        pylab.subplot(2, 3, 4 + i, sharex=ax, sharey=ax)
        pylab.imshow(diffNewIm[:,:,i], vmin=-0.15, vmax=.15)

# spotDs = gdal.Open(spotFileName, gdal.GA_ReadOnly)
# ngiDs = gdal.Open(ngiCalibFileName, gdal.GA_ReadOnly)
#
# spotIm = np.zeros((spotDs.RasterYSize, spotDs.RasterXSize, 3), dtype = np.float32)
# ngiIm = np.zeros_like(spotIm)
#
# for b in range(1,numBands):
#     spotIm[:,:,b-1] = spotDs.GetRasterBand(b).ReadAsArray()
#     ngiIm[:,:,b-1] = ngiDs.GetRasterBand(b).ReadAsArray()
#
# pylab.figure()
# pylab.subplot(1, 2, 1)
# pylab.imshow(spotIm/spotIm.max())
# pylab.subplot(1, 2, 2)
# pylab.imshow(ngiIm/ngiIm.max())
#
#
# mask = np.bool8(spotDs.GetRasterBand(1).GetMaskBand().ReadAsArray())
# mask2 = np.bool8(ngiDs.GetRasterBand(1).GetMaskBand().ReadAsArray())
# mask = mask&mask2
pylab.close('all')
bands = (0, 1, 2)
colors = ( 'orange', 'r', 'g')
legend = ('NIR', 'Red', 'Green')
step = 1000
fn = pylab.gcf().number
pn = 1
#
#
# fontSize = 24.
# mpl.rcParams.update({'font.size': fontSize})

f1 = pylab.figure('Uncalibrated')
f1.set_size_inches(16., 5.25, forward=True)
f2 = pylab.figure('Calibrated')
f2.set_size_inches(16., 5.25, forward=True)

## NB
spotIm = spotNewIm
ePixels = []
slopes = []
intercepts = []
for b in bands:
    sPixels = spotIm[:, :, b][mask]
    dRawPixels = ngiRawIm[:, :, b][mask]
    dCalibPixels = 0.96*ngiCalibIm[:, :, b][mask]
    e = sPixels -dCalibPixels
    ePixels.append(e)

    (slope, intercept, r, p, stde) = scipy.stats.linregress(sPixels.flatten(), dRawPixels.flatten())
    pylab.figure('Uncalibrated')
    pylab.subplot(1, 3, pn)
    pylab.plot(sPixels[::step], dRawPixels[::step], color='k', marker='.', linestyle='', markersize=.5)
    xl = pylab.gca().get_xlim()
    yl = pylab.gca().get_ylim()
    pylab.text((xl[0] + np.diff(xl)*0.05)[0], (yl[0] + np.diff(yl)*0.8)[0], str.format('$R^2$ = {0:.2f}',
                                                                                       np.round(r**2, 2)))
    # pylab.xlabel(r'SPOT5 $\rho_t$')
    pylab.xlabel(r'SPOT5 surface refl.')
    pylab.ylabel('DMC DN')
    pylab.title(legend[b])
    pylab.grid('on')
    pylab.tight_layout()
    pyplot.locator_params(axis='y', nbins=5)#to specify number of ticks on both or any single axes
    pyplot.locator_params(axis='x', nbins=6)

    #pylab.hold('on')

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
    pylab.title(legend[b])
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

f1.savefig('C:/Data/Development/Projects/MSc GeoInformatics/Docs/My Docs/Thesis/Retrieval of Surface Reflectance '
           'from Aerial Imagery/Figure 12 - DMC DN and SPOT5 Surface Reflectance Correlation.eps', dpi=1200)
f2.savefig('C:/Data/Development/Projects/MSc GeoInformatics/Docs/My Docs/Thesis/Retrieval of Surface Reflectance '
           'from Aerial Imagery/Figure 13 - DMC Homogenised and SPOT5 Surface Reflectance Correlation.eps', dpi=1200)

print str.format("N (num points plotted): {0}, step: {1}", sPixels[::step].__len__(), step)


# errors
# std dev(abs error), std dev(rms) etc dont really make sense and one measure of error variation is good enough i.e. std dev
# so i.e. use dont use the rhs () bracket vals but rather the third and last rows for error variation
ePixels = np.array(ePixels)*100
print str.format("Per band Abs error (%): {0} ({1})", np.abs(ePixels).mean(axis=1), np.abs(ePixels).std(axis=1))
print str.format("Per band RMS error (%): {0} ({1})", np.sqrt((ePixels**2).mean(axis=1)), np.sqrt((ePixels**2).std(axis=1)))  # does root-std-sqr make any sense?
print str.format("Per band Std Dev error (%): {0}", ePixels.std(axis=1))
print str.format("Mean Abs error (%): {0} ({1})", np.abs(ePixels).mean(), np.abs(ePixels).std())
print str.format("Mean RMS error (%): {0} ({1})", np.sqrt((ePixels**2).mean()), np.sqrt((ePixels**2).std()))  # does root-std-sqr make any sense?
print str.format("Std Dev error (%): {0}", ePixels.std())



#pylab.close('all')
colors = ('r','g','b','m')
step = 1000

pylab.figure()

for b in range(0, 3):
    sPixels = spotIm[:,:,b]
    sPixels = sPixels[mask].flatten()
    dCalibPixels = ngiCalibIm[:,:,b]
    dCalibPixels = dCalibPixels[mask].flatten()

    (slope, intercept, r, p, stde) = scipy.stats.linregress(sPixels, dCalibPixels)
    pylab.subplot(3, 1, b+1)
    pylab.plot(sPixels[::step], dCalibPixels[::step], colors[b-1] + 'x')
    pylab.title(str.format('r = {}', r**2))

#
## 
##
#make a scatter plot of SPOT vs MODIS
import numpy as np
import pylab
import os
import scipy.stats
import gdal
import scipy.ndimage
# from mayavi import mlab
import matplotlib.image
from matplotlib import pyplot
reflScale = 5000.

spotFileName = "D:\Data\Development\Projects\MSc GeoInformatics\Data\NGI\My Rectified\SpotInt16CsDs.tif"
modisFileName = "D:\Data\Development\Projects\MSc GeoInformatics\Data\NGI\Cross Calibration\Mosaics\ModisSpotExtent.tif"

spotDs = gdal.Open(spotFileName, gdal.GA_ReadOnly)
modisDs = gdal.Open(modisFileName, gdal.GA_ReadOnly)

spotIm = np.zeros((spotDs.RasterYSize, spotDs.RasterXSize, 3), dtype = np.float32)
modisIm = np.zeros((modisDs.RasterYSize, modisDs.RasterXSize, 3), dtype = np.float32)

modisMask = np.bool8(modisDs.GetRasterBand(1).GetMaskBand().ReadAsArray())

mask = np.bool8(spotDs.GetRasterBand(1).GetMaskBand().ReadAsArray())
mask2 = np.bool8(modisDs.GetRasterBand(1).GetMaskBand().ReadAsArray())
mask = mask & mask2

for b in range(1,4):
    sb = spotDs.GetRasterBand(b).ReadAsArray()
    sb[~mask] = 0
    spotIm[:,:,b-1] = np.float32(sb)/reflScale
    mb = modisDs.GetRasterBand(b).ReadAsArray()
    mb[~mask] = 0
    modisIm[:,:,b-1] = np.float32(mb)/reflScale


spotDs = None
modisDs = None

if False:
    pylab.figure()
    pylab.subplot(1, 2, 1)
    pylab.imshow(spotIm/spotIm.max())
    pylab.subplot(1, 2, 2)
    pylab.imshow(modisIm/modisIm.max())

# spotDs = gdal.Open(spotFileName, gdal.GA_ReadOnly)
# ngiDs = gdal.Open(ngiCalibFileName, gdal.GA_ReadOnly)
#
# spotIm = np.zeros((spotDs.RasterYSize, spotDs.RasterXSize, 3), dtype = np.float32)
# ngiIm = np.zeros_like(spotIm)
#
# for b in range(1,numBands):
#     spotIm[:,:,b-1] = spotDs.GetRasterBand(b).ReadAsArray()
#     ngiIm[:,:,b-1] = ngiDs.GetRasterBand(b).ReadAsArray()
#
# pylab.figure()
# pylab.subplot(1, 2, 1)
# pylab.imshow(spotIm/spotIm.max())
# pylab.subplot(1, 2, 2)
# pylab.imshow(ngiIm/ngiIm.max())
#
#
# mask = np.bool8(spotDs.GetRasterBand(1).GetMaskBand().ReadAsArray())
# mask2 = np.bool8(ngiDs.GetRasterBand(1).GetMaskBand().ReadAsArray())
# mask = mask&mask2
pylab.close('all')
bands = (0, 1, 2)
colors = ( 'orange', 'r', 'g')
legend = ('NIR', 'Red', 'Green')
step = 2
fn = pylab.gcf().number
pn = 1
#
#
# fontSize = 24.
# mpl.rcParams.update({'font.size': fontSize})

fontSize = 24.
import matplotlib as mpl
from matplotlib import pyplot
mpl.rcParams.update({'font.size': fontSize})

f1 = pylab.figure('MODIS vs SPOT')
f1.set_size_inches(16., 5.25, forward=True)

ePixels = []
spotModisSlope = []
spotModisIntercept = []

for b in bands:
    sPixels = spotIm[:, :, b][mask]
    mPixels = modisIm[:, :, b][mask]
    ePixels.append(sPixels - mPixels)


    (slope, intercept, r, p, stde) = scipy.stats.linregress(sPixels.flatten(), mPixels.flatten())
    spotModisSlope.append(slope)
    spotModisIntercept.append(intercept)
    pylab.figure('MODIS vs SPOT')
    pylab.subplot(1, 3, pn)
    pylab.plot(sPixels[::step], mPixels[::step], color='k', marker='.', linestyle='', markersize=.5)
    pylab.hold(True)
    m = 1. #np.min([np.max([np.percentile(sPixels[::step], 99.99), np.percentile(mPixels[::step], 99.99)]), 1.])
    oneLineH, = pylab.plot([0, m], [0, m], color='k', marker='', linestyle='--', label = '1:1')

    pylab.gca().set_xlim([0, m])
    pylab.gca().set_ylim([0, m])
    xl = pylab.gca().get_xlim()
    yl = pylab.gca().get_ylim()
    pylab.text((xl[0] + np.diff(xl)*0.05)[0], (yl[0] + np.diff(yl)*0.8)[0], str.format('$R^2$ = {0:.2f}',
                                                                                       np.round(r**2, 2)))
    # pylab.text((xl[0] + np.diff(xl)*0.05)[0], (yl[0] + np.diff(yl)*0.9)[0], str.format('y = {0:.2f} x + {1:.3f}',
    #                                                                                     slope, intercept))
    pylab.title(legend[b])
    # pylab.xlabel(r'SPOT5 $\rho_t$')
    # pylab.ylabel(r'DMC $\rho_t$')
    pylab.xlabel(r'SPOT5 surface refl.')
    pylab.ylabel(r'MODIS surface refl.')
    pylab.grid('on')
    pylab.legend(handles=[oneLineH], prop={'size':20}, loc=4)
    pylab.tight_layout()
    pyplot.locator_params(axis='y', nbins=5)#to specify number of ticks on both or any single axes
    pyplot.locator_params(axis='x', nbins=5)
    #pylab.hold('on')
    pn += 1

f1.savefig('C:/Data/Development/Projects/MSc GeoInformatics/Docs/My Docs/Thesis/Retrieval of Surface Reflectance '
           'from Aerial Imagery/Figure 14 - MODIS and SPOT5 Surface Reflectance Correlation.eps', dpi=1200)
# f2.savefig('C:/Data/Development/Projects/MSc GeoInformatics/Docs/My Docs/Thesis/Retrieval of Surface Reflectance from Aerial Imagery/dmc_surf_refl_spot_surf_refl_scatter_v2.png', dpi=300)

print str.format("N (num points plotted): {0}, step: {1}", sPixels[::step].__len__(), step)


# errors
# std dev(abs error), std dev(rms) etc dont really make sense and one measure of error variation is good enough i.e. std dev
# so i.e. use dont use the rhs () bracket vals but rather the third and last rows for error variation
ePixels = np.array(ePixels)*100
print str.format("Per band Abs error (%): {0} ({1})", np.abs(ePixels).mean(axis=1), np.abs(ePixels).std(axis=1))
print str.format("Per band RMS error (%): {0} ({1})", np.sqrt((ePixels**2).mean(axis=1)), np.sqrt((ePixels**2).std(axis=1)))  # does root-std-sqr make any sense?
print str.format("Per band Std Dev error (%): {0}", ePixels.std(axis=1))
print str.format("Mean Abs error (%): {0} ({1})", np.abs(ePixels).mean(), np.abs(ePixels).std())
print str.format("Mean RMS error (%): {0} ({1})", np.sqrt((ePixels**2).mean()), np.sqrt((ePixels**2).std()))  # does root-std-sqr make any sense?
print str.format("Std Dev error (%): {0}", ePixels.std())


## calculation of actual viewing geometry variations for the case study geometry
import numpy as np
altitude = 5000.
pixelWidth = 500.
dmcFov = 69.3 #42

#at nadir
dmcAngle = np.arctan(0.5*pixelWidth/altitude)*180./np.pi

#at extreme end of FOV
dmcExtent = 5000*np.tan(dmcFov*np.pi/(2*180.))
tmp = np.arctan((dmcExtent-pixelWidth)/altitude)*180./np.pi
dmcAngle = dmcFov/2. - tmp