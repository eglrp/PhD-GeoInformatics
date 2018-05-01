import gdal
import numpy as np
import os
import pylab
import scipy.stats
from matplotlib import pyplot
import matplotlib as mpl

spotFileName = "D:\Data\Development\Projects\PhD GeoInformatics\Data\SPOT\S131022114824832\Orthorectification\oATCORCorrected_METADATA_00812098_AutoGCP_NgiFormat.tif"  # improved orthorect
modisFileName = 'D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.A2010025.h19v12.005.2010043064233.Lo23.RGBN.tif'
calibRootDir = 'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\XCalib Experiments\Calibrated2\\'

mosaicFileName = "D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\XCalib Experiments\Calibrated\w11p1\w33p1_Mosaic10m.tif"
mosaicFileName = "D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\XCalib Experiments\Calibrated\w44p2\w44p2_Mosaic10m.tif"
reflScale = 5000.
plotFigures = True

spotDs = gdal.Open(spotFileName, gdal.GA_ReadOnly)

spotIm = np.zeros((spotDs.RasterYSize, spotDs.RasterXSize, 3), dtype=np.float32)

spotMask = np.bool8(spotDs.GetRasterBand(1).GetMaskBand().ReadAsArray())
for b in range(1, 4):
    sb = spotDs.GetRasterBand(b).ReadAsArray()
    sb[~spotMask] = 0
    spotIm[:, :, b - 1] = np.float32(sb) / reflScale

spotGeoTransform = spotDs.GetGeoTransform()
spotInvGeoTransform = gdal.InvGeoTransform(spotGeoTransform)

spotDs = None

# loop though subdirs and make a dict of what is there
results = {}
for subDir in os.listdir(calibRootDir):
    if os.path.isdir(os.path.join(calibRootDir, subDir)):
        print subDir
        res = {}
        res['name'] = subDir
        res['winSize'] = (np.int(subDir[1]), np.int(subDir[2]))
        res['winArea'] = np.int(subDir[1])*np.int(subDir[2])
        res['model'] = np.int(subDir[4])
        res['ngiCalibFileName'] = os.path.join(calibRootDir, subDir, subDir + '_Mosaic10m.tif')
        results[subDir] = res


models = np.array([r['model'] for r in results.values()])
winSizes = np.array([r['winSize'] for r in results.values()])
winAreas = np.array([r['winArea'] for r in results.values()])

spotBands = (0, 1, 2)
ngiBands = (3, 0, 1)


fontSize = 12.
mpl.rcParams.update({'font.size': fontSize})

# loop though dict and visualise results by model and win size
for mi, model in enumerate(np.unique(models)):
    if plotFigures:
        pylab.close('Model ' + str(model))
        f = pylab.figure('Model ' + str(model))
        f.set_size_inches(16., 9., forward=True)

    modelWinAreas = np.sort(winAreas[models == model])
    #winAreasUnique = np.unique(winAreas)
    for wi, modelWinArea in enumerate(modelWinAreas):
        idx =  np.argwhere(np.logical_and(models == model, winAreas == modelWinArea))[0]
        if idx.size > 1:
            raise "len(idx = models == model & winAreas == winArea) > 1"
        res = results.values()[idx[0]]
        print res['name']
        print 'Opening:' + res['ngiCalibFileName']
        ngiCalibDs = gdal.Open(res['ngiCalibFileName'], gdal.GA_ReadOnly)
        ngiCalibIm = np.zeros((ngiCalibDs.RasterYSize, ngiCalibDs.RasterXSize, 4), dtype=np.float32)

        ngiCalibMask = np.bool8(ngiCalibDs.GetRasterBand(1).GetMaskBand().ReadAsArray())
        for b in range(1, 5):
            ncd = ngiCalibDs.GetRasterBand(b).ReadAsArray()
            ngiCalibMask = (ngiCalibMask & (ncd > 0))

        for b in range(1, 5):
            ncd = ngiCalibDs.GetRasterBand(b).ReadAsArray()
            ncd[~ngiCalibMask] = 0
            ngiCalibIm[:, :, b - 1] = np.float32(ncd) / reflScale

        # find ngi extents in spot
        ngiCalibGeoTransform = ngiCalibDs.GetGeoTransform()

        (ulx, uly) = gdal.ApplyGeoTransform(ngiCalibGeoTransform, 0, 0)
        (brx, bry) = gdal.ApplyGeoTransform(ngiCalibGeoTransform, ngiCalibDs.RasterXSize, ngiCalibDs.RasterYSize)
        (refUlx, refUly) = np.int32(gdal.ApplyGeoTransform(spotInvGeoTransform, ulx, uly))
        (refBrx, refBry) = np.int32(gdal.ApplyGeoTransform(spotInvGeoTransform, brx, bry))

        spotSubIm = spotIm[refUly:refBry, refUlx:refBrx, :]
        spotSubMask = spotMask[refUly:refBry, refUlx:refBrx]
        mask = ngiCalibMask & spotSubMask

        ngiCalibDs = None

        ePixels = []
        slopes = []
        intercepts = []
        rs = []
        step = 1000
        diffIm = np.zeros_like(spotSubIm)
        colors = ('orange', 'r', 'g')
        legend = ('NIR', 'Red', 'Green')
        for (bi, (sb, nb)) in enumerate(zip(spotBands, ngiBands)):
            diffIm[:, :, sb] =  spotSubIm[:, :, sb] - ngiCalibIm[:, :, nb]
            diffIm[:, :, sb][~mask] =  0
            sPixels = spotSubIm[:, :, sb][mask]
            dCalibPixels = ngiCalibIm[:, :, nb][mask]
            e = sPixels - dCalibPixels  #  - (sPixels.min() - dCalibPixels.min())   # naughty hack
            # e = e-np.mean(e)
            ePixels.append(e)

            (slope, intercept, r, p, stde) = scipy.stats.linregress(sPixels.flatten(), dCalibPixels.flatten())
            slopes.append(slope)
            intercepts.append(intercept)
            rs.append(r)

            if plotFigures:
                pylab.subplot(3, modelWinAreas.size, sb*modelWinAreas.size + wi + 1)
                pylab.plot(sPixels[::step], dCalibPixels[::step], color='k', marker='.', linestyle='', markersize=.5)
                pylab.hold(True)

                pylab.plot([0., 1.], [intercept, slope+intercept], 'g:')
                m = 1.  # np.min([np.max([np.percentile(sPixels[::step], 99.99), np.percentile(dCalibPixels[::step], 99.99)]), 1.])
                oneLineH, = pylab.plot([0, m], [0, m], color='k', marker='', linestyle='--', label='1:1')
                pylab.gca().set_xlim([0, m])
                pylab.gca().set_ylim([0, m])
                xl = pylab.gca().get_xlim()
                yl = pylab.gca().get_ylim()
                pylab.text((xl[0] + np.diff(xl) * 0.05)[0], (yl[0] + np.diff(yl) * 0.8)[0], str.format('$R^2$ = {0:.2f}',
                                                                                                       np.round(r ** 2, 3)))
                # pylab.text((xl[0] + np.diff(xl)*0.05)[0], (yl[0] + np.diff(yl)*0.9)[0], str.format('y = {0:.2f} x + {1:.3f}',
                #                                                                                    slope, intercept))
                pylab.title(legend[sb] + ' - Win: ({0}x{1})'.format(*res['winSize']))
                # pylab.xlabel(r'SPOT5 $\rho_t$')
                # pylab.ylabel(r'DMC $\rho_t$')
                pylab.xlabel(r'SPOT5 surface refl.')
                pylab.ylabel(r'DMC surface refl.')
                pylab.grid('on')
                # pylab.legend(handles=[oneLineH], prop={'size': 20}, loc=4)
                # pylab.tight_layout()
                pyplot.locator_params(axis='y', nbins=5)  # to specify number of ticks on both or any single axes
                pyplot.locator_params(axis='x', nbins=5)
        # pylab.hold('on')
        # pn += 1

        ePixels = np.array(np.float64(ePixels)) * 100
        res['slopes'] = np.array(slopes)
        res['intercepts'] = np.array(intercepts)
        res['rs'] = np.array(rs)
        res['maeb'] = np.abs(ePixels).mean(axis=1)
        res['mae'] = np.abs(ePixels).mean()
        res['rmsb'] = np.sqrt((ePixels ** 2).mean(axis=1))
        res['rms'] = np.sqrt((ePixels ** 2).mean())
        res['rs2'] = np.array(rs)**2
        res['mean(rs2)'] = np.mean(res['rs2'])
        # if plotFigures:
        #     res['dPixels'] = np.float32(dPixels)
        #     res['diffIm'] = np.float32(diffIm)
        results[res['name']] = res  # update results

    if plotFigures:
        pylab.tight_layout()


maes = np.array([r['mae'] for r in results.values()])
r2s = np.array([r['mean(rs2)'] for r in results.values()])
rmss = np.array([r['rms'] for r in results.values()])

pylab.figure()
for model in np.unique(models):
    modelIdx = models == model
    modelWinAreas = winAreas[modelIdx]
    sortIdx = np.argsort(modelWinAreas)
    pylab.subplot(3, 1, 1)
    pylab.plot(modelWinAreas[sortIdx], (maes[modelIdx])[sortIdx], 'x-')
    pylab.subplot(3, 1, 2)
    pylab.plot(modelWinAreas[sortIdx], (r2s[modelIdx])[sortIdx], 'x-')
    pylab.subplot(3, 1, 3)
    pylab.plot(modelWinAreas[sortIdx], (rmss[modelIdx])[sortIdx], 'x-')

pylab.subplot(3, 1, 1)
pylab.xlabel('Win. area (pixels)')
pylab.ylabel('Mean Abs. error (%)')
pylab.title('Effect of sliding window size on SPOT comparison')
pylab.grid()
pylab.legend(np.unique(models))
# pylab.legend('Gain only model')
pylab.subplot(3, 1, 2)
pylab.xlabel('Win. area (pixels)')
pylab.ylabel('$R^2$')
pylab.grid()
# pylab.title('$R^2$ vs window size')
pylab.subplot(3, 1, 3)
pylab.xlabel('Win. area (pixels)')
pylab.ylabel('RMS error (%)')
pylab.grid()
# pylab.title('RMS error vs window size')
pylab.tight_layout()



# for avn meet:
pylab.figure()
for model in [1]:
    modelIdx = models == model
    modelWinAreas = winAreas[modelIdx]
    sortIdx = np.argsort(modelWinAreas)
    pylab.subplot(3, 1, 1)
    pylab.plot(modelWinAreas[sortIdx], (maes[modelIdx])[sortIdx], 'x-')
    pylab.subplot(3, 1, 2)
    pylab.plot(modelWinAreas[sortIdx], (r2s[modelIdx])[sortIdx], 'x-')
    pylab.subplot(3, 1, 3)
    pylab.plot(modelWinAreas[sortIdx], (rmss[modelIdx])[sortIdx], 'x-')

pylab.subplot(3, 1, 1)
pylab.xlabel('Win. area (pixels)')
pylab.ylabel('Mean Abs. error (%)')
pylab.title('Effect of sliding window size on SPOT comparison (gain only model)')
pylab.grid()
# pylab.legend('Gain only model')
pylab.subplot(3, 1, 2)
pylab.xlabel('Win. area (pixels)')
pylab.ylabel('$R^2$')
pylab.grid()
# pylab.title('$R^2$ vs window size')
pylab.subplot(3, 1, 3)
pylab.xlabel('Win. area (pixels)')
pylab.ylabel('RMS error (%)')
pylab.grid()
# pylab.title('RMS error vs window size')
pylab.tight_layout()

