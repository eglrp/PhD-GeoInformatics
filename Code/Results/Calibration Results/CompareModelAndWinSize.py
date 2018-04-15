import gdal
import numpy as np
import os
import pylab
import scipy.stats
from matplotlib import pyplot

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

# loop though subdirs
results = {}
for subDir in os.listdir(calibRootDir):
    if os.path.isdir(os.path.join(calibRootDir, subDir)):
        print subDir
        mosaicFileName = os.path.join(calibRootDir, subDir, subDir + '_Mosaic10m.tif')
        ngiCalibDs = gdal.Open(mosaicFileName, gdal.GA_ReadOnly)
        ngiCalibIm = np.zeros((ngiCalibDs.RasterYSize, ngiCalibDs.RasterXSize, 4), dtype=np.float32)

        ngiCalibMask = np.bool8(ngiCalibDs.GetRasterBand(1).GetMaskBand().ReadAsArray())
        for b in range(1, 5):
            ncd = ngiCalibDs.GetRasterBand(b).ReadAsArray()
            ngiCalibMask = (ngiCalibMask & (ncd > 0))

        if False:
            pylab.figure()
            ax = pylab.subplot(1, 2, 1)
            pylab.imshow(spotMask)
            pylab.subplot(1, 2, 2)
            pylab.imshow(ngiCalibMask)

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
        res = {}
        res['name'] = subDir
        res['winSize'] = (np.int(subDir[1]), np.int(subDir[2]))
        res['winArea'] = np.int(subDir[1])*np.int(subDir[2])
        res['model'] = np.int(subDir[4])
        res['ngiCalibFileName'] = mosaicFileName

        spotBands = (0, 1, 2)
        ngiBands = (3, 0, 1)
        step = 1000
        pn = 1
        if plotFigures:
            colors = ('orange', 'r', 'g')
            legend = ('NIR', 'Red', 'Green')
            # pylab.close('all')
            # fn = pylab.gcf().number
            # # fontSize = 24.
            # # mpl.rcParams.update({'font.size': fontSize})
            # f2 = pylab.figure('Calibrated')
            # f2.set_size_inches(16., 5.25, forward=True)

        ## NB
        ePixels = []
        slopes = []
        intercepts = []
        rs = []
        diffIm = np.zeros_like(spotSubIm)
        for sb, nb in zip(spotBands, ngiBands):
            # diffIm[:, :, sb] =  spotSubIm[:, :, sb] - ngiCalibIm[:, :, nb]
            # diffIm[:, :, ~mask] =  0
            sPixels = spotSubIm[:, :, sb][mask]
            dCalibPixels = ngiCalibIm[:, :, nb][mask]
            e = sPixels - dCalibPixels
            ePixels.append(e)

            (slope, intercept, r, p, stde) = scipy.stats.linregress(sPixels.flatten(), dCalibPixels.flatten())
            slopes.append(slope)
            intercepts.append(intercept)
            rs.append(r)
            if plotFigures:
                pylab.figure('Model - ' + str(res['model']))
                pylab.subplot(3, 5, 1+pn*5)
                pylab.plot(sPixels[::step], dCalibPixels[::step], color='k', marker='.', linestyle='', markersize=.5)
                pylab.hold(True)
                m = 1.  # np.min([np.max([np.percentile(sPixels[::step], 99.99), np.percentile(dCalibPixels[::step], 99.99)]), 1.])
                oneLineH, = pylab.plot([0, m], [0, m], color='k', marker='', linestyle='--', label='1:1')
                pylab.gca().set_xlim([0, m])
                pylab.gca().set_ylim([0, m])
                xl = pylab.gca().get_xlim()
                yl = pylab.gca().get_ylim()
                pylab.text((xl[0] + np.diff(xl) * 0.05)[0], (yl[0] + np.diff(yl) * 0.8)[0], str.format('$R^2$ = {0:.2f}',
                                                                                                       np.round(r ** 2, 2)))
                # pylab.text((xl[0] + np.diff(xl)*0.05)[0], (yl[0] + np.diff(yl)*0.9)[0], str.format('y = {0:.2f} x + {1:.3f}',
                #                                                                                    slope, intercept))
                pylab.title(legend[sb])
                # pylab.xlabel(r'SPOT5 $\rho_t$')
                # pylab.ylabel(r'DMC $\rho_t$')
                pylab.xlabel(r'SPOT5 surface refl.')
                pylab.ylabel(r'DMC surface refl.')
                pylab.grid('on')
                pylab.legend(handles=[oneLineH], prop={'size': 20}, loc=4)
                pylab.tight_layout()
                pyplot.locator_params(axis='y', nbins=5)  # to specify number of ticks on both or any single axes
                pyplot.locator_params(axis='x', nbins=5)
            # pylab.hold('on')
            pn += 1

        #store results in dict
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
        results[subDir] = res

        #
        # f1.savefig('C:/Data/Development/Projects/PhD GeoInformatics/Docs/My Docs/Thesis/Retrieval of Surface Reflectance '
        #            'from Aerial Imagery/Figure 12 - DMC DN and SPOT5 Surface Reflectance Correlation V2.eps', dpi=1200)
        # f2.savefig('C:/Data/Development/Projects/PhD GeoInformatics/Docs/My Docs/Thesis/Retrieval of Surface Reflectance '
        #            'from Aerial Imagery/Figure 13 - DMC Homogenised and SPOT5 Surface Reflectance Correlation V2.eps', dpi=1200)

        # print str.format("N (num points plotted): {0}, step: {1}", sPixels[::step].__len__(), step)
        #
        # # errors
        # # std dev(abs error), std dev(rms) etc dont really make sense and one measure of error variation is good enough i.e. std dev
        # # so i.e. use dont use the rhs () bracket vals but rather the third and last rows for error variation
        # ePixels = np.array(np.float64(ePixels)) * 100
        # print str.format("Per band Abs error (%): {0} ({1})", np.abs(ePixels).mean(axis=1).round(decimals=2),
        #                  np.abs(ePixels).std(axis=1).round(decimals=2))
        # print str.format("Per band RMS error (%): {0} ({1})", np.sqrt((ePixels ** 2).mean(axis=1)).round(decimals=2),
        #                  np.sqrt((ePixels ** 2).std(axis=1)).round(decimals=2))  # does root-std-sqr make any sense?
        # print str.format("Per band Std Dev error (%): {0}", ePixels.std(axis=1).round(decimals=2))
        # print str.format("Mean Abs error (%): {0:.2f} ({1:.2f})", np.abs(ePixels).mean(), np.abs(ePixels).std())
        # print str.format("Mean RMS error (%): {0:.2f} ({1:.2f})", np.sqrt((ePixels ** 2).mean()),
        #                  np.sqrt((ePixels ** 2).std()))  # does root-std-sqr make any sense?
        # print str.format("Std Dev error (%): {0:.2f}", ePixels.std())
models = np.array([r['model'] for r in results.values()])
winSizes = np.array([r['winSize'] for r in results.values()])
winAreas = np.array([r['winArea'] for r in results.values()])
maes = np.array([r['mae'] for r in results.values()])
r2s = np.array([r['mean(rs2)'] for r in results.values()])

pylab.figure()
for model in np.unique(models):
    modelIdx = models == model
    modelWinAreas = winAreas[modelIdx]
    sortIdx = np.argsort(modelWinAreas)
    pylab.subplot(2, 1, 1)
    pylab.plot(modelWinAreas[sortIdx], (maes[modelIdx])[sortIdx], 'x-')
    pylab.subplot(2, 1, 2)
    pylab.plot(modelWinAreas[sortIdx], (r2s[modelIdx])[sortIdx], 'x-')

pylab.subplot(2, 1, 1)
pylab.xlabel('Win area')
pylab.ylabel('MAE')
pylab.grid()
pylab.legend(np.unique(models))
pylab.subplot(2, 1, 2)
pylab.xlabel('Win area')
pylab.ylabel('R2')
pylab.grid()
pylab.legend(np.unique(models))

spotBands = (0, 1, 2)
ngiBands = (3, 0, 1)

if plotFigures:  # detailed scatter and diff im plots
    import matplotlib as mpl
    fontSize = 12.
    mpl.rcParams.update({'font.size': fontSize})
    for mi, model in enumerate(np.unique(models)):
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
                e = sPixels - dCalibPixels
                ePixels.append(e)

                (slope, intercept, r, p, stde) = scipy.stats.linregress(sPixels.flatten(), dCalibPixels.flatten())
                slopes.append(slope)
                intercepts.append(intercept)
                rs.append(r)

                pylab.subplot(3, modelWinAreas.size, sb*modelWinAreas.size + wi + 1)
                pylab.plot(sPixels[::step], dCalibPixels[::step], color='k', marker='.', linestyle='', markersize=.5)
                pylab.hold(True)

                pylab.plot([0., 1.], [intercept, slope+intercept], 'g-')
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
            pn += 1

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
            results[res['name']] = res

        pylab.tight_layout()


