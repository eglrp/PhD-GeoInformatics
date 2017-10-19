#setting of nodata in nir
#arcmap xres,yres differs to qgis
#arcmap has nodata set but qgis not ?
#does a ok image also have geoxform ul cnr excl nodata
#the rasters themselves seem to have the right cnr co-ords i.e. the ones excl nodata so we dont need shapefiles ?
import gdal
import osr
import ogr
import pylab
import numpy as np

nirFile = "E:/NIR/3323D_2015_1001/NIR/3323d_2015_1001_01~0002_n - scratch.tif"
rgbFile = "E:/Unrectified_Aerials/3323D_2015_1001/3323D_2015_1001_01_0002_RGB.tif"
combFile = "E:/NIR/3323D_2015_1001/NIR/RGBNtest.tif"

def UpdateGeotransform(fileName, transform):
    ds = gdal.OpenEx(fileName, gdal.OF_RASTER)  # | gdal.GA_Update)
    if ds is None:
        print "Open failed."
        return

    print 'Driver: ', ds.GetDriver().ShortName,'/', \
          ds.GetDriver().LongName
    print 'Size is ',ds.RasterXSize,'x',ds.RasterYSize, \
          'x',ds.RasterCount
    print 'Projection is ',ds.GetProjection()
    geotransform = ds.GetGeoTransform()
    if not geotransform is None:
        print 'Origin = (',geotransform[0], ',',geotransform[3],')'
        print 'Pixel Size = (',geotransform[1], ',',geotransform[5],')'
    # if lry is not None:
    #     gt = [ ulx, (lrx - ulx) / ds.RasterXSize, 0,
    #            uly, 0, (lry - uly) / ds.RasterYSize ]
    #     ds.SetGeoTransform(gt)
    #
    # if yres is not None:
    #     gt = ds.GetGeoTransform()
    #     # Doh ! why is gt a tuple and not an array...
    #     gt = [ gt[j] for j in range(6) ]
    #     gt[1] = xres
    #     gt[5] = yres
    #     ds.SetGeoTransform(gt)
    ulx = geotransform[0]
    uly = geotransform[3]
    lrx = ulx + np.float64(ds.RasterXSize) * geotransform[1]
    lry = uly + np.float64(ds.RasterYSize) * geotransform[5]
    print type(geotransform[0])
    print geotransform
    print "UL: (", ulx, ", ", uly, ")"
    print "LR: (", lrx, ", ", lry, ")"
    pts = np.array([[ulx, uly], [lrx, lry]])
    transPts = np.dot(np.column_stack((pts, np.ones((2,1)))), transform)
    print transPts

    pylab.figure()
    pylab.plot(pts[:, 0], pts[:, 1], 'b')
    pylab.plot(transPts[:, 0], transPts[:, 1], 'r')
    pylab.axis('square')

    gt = [transPts[0,0], (transPts[1,0] - transPts[0,0]) / ds.RasterXSize, 0,
          transPts[0, 1], 0, (transPts[1, 1] - transPts[0, 1]) / ds.RasterYSize]
    # ds.SetGeoTransform(gt)

    ds = None
    return pts, transPts


def CopyGeotransform(srcFile, destFile):
    ds = gdal.OpenEx(srcFile, gdal.OF_RASTER)  # | gdal.GA_Update)
    if ds is None:
        print "Open failed."
        return

    print 'Driver: ', ds.GetDriver().ShortName,'/', \
          ds.GetDriver().LongName
    print 'Size is ',ds.RasterXSize,'x',ds.RasterYSize, \
          'x',ds.RasterCount
    print 'Projection is ',ds.GetProjection()
    srcGt = ds.GetGeoTransform()
    if not srcGt is None:
        print 'Origin = (',srcGt[0], ',',srcGt[3],')'
        print 'Pixel Size = (',srcGt[1], ',',srcGt[5],')'
    ds = None

    ds = gdal.OpenEx(destFile, gdal.OF_RASTER | gdal.GA_Update)
    if ds is None:
        print "Open failed."
        return

    print 'Driver: ', ds.GetDriver().ShortName,'/', \
          ds.GetDriver().LongName
    print 'Size is ',ds.RasterXSize,'x',ds.RasterYSize, \
          'x',ds.RasterCount
    print 'Projection is ',ds.GetProjection()
    destGt = ds.GetGeoTransform()
    if not destGt is None:
        print 'Origin = (',destGt[0], ',',destGt[3],')'
        print 'Pixel Size = (',destGt[1], ',',destGt[5],')'
    ds.SetGeoTransform(srcGt)
    ds = None


    # if lry is not None:
    #     gt = [ ulx, (lrx - ulx) / ds.RasterXSize, 0,
    #            uly, 0, (lry - uly) / ds.RasterYSize ]
    #     ds.SetGeoTransform(gt)
    #
    # if yres is not None:
    #     gt = ds.GetGeoTransform()
    #     # Doh ! why is gt a tuple and not an array...
    #     gt = [ gt[j] for j in range(6) ]
    #     gt[1] = xres
    #     gt[5] = yres
    #     ds.SetGeoTransform(gt)
    return srcGt, destGt

def GetCornerCoords(fileName):
    ds = gdal.OpenEx(fileName, gdal.OF_VECTOR)
    if ds is None:
        print "Open failed"
        return

    lyr = ds.GetLayerByIndex(0)
    lyr.ResetReading()
    nir_spatial_ref = lyr.GetSpatialRef()

    pts = np.zeros((4,2))
    for (i, feat) in enumerate(lyr):
        print '.',
        # feat_defn = lyr.GetLayerDefn()
        # f = {}
        # for i in range(feat_defn.GetFieldCount()):
        #     field_defn = feat_defn.GetFieldDefn(i)
        #     f[field_defn.GetName()] = feat.GetField(i)
        poly = feat.GetGeometryRef()
        if poly is not None and (poly.GetGeometryType() == ogr.wkbPolygon):
            boundary = poly.GetBoundary()
            print boundary.GetPointCount()
            for ip in xrange(boundary.GetPointCount()-1):
                point = boundary.GetPoint(ip)
                print "%.6f, %.6f" % (point[0], point[1])
                pts[ip,:] =[point[0], point[1]]
        else:
            print "no point geometry/n"
        # gcpList.append(f)
    print ' '

    ds = None
    return pts

nirPts = GetCornerCoords(nirFile + ".shp")
rgbPts = GetCornerCoords(rgbFile + ".shp")

pylab.figure()
pylab.plot(nirPts[:,0], nirPts[:,1],'b')
pylab.plot(rgbPts[:,0], rgbPts[:,1],'r')
pylab.axis('square')

# nirPts * trans = rgbPts
# n x 3 * 3 x 2 = n x 2
# linalg.lstsq(a,b) : a x = b for x
trans, res, rank, s = np.linalg.lstsq(np.column_stack((nirPts, np.ones((4,1)))), rgbPts)
np.dot(np.column_stack((nirPts, np.ones((4,1)))), trans) - rgbPts

pts, transPts = UpdateGeotransform(nirFile, trans)

rgbPts2, dum = UpdateGeotransform(rgbFile, trans)

pylab.figure()
pylab.plot(nirPts[:,0], nirPts[:,1], 'b')
pylab.plot(rgbPts[:,0], rgbPts[:,1], 'r')
pylab.plot(rgbPts2[:,0], rgbPts2[:,1], 'g-x')
# pylab.plot(pts[:,0], pts[:,1], 'g-x')
# pylab.plot(transPts[:,0], transPts[:,1], 'k-x')
pylab.axis('square')

#
# # construct intermediate matrix
# Q = p[1:] - p[0]
# Q_prime = p_prime[1:] - p_prime[0]
#
# # calculate rotation matrix
# R = np.dot(np.linalg.inv(np.row_stack((Q, np.cross(*Q)))),
#            np.row_stack((Q_prime, np.cross(*Q_prime))))
#
# # calculate translation vector
# t = p_prime[0] - np.dot(p[0], R)
#
# # calculate affine transformation matrix
# return np.column_stack((np.row_stack((R, t)),
#                         (0, 0, 0, 1)))
# trans t = p_prime[0] - np.dot(p[0], R)
#
#     # calculate affine transformation matrix
#     return np.column_stack((np.row_stack((R, t)),
#                             (0, 0, 0, 1)))

rgbGt, nirGt = CopyGeotransform(rgbFile, nirFile)


rgbGt, combGt = CopyGeotransform(rgbFile, combFile)