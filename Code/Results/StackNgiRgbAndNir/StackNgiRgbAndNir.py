import gdal
import osr
import ogr
import pylab
import numpy as np

nirFile = "E:/NIR/3323D_2015_1001/NIR/3323d_2015_1001_01~0002_n.tif"
rgbFile = "E:/Unrectified_Aerials/3323D_2015_1001/3323D_2015_1001_01_0002_RGB.tif"

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