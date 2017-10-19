#hacked from gdal_merge to copy geotransform from one file to another
import math
import sys
import time

from osgeo import gdal

try:
    progress = gdal.TermProgress_nocb
except:
    progress = gdal.TermProgress

__version__ = '$id$'[5:-1]
verbose = 0
quiet = 0


def CopyGeotransform(srcFile, destFile):
    ds = gdal.OpenEx(srcFile, gdal.OF_RASTER)  # | gdal.GA_Update)
    if ds is None:
        print "Open failed."
        return
    srcGt = ds.GetGeoTransform()

    if False:
        print 'Driver: ', ds.GetDriver().ShortName,'/', \
              ds.GetDriver().LongName
        print 'Size is ',ds.RasterXSize,'x',ds.RasterYSize, \
              'x',ds.RasterCount
        print 'Projection is ',ds.GetProjection()
        if not srcGt is None:
            print 'Origin = (',srcGt[0], ',',srcGt[3],')'
            print 'Pixel Size = (',srcGt[1], ',',srcGt[5],')'
    ds = None

    ds = gdal.OpenEx(destFile, gdal.OF_RASTER | gdal.GA_Update)
    if ds is None:
        print "Open failed."
        return
    destGt = ds.GetGeoTransform()

    if False:
        print 'Driver: ', ds.GetDriver().ShortName,'/', \
              ds.GetDriver().LongName
        print 'Size is ',ds.RasterXSize,'x',ds.RasterYSize, \
              'x',ds.RasterCount
        print 'Projection is ',ds.GetProjection()
        if not destGt is None:
            print 'Origin = (',destGt[0], ',',destGt[3],')'
            print 'Pixel Size = (',destGt[1], ',',destGt[5],')'
        # ds.SetGeoTransform(srcGt)
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


# =============================================================================
def Usage():
    print('USE WITH CAUTION! - Will edit in place the geotransform of dest_file with that of src_file')
    print('Usage: copy_geotransform.py [src_file] [dest_file]')
    print('')

# =============================================================================
#
# Program mainline.
#

def main( argv=None ):

    global verbose, quiet
    verbose = 0
    quiet = 0
    names = []


    gdal.AllRegister()
    if argv is None:
        argv = sys.argv
    argv = gdal.GeneralCmdLineProcessor( argv )
    if argv is None:
        sys.exit( 0 )

    # Parse command line arguments.
    i = 1
    while i < len(argv):
        arg = argv[i]

        if arg[:1] == '-':
            print('Unrecognized command option: %s' % arg)
            Usage()
            sys.exit( 1 )
        else:
            names.append(arg)

        i = i + 1

    if len(names) == 0:
        print('No files selected.')
        Usage()
        sys.exit( 1 )

    if len(names) >= 2:
        print('Use only 2 files: src and dest.')
        Usage()
        sys.exit( 1 )

    print "Source file: ", names[0]
    print "Dest file: ", names[1]
    CopyGeotransform(names[0], names[1])

if __name__ == '__main__':
    sys.exit(main())
