import numpy as np
import pylab
import os
import gdal
import subprocess

sourceDir = "D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\XCalib Experiments\Calibrated2\w11p1"
gdaltindexExe = "C:\ProgramData\Anaconda3\envs\py27\Library\\bin\gdaltindex.exe"

print 'Processing {0:%s}' % (sourceDir)
tileIndexFileName = os.path.join(sourceDir, os.path.basename(sourceDir) + '_TileIndex.shp')
subprocess.call('{0} "{1}" "{2}"'.format(gdaltindexExe, tileIndexFileName,
                                           os.path.join(sourceDir, "*_CMP_XCALIB.tif")), shell=True)
