import pylab
import os
import subprocess
import glob
# import gdal
# import ogr
# import rasterio
import numpy as np
from rasterio.transform import Affine

# NB setup env to run Anaconda GDAL
os.environ['PATH'] += "C:\ProgramData\Anaconda3\envs\py27\Library\\bin"
os.environ['GDAL_DATA'] = "C:\ProgramData\Anaconda3\envs\py27\Library\share\gdal"
os.environ['GDAL_DRIVER_PATH']="C:\ProgramData\Anaconda3\envs\py27\Library\\bin\gdalplugins"
# os.environ['PATH'] += "C:\OSGeo4W64\\bin"
# os.environ['GDAL_DATA'] = "C:\OSGeo4W64\share\gdal"
# os.environ['GDAL_DRIVER_PATH']="C:\OSGeo4W64\\bin\gdalplugins"


gdalwarpExe = "C:\ProgramData\Anaconda3\envs\py27\Library\\bin\gdalwarp.exe"  # Use anaconda version as osgeo4w wont run under Anaconda
gdaltranslateExe = r"C:\OSGeo4W64\bin\gdal_translate.exe"  # Use anaconda version as osgeo4w wont run under Anaconda
gdalmergeExe = "C:\ProgramData\Anaconda3\envs\py27\Scripts\gdal_merge.py"  # Use anaconda version as osgeo4w wont run under Anaconda
gdaladdoExe = "C:\ProgramData\Anaconda3\envs\py27\Library\\bin\gdaladdo.exe"  # Use anaconda version as osgeo4w wont run under Anaconda

os.environ['PATH'] += "C:\Data\Development\Toolboxes\OpenCV-2.4.8\\build\\x64\\vc12\\bin"

# xcalibExe = "C:\Data\Development\Projects\PhD GeoInformatics\Code\Cross Calibration\\x64\Release\CrossCalibration.exe"
# modisFileName = "D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.A2016298.h19v12.006.2016307063215\MCD43A4.A2016302.h19v12.006.2016311094217.NgiBandOrder.TmLo19.tif"
imDir = r'V:\Data\NGI\UnRectified\3323D_2015_1001\Photoscan'

# make the working dir on the ssd drive, this speeds things up significantly
# os.chdir("C:\Data\Development\Projects\PhD GeoInformatics\Code\Results\NGI Radiometric Homogenisation")
for imFileName in glob.glob(os.path.join(imDir, '*RGBN.tif')):
    baseFileName = os.path.split(imFileName)[-1]
    outFileName = baseFileName[:-4] + '_DCMP.tif'
    print 'Recompressing {0}'.format(os.path.join(imDir, baseFileName))
    cmdString = '"{0}" -co COMPRESS=DEFLATE -co NBITS=16 "{1}" "{2}"'.format(gdaltranslateExe,
        os.path.join(imDir, baseFileName), os.path.join(imDir, outFileName))
    print cmdString
    subprocess.call(cmdString, shell=True, env=os.environ)

# os.chdir("C:\Data\Development\Projects\PhD GeoInformatics\Code\Results\NGI Radiometric Homogenisation")
for imFileName in glob.glob(os.path.join(imDir, '*RGBN_DCMP.tif')):
    baseFileName = os.path.split(imFileName)[-1]
    print 'Adding overviews for {0}'.format(os.path.join(imDir, baseFileName))
    cmdString = '"{0}" -ro -r average --config COMPRESS_OVERVIEW DEFLATE "{1}" 2 4 8 16 32 64'.format(gdaladdoExe,
        os.path.join(imDir, baseFileName))
    subprocess.call(cmdString, shell=True, env=os.environ)


#
# for file in glob.glob(os.path.join(orthoDir, 'o*_RGB_XCALIB.tif')):
#     baseFileName = os.path.split(file)[-1]  # 20 or 60m file
#     outFileName = baseFileName + '.ovr'
#     if os.path.exists(os.path.join(orthoDir, outFileName)):
#         print '{0} exists'.format(os.path.join(orthoDir, outFileName))
#     else:
#         print 'Generating overviews for {0}'.format(os.path.join(orthoDir, outFileName))
#         cmdString = '"{0}" -ro -r average --config COMPRESS_OVERVIEW DEFLATE "{1}" 2 4 6 8 16 32 64'.format(gdaladdoExe, file)
#         subprocess.call(cmdString, shell=True, env=os.environ)
#
