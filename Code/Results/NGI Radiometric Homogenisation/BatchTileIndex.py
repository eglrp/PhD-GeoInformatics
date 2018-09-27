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
gdaltranslateExe = "C:\ProgramData\Anaconda3\envs\py27\Library\\bin\gdal_translate.exe"  # Use anaconda version as osgeo4w wont run under Anaconda
gdalmergeExe = "C:\ProgramData\Anaconda3\envs\py27\Scripts\gdal_merge.py"  # Use anaconda version as osgeo4w wont run under Anaconda
gdaltindexExe = "C:\ProgramData\Anaconda3\envs\py27\Library\\bin\gdaltindex.exe"  # Use anaconda version as osgeo4w wont run under Anaconda

os.environ['PATH'] += "C:\Data\Development\Toolboxes\OpenCV-2.4.8\\build\\x64\\vc12\\bin"

xcalibExe = "C:\Data\Development\Projects\PhD GeoInformatics\Code\Cross Calibration\\x64\Release\CrossCalibration.exe"
# modisFileName = "D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.A2016298.h19v12.006.2016307063215\MCD43A4.A2016302.h19v12.006.2016311094217.NgiBandOrder.TmLo19.tif"
batchDir = r"E:\Homogenised\3318B_2016_1142"
batchDir = r"V:\Data\NGI\Rectified\3324C_2015_1004\RGBN"

tileIndexFileName = r"E:\Homogenised\3318B_2016_1142\3318B_2016_1142_RGB_XCALIB_TileIndex.shp"
tileIndexFileName = r"E:\Homogenised\3324C_2015_1004\RGBN\3324C_2015_1004_RGBN_TileIndex.shp"
# make the working dir on the ssd drive, this speeds things up significantly
os.chdir("C:\Data\Development\Projects\PhD GeoInformatics\Code\Results\NGI Radiometric Homogenisation")
for tifFileName in glob.glob(os.path.join(batchDir, '*RGBN.tif')):
    baseFileName = os.path.split(tifFileName)[-1]
    print 'Adding tile index for {0}'.format(os.path.join(batchDir, baseFileName))
    cmdString = '"{0}" "{1}" "{2}"'.format(gdaltindexExe, tileIndexFileName,
        os.path.join(batchDir, baseFileName))
    subprocess.call(cmdString, shell=True, env=os.environ)


#
# for file in glob.glob(os.path.join(batchDir, 'o*_RGB_XCALIB.tif')):
#     baseFileName = os.path.split(file)[-1]  # 20 or 60m file
#     outFileName = baseFileName + '.ovr'
#     if os.path.exists(os.path.join(batchDir, outFileName)):
#         print '{0} exists'.format(os.path.join(batchDir, outFileName))
#     else:
#         print 'Generating overviews for {0}'.format(os.path.join(batchDir, outFileName))
#         cmdString = '"{0}" -ro -r average --config COMPRESS_OVERVIEW DEFLATE "{1}" 2 4 6 8 16 32 64'.format(gdaladdoExe, file)
#         subprocess.call(cmdString, shell=True, env=os.environ)
#
