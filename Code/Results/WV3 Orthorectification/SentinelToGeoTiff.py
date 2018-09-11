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


# if not os.path.exists(tileIndexFileName):
#     for file in glob.glob(os.path.join(sourceDir, tileIndexWildCard)):
#         subprocess.call('{0} "{1}" "{2}"'.format(gdaltindexExe, tileIndexFileName, file), shell=True, env=os.environ)
#         print 'Adding ' + file


##############################################################################################
# So the mapping should look something like
#
# WV3                           Sentinel
# 1 Coastal: 400 - 450 nm       b1 Coastal aerosol
# 2 Blue: 450 - 510 nm          b2 Blue
# 3 Green: 510 - 580 nm         b3 Green
# 4 Yellow: 585 - 625 nm        ? No corresponding band -  b3 is closest but does not overlap (checked vis in qgis)
# 5 Red: 630 - 690 nm           b4 Red
# 6 Red Edge: 705 - 745 nm      (b5+b6)/2  (Vegetation red edge + Vegetation red edge)/2 (checked vis in qgis)
# 7 Near-IR1: 770 - 895 nm      b8 NIR
# 8 Near-IR2: 860 - 1040 nm     (b8a + b9)/2 (Narrow NIR + Water Vapour)/2  or just b8a looks ok in qgis too
#
############################################################################
# upsample 20 and 60m bands to 10m so that these can be merged into one file

sentinelRootBandDir = "D:\Data\Development\Projects\PhD GeoInformatics\Data\Sentinel\S2A_MSIL2A_20171002T075741_N0205_R035_T34HGH_20171002T081741.SAFE - with SRTM90 all res\GRANULE\L2A_T34HGH_A011901_20171002T081741\IMG_DATA"
sentinel10mBandDir = os.path.join(sentinelRootBandDir, 'R10m')
resBandFileWildCard = 'L2A_T34HGH_20171002T075741_B*m.jp2'

# start with 20m directory as it is preferred to upsample these bands rather than their 60m counterparts
# then in 60m directory, check if the target file (from 20m or orig 10m) already exists before upsampling

# note that anaconda gdal does not support writing of jp2, so we make the upsampled files tiffs
for res in (20, 60):
    resDir = '{0}\R{1}m'.format(sentinelRootBandDir, res)
    for file in glob.glob(os.path.join(resDir, resBandFileWildCard)):
        baseFileName = os.path.split(file)[-1]                              # 20 or 60m file
        targetFileName = baseFileName.replace('_{0}m'.format(res), '_10m')  # 10m .jp2 filename
        outFileName = os.path.splitext(targetFileName)[0].replace('_10m', '_{0}to10m'.format(res)) + '.tif'              # 10m .tif filename
        if os.path.exists(os.path.join(sentinel10mBandDir, targetFileName)):
            print '{0} exists'.format(targetFileName)
        elif os.path.exists(os.path.join(sentinel10mBandDir, outFileName)):
            print '{0} exists'.format(outFileName)
        else:
            print 'Upsampling {0} to {1}'.format(os.path.split(file)[-1], os.path.split(outFileName)[-1])
            cmdString = '"{0}" -of GTiff -co "COMPRESS=DEFLATE" -multi -wo NUM_THREADS=ALL_CPUS -overwrite -srcnodata 0 -dstnodata 0 -wm 2048 -r cubicspline -tr 10 10 "{1}" "{2}"'.format(
                gdalwarpExe, file, os.path.join(sentinel10mBandDir, outFileName))
            subprocess.call(cmdString, shell=True, env=os.environ)
            # break

# after this has run, you can delete the 60to10m bands that have 20to10m equivalents.
# then use qgis to make (b5+b6)/2 and (b8a + b9)/2

#
# sprintf_s(gdalString, MAX_PATH,
#           "gdalwarp~-multi~-wo~NUM_THREADS=ALL_CPUS~-co~BIGTIFF=YES~-multi~-wo~NUM_THREADS=ALL_CPUS~-overwrite~-srcnodata~0~-dstnodata~0~-wm~2048~-r~cubicspline~-tr~%f~%f~%s~%s~",
#           fabs(srcGeoTransform[1]), fabs(srcGeoTransform[5]), paramDsFileName.c_str(), paramUsFileName.c_str());
#


sentinelBandDir = "D:\Data\Development\Projects\PhD GeoInformatics\Data\Sentinel\S2A_MSIL2A_20171002T075741_N0205_R035_T34HGH_20171002T081741.SAFE - No DEM all res\GRANULE\L2A_T34HGH_A011901_20171002T081741\IMG_DATA\R60m"
sentinelBandDir = "D:\Data\Development\Projects\PhD GeoInformatics\Data\Sentinel\S2A_MSIL2A_20171002T075741_N0205_R035_T34HGH_20171002T081741.SAFE - with SRTM90 all res\GRANULE\L2A_T34HGH_A011901_20171002T081741\IMG_DATA\R60m"
# L2A_T34HGH_20171002T075741_B01_60m.jp2
sentinelOutFileName = sentinelBandDir + "\L2A_T34HGH_20171002T075741_Wv3Bands_10m.tif"

# CRS EPSG:32734 - WGS 84 / UTM zone 34S - Projected
# USER:100027 - * Generated CRS (+proj=tmerc +lat_0=0 +lon_0=23 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs) - Projected

# gdal_merge -of Gtiff -separate -a_nodata 0 -o "MCD43A4.A2015241.h20v12.005.2015260133130.NgiBandOrder.tif" ModisBand_1.tif ModisBand_4.tif ModisBand_3.tif ModisBand_2.tif
senBandMap = ['B01_60to10m.tif', 'B02_10m.jp2', 'B03_10m.jp2', 'B03_10m.jp2', 'B04_10m.jp2',
              'MeanOfB05AndB06_20to10m.tif', 'B08_10m.jp2', 'MeanOfB8AandB09_20to10m.tif']   # how sentinel bands will be ordered to macth wv3

senBandMap = ['B01_60m.jp2', 'B02_10m.jp2', 'B03_60m.jp2', 'B03_60m.jp2', 'B04_60m.jp2',
              'MeanOfB05AndB06_60m.tif', 'B08_60m.jp2', 'MeanOfB8AandB09_60m.tif']   # how sentinel bands will be ordered to macth wv3

mergeCmdStr = '{0} -of Gtiff -separate -a_nodata 0 -o "{1}" '.format(gdalmergeExe, sentinelOutFileName)

for senBand in senBandMap:
    mergeCmdStr += '"{0}\{1}" '.format(sentinelBandDir, 'L2A_T34HGH_20171002T075741_{0}'.format(senBand))

if os.path.exists(sentinelOutFileName):
        os.remove(sentinelOutFileName)

subprocess.call(mergeCmdStr, shell=True, env=os.environ)

sentinelOutFileName2 = os.path.splitext(sentinelOutFileName)[0] + '_Cmp.tif'

# subprocess.call('{0} -co "COMPRESS=DEFLATE" "{1}" "{2}"'.format(gdaltranslateExe, sentinelOutFileName, sentinelOutFileName2),
#                 shell=True, env=os.environ)
# project to wv3 and compress
subprocess.call('{0} -t_srs "+proj=tmerc +lat_0=0 +lon_0=23 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" -tr 10.0 10.0 '
                '-r cubic -tap -co "COMPRESS=DEFLATE" "{1}" "{2}"'.format(gdalwarpExe, sentinelOutFileName, sentinelOutFileName2),
                shell=True, env=os.environ)


#gdalwarp -t_srs "+proj=tmerc +lat_0=0 +lon_0=23 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" -tr 10.0 10.0 -r cubic -tap -co "COMPRESS=DEFLATE" L2A_T34HGH_20171002T075741_Wv3Bands_10m_Cmp.tif L2A_T34HGH_20171002T075741_Wv3Bands_10m_CmpWv3Proj.tif


###################################################################################
# call xcalib for wv3 image
refFile = "D:\Data\Development\Projects\PhD GeoInformatics\Data\Sentinel\S2A_MSIL2A_20171002T075741_N0205_R035_T34HGH_20171002T081741.SAFE - with SRTM90 all res\GRANULE\L2A_T34HGH_A011901_20171002T081741\IMG_DATA\R10m\L2A_T34HGH_20171002T075741_Wv3Bands_10m_CmpWv3Proj.tif"
inputFile = "D:\Data\Development\Projects\PhD GeoInformatics\Data\Sentinel\S2A_MSIL2A_20171002T075741_N0205_R035_T34HGH_20171002T081741.SAFE - with SRTM90 all res\GRANULE\L2A_T34HGH_A011901_20171002T081741\IMG_DATA\R10m\L2A_T34HGH_20171002T075741_Wv3Bands_10m_CmpWv3Proj.tif"

xcalibExe = 'C:/Data/Development/Projects/PhD GeoInformatics/Code/Cross Calibration//x64/Release/CrossCalibration'

subprocess.call('"{0}" -o -w 3 3 -p 1 "{1}" "{2}"'.format(xcalibExe, refFile, inputFile), shell=True)
