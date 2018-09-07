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

sentinelBandDir = "D:\Data\Development\Projects\PhD GeoInformatics\Data\Sentinel\S2A_MSIL2A_20171002T075741_N0205_R035_T34HGH_20171002T081741.SAFE - No DEM all res\GRANULE\L2A_T34HGH_A011901_20171002T081741\IMG_DATA\R60m"
sentinelBandDir = "D:\Data\Development\Projects\PhD GeoInformatics\Data\Sentinel\S2A_MSIL2A_20171002T075741_N0205_R035_T34HGH_20171002T081741.SAFE - with SRTM90 all res\GRANULE\L2A_T34HGH_A011901_20171002T081741\IMG_DATA\R60m"
# L2A_T34HGH_20171002T075741_B01_60m.jp2
sentinelOutFileName = sentinelBandDir + "\L2A_T34HGH_20171002T075741_Wv3Bands_60m.tif"

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
# 4 Yellow: 585 - 625 nm        ? No corresponding band -  b3 is closest but does not overlap
# 5 Red: 630 - 690 nm           b4 Red
# 6 Red Edge: 705 - 745 nm      (b5+b6)/2  (Vegetation red edge + Vegetation red edge)/2
# 7 Near-IR1: 770 - 895 nm      b8 NIR
# 8 Near-IR2: 860 - 1040 nm     (b8a + b9)/2 (Narrow NIR + Water Vapour)/2


# gdal_merge -of Gtiff -separate -a_nodata 0 -o "MCD43A4.A2015241.h20v12.005.2015260133130.NgiBandOrder.tif" ModisBand_1.tif ModisBand_4.tif ModisBand_3.tif ModisBand_2.tif
senBandMap = ['01', '02', '03', '03', '04', '05', '8A', '09']   # how sentinel bands will be ordered to macth wv3

mergeCmdStr = '{0} -of Gtiff -separate -a_nodata 0 -o "{1}" '.format(gdalmergeExe, sentinelOutFileName)

for senBand in senBandMap:
    mergeCmdStr += '"{0}\{1}" '.format(sentinelBandDir, '\L2A_T34HGH_20171002T075741_B{0}_60m.jp2'.format(senBand))

if os.path.exists(sentinelOutFileName):
        os.remove(sentinelOutFileName)

subprocess.call(mergeCmdStr, shell=True, env=os.environ)


##########################################################################
# upsample 20 and 60m bands to 10m

sentinelRootBandDir = "D:\Data\Development\Projects\PhD GeoInformatics\Data\Sentinel\S2A_MSIL2A_20171002T075741_N0205_R035_T34HGH_20171002T081741.SAFE - with SRTM90 all res\GRANULE\L2A_T34HGH_A011901_20171002T081741\IMG_DATA"
sentinel10mBandDir = os.path.join(sentinelRootBandDir, 'R10m')
resBandFileWildCard = 'L2A_T34HGH_20171002T075741_B*m.jp2'

for res in (20, 60):
    resDir = '{0}\R{1}m'.format(sentinelRootBandDir, res)
    for file in glob.glob(os.path.join(resDir, resBandFileWildCard)):
        outFile = os.path.join(sentinel10mBandDir, os.path.split(os.path.splitext(file)[0])[-1].replace(str(res),
            'UsTo{1}'.format(str(res),'10'))) + '.tif'
        if os.path.exists(outFile):
            print '{0} exists'.format(os.path.split(outFile)[-1])
        else:
            print 'Upsampling {0} to {1}'.format(os.path.split(file)[-1], os.path.split(outFile)[-1])
            cmdString = '"{0}" -of GTiff -co "COMPRESS=DEFLATE" -multi -wo NUM_THREADS=ALL_CPUS -overwrite -srcnodata 0 -dstnodata 0 -wm 2048 -r cubicspline -tr 10 10 "{1}" "{2}"'.format(
                gdalwarpExe, file, outFile)
            subprocess.call(cmdString, shell=True, env=os.environ)
            # break

#
# sprintf_s(gdalString, MAX_PATH,
#           "gdalwarp~-multi~-wo~NUM_THREADS=ALL_CPUS~-co~BIGTIFF=YES~-multi~-wo~NUM_THREADS=ALL_CPUS~-overwrite~-srcnodata~0~-dstnodata~0~-wm~2048~-r~cubicspline~-tr~%f~%f~%s~%s~",
#           fabs(srcGeoTransform[1]), fabs(srcGeoTransform[5]), paramDsFileName.c_str(), paramUsFileName.c_str());
#
