echo off
SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

CALL %OSGEO4W_ROOT%bin\o4w_env.bat

REM SET PATH=%PATH%;%OTB%\applications
REM SET PYTHONPATH=%QGIS%\python;%OTB%\python;%PYTHONPATH%
setlocal enabledelayedexpansion

REM Mosaic the 1st 2 tiles

set nirFile="E:\NIR\3323D_2015_1001\NIR\3323d_2015_1001_01~0002_n - scratch.tif"
set rgbFile="E:\Unrectified_Aerials\3323D_2015_1001\3323D_2015_1001_01_0002_RGB.tif"
REM gdal_merge.py -o %outFile% %mulR1C1% %mulR1C2%


echo %rgbFile%.vrt
REM gdal_translate -b mask -of vrt -a_nodata 0 %rgbFile% %rgbFile%.vrt
REM Note the  -a_nodata 0 doesn't seem to work when the mask is input, so do another pass
echo %rgbFile%2.vrt
REM gdal_translate -b 1 -of vrt -a_nodata 0 %rgbFile%.vrt %rgbFile%2.vrt
echo %rgbFile%.shp
REM gdal_polygonize.py -8 %rgbFile%2.vrt -b 1 -f "ESRI Shapefile" %rgbFile%.shp
pause
REM echo %nirFile%.vrt
REM gdal_translate -b mask -of vrt -a_nodata 0 %nirFile% %nirFile%.vrt
REM Note the  -a_nodata 0 doesn't seem to work when the mask is input, so do another pass
echo %nirFile%2.vrt
REM gdal_translate -b 1 -of vrt -a_nodata 0 %nirFile%.vrt %nirFile%2.vrt
echo %nirFile%.shp
REM gdal_polygonize.py -8 %nirFile%2.vrt -b 1 -f "ESRI Shapefile" %nirFile%.shp
pause

REM add NIR and RGB bands into one file to see if it works
REM gdal_merge.py -o E:\NIR\3323D_2015_1001\NIR\RGBNtest.tif -separate %rgbFile% %nirFile%
pause
REM gdalbuildvrt E:\NIR\3323D_2015_1001\NIR\RGBNtest.vrt -separate %rgbFile% %nirFile%
REM gdal_translate E:\NIR\3323D_2015_1001\NIR\RGBNtest.vrt E:\NIR\3323D_2015_1001\NIR\RGBNtest.vrt.tif
pause
gdal_merge.py -o E:\NIR\3323D_2015_1001\NIR\RGBNtest2.tif -co "COMPRESS=JPEG" -co "JPEG_QUALITY=85" -co "NBITS=12" -a_nodata 0 -separate %rgbFile% %nirFile%
gdaladdo --config COMPRESS_OVERVIEW JPEG -ro -r average RGBNtest2.tif 2 4 8 16 32
pause