echo off
SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

CALL %OSGEO4W_ROOT%bin\o4w_env.bat

REM SET PATH=%PATH%;%OTB%\applications
REM SET PYTHONPATH=%QGIS%\python;%OTB%\python;%PYTHONPATH%
setlocal enabledelayedexpansion

Source file:
Dest file:
Source Geotransform:  (97406.2828861553, -0.5032597025761785, 0.001525988596180916, -3722620.2968298183, 0.0015259885961995678, 0.50325970257623)
Dest Geotransform:  (97411.63722928632, -0.5025025109868091, 0.00039498443239160963, -3722477.0546218366, 0.00041409720597584965, 0.5025863069134289)
0...10...20...30...40...50...60...70...80...90...100 - done.
"E:\NIR\3323D_2015_1001\NIR\3323d_2015_1001_03~0086_n.tif"


set nirFile="E:\NIR\3323D_2015_1001\NIR\3323d_2015_1001_03~0085_n.tif"
set rgbFile="E:\Unrectified_Aerials\3323D_2015_1001\3323d_2015_1001_03_0085_RGB.tif"

echo "copying"
"%~dp0\copy_geotransform.py" %rgbFile% %nirFile%
echo "copy done"
echo "merging1"
"C:\OSGeo4W64\bin\gdal_merge.py" -o "E:\Unrectified_Aerials\3323D_2015_1001\RGBN\85.tif" -separate %rgbFile% %nirFile%
echo "merge done"
pause
exit

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
REM gdal_merge.py -o E:\NIR\3323D_2015_1001\NIR\RGBNtest2.tif -co "COMPRESS=JPEG" -co "JPEG_QUALITY=85" -co "NBITS=12" -a_nodata 0 -separate %rgbFile% %nirFile%
REM gdaladdo --config COMPRESS_OVERVIEW JPEG -ro -r average RGBNtest2.tif 2 4 8 16 32
pause