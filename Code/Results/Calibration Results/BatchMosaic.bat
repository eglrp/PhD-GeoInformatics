echo off
SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

CALL %OSGEO4W_ROOT%bin\o4w_env.bat

REM SET PATH=%PATH%;%OTB%\applications
REM SET PYTHONPATH=%QGIS%\python;%OTB%\python;%PYTHONPATH%
setlocal enabledelayedexpansion

set PATH=%PATH%;"C:\Data\Development\Projects\PhD GeoInformatics\Code\Misc Tools"

goto :master
REM call mosaic.bat "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3321B_2010_317\RGBN\*CMP_XCALIB.tif" "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3321B_2010_317\RGBN\3321B_Mosaic10m.tif" 10
call mosaic.bat "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3321B_2010_317\RGBN\*CMP_XCALIB.tif" "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3321B_2010_317\RGBN\3321B_Mosaic30m.tif" 30
call mosaic.bat "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3321B_2010_317\RGBN\*CMP_XCALIB.tif" "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3321B_2010_317\RGBN\3321B_Mosaic500m.tif" 500

call mosaic.bat "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3321D_2010_319\RGBN\*CMP_XCALIB.tif" "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3321D_2010_319\RGBN\3321D_Mosaic10m.tif" 10
call mosaic.bat "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3321D_2010_319\RGBN\*CMP_XCALIB.tif" "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3321D_2010_319\RGBN\3321D_Mosaic30m.tif" 30
call mosaic.bat "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3321D_2010_319\RGBN\*CMP_XCALIB.tif" "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3321D_2010_319\RGBN\3321D_Mosaic500m.tif" 500

call mosaic.bat "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3322A_2010_320\RGBN\*CMP_XCALIB.tif" "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3322A_2010_320\RGBN\3322A_Mosaic10m.tif" 10
call mosaic.bat "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3322A_2010_320\RGBN\*CMP_XCALIB.tif" "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3322A_2010_320\RGBN\3322A_Mosaic30m.tif" 30
call mosaic.bat "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3322A_2010_320\RGBN\*CMP_XCALIB.tif" "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3322A_2010_320\RGBN\3322A_Mosaic500m.tif" 500

call mosaic.bat "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3322C_2010_322\RGBN\*CMP_XCALIB.tif" "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3322C_2010_322\RGBN\3322C_Mosaic10m.tif" 10
call mosaic.bat "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3322C_2010_322\RGBN\*CMP_XCALIB.tif" "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3322C_2010_322\RGBN\3322C_Mosaic30m.tif" 30
call mosaic.bat "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3322C_2010_322\RGBN\*CMP_XCALIB.tif" "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3322C_2010_322\RGBN\3322C_Mosaic500m.tif" 500


copy "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3321B_2010_317\RGBN\*_Mosaic*.tif" "W:\PhD GeoInformatics\Data\NGI\Cross Calibration"
copy "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3321D_2010_319\RGBN\*_Mosaic*.tif" "W:\PhD GeoInformatics\Data\NGI\Cross Calibration"
copy "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3322A_2010_320\RGBN\*_Mosaic*.tif" "W:\PhD GeoInformatics\Data\NGI\Cross Calibration"
copy "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\3322C_2010_322\RGBN\*_Mosaic*.tif" "W:\PhD GeoInformatics\Data\NGI\Cross Calibration"

:master
W:
cd W:\PhD GeoInformatics\Data\NGI\Cross Calibration\

dir /b *_mosaic10m*.tif > mtiff_list.txt

gdal_merge -o XCalibMosaic10m.tif -v -n 0 -a_nodata 0 -co "TILED=YES" -co "BLOCKXSIZE=512" -co "BLOCKYSIZE=512" -co "COMPRESS=DEFLATE" --optfile mtiff_list.txt

dir /b *_mosaic30m*.tif > mtiff_list.txt

gdal_merge -o XCalibMosaic30m.tif -v -n 0 -a_nodata 0 -co "TILED=YES" -co "BLOCKXSIZE=512" -co "BLOCKYSIZE=512" -co "COMPRESS=DEFLATE" --optfile mtiff_list.txt

dir /b *_mosaic500m*.tif > mtiff_list.txt

gdal_merge -o XCalibMosaic500m.tif -v -n 0 -a_nodata 0 -co "TILED=YES" -co "BLOCKXSIZE=512" -co "BLOCKYSIZE=512" -co "COMPRESS=DEFLATE" --optfile mtiff_list.txt

pause