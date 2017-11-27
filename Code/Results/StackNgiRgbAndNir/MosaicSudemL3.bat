echo off
SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

CALL %OSGEO4W_ROOT%bin\o4w_env.bat

REM SET PATH=%PATH%;%OTB%\applications
REM SET PYTHONPATH=%QGIS%\python;%OTB%\python;%PYTHONPATH%
setlocal enabledelayedexpansion
cd D:\Data\Development\Projects\PhD GeoInformatics\Data\CGA\SUDEM L3 Unedited
D:

gdal_merge.py -o x3323da-b_2015_L3a.tif -n -150 -co "COMPRESS=LZW" -co "BIGTIFF=YES" -a_nodata -150 x3323da_2015_L3a.tif x3323db_2015_L3a.tif
gdaladdo -ro --config COMPRESS_OVERVIEW LZW –-config BIGTIFF_OVERVIEW YES -r average x3323da-b_2015_L3a.tif 2 4 8 16 32

gdal_merge.py -o x3324c_2015_L3a.tif -n -150 -co "COMPRESS=DEFLATE" -co "BIGTIFF=YES" -a_nodata -150 x3324ca_2015_L3a.tif x3324cb_2015_L3a.tif x3324cd_2015_L3a.tif
gdaladdo -ro --config COMPRESS_OVERVIEW LZW –-config BIGTIFF_OVERVIEW YES -r average x3324c_2015_L3a.tif 2 4 8 16 32
