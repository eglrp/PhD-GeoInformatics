echo off
SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

CALL %OSGEO4W_ROOT%bin\o4w_env.bat

REM SET PATH=%PATH%;%OTB%\applications
REM SET PYTHONPATH=%QGIS%\python;%OTB%\python;%PYTHONPATH%
setlocal enabledelayedexpansion

SET MODIS_FILE1="D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.A2015241.h19v12.005.2015260132233\MCD43A4.A2015241.h19v12.005.2015260132233.NgiBandOrder.tif"
SET MODIS_FILE2="D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.A2015241.h20v12.005.2015260133130\MCD43A4.A2015241.h20v12.005.2015260133130.NgiBandOrder.tif"

gdalwarp -overwrite -r cubic -srcnodata 0 -dstnodata 0 -t_srs "+proj=tmerc +lat_0=0 +lon_0=25 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" -tap -tr 500 500 %MODIS_FILE1% %MODIS_FILE2% "D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.3324C_Mosaic.NgiBandOrder.Lo25Wgs84.tif"

gdalwarp -overwrite -r cubic -srcnodata 0 -dstnodata 0 -t_srs "+proj=tmerc +lat_0=0 +lon_0=23 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" -tap -tr 500 500 %MODIS_FILE1% %MODIS_FILE2% "D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.Mosaic.NgiBandOrder.Lo23Wgs84.tif"


pause
REM gdal_merge.py -o x3324c_2015_L3a.tif -n -150 -co "COMPRESS=DEFLATE" -co "BIGTIFF=YES" -a_nodata -150 x3324ca_2015_L3a.tif x3324cb_2015_L3a.tif x3324cd_2015_L3a.tif
