echo off
SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

CALL %OSGEO4W_ROOT%bin\o4w_env.bat

REM SET PATH=%PATH%;%OTB%\applications
REM SET PYTHONPATH=%QGIS%\python;%OTB%\python;%PYTHONPATH%
setlocal enabledelayedexpansion
cd D:\Data\Development\Projects\PhD GeoInformatics\Data\CGA\SUDEM L3 Unedited
D:


gdalwarp -overwrite -t_srs "+proj=tmerc +lat_0=0 +lon_0=23 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" -r cubic -multi -srcnodata -150 -dstnodata -150 -of GTiff  -co "COMPRESS=LZW" -co "BIGTIFF=YES" "D:/Data/Development/Projects/PhD GeoInformatics/Data/CGA/SUDEM L3 Unedited/x3323db_2015_L3a.tif" "D:/Data/Development/Projects/PhD GeoInformatics/Data/CGA/SUDEM L3 Unedited/x3323db_2015_L3a_Lo23Wgs84.tif"

gdaladdo -ro --config COMPRESS_OVERVIEW LZW –-config BIGTIFF_OVERVIEW YES -r average x3323db_2015_L3a_Lo23Wgs84.tif 2 4 8 16 32

pause

REM make tile index (polygons of extents)

REM gdaltindex V:\Data\NGI\UnRectified\3323D_2015_1001\3323D_Extents.shp V:\Data\NGI\UnRectified\3323D_2015_1001\RGBN\*.tif