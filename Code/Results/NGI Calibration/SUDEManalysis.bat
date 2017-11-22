REM Generate hillshade etc SUDEM in GEF site.  To assist with stratification and sample site generation.
REM see http://www.gdal.org/gdaldem.html for other options 
REM echo off
SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

CALL %OSGEO4W_ROOT%bin\o4w_env.bat

setlocal EnableDelayedExpansion

SET DEMFILE="D:\Data\Development\Projects\PhD GeoInformatics\Data\CGA\SUDEM L3 Unedited\x3323db_2015_L3a.tif"
SET HSFILE="D:\Data\Development\Projects\PhD GeoInformatics\Data\CGA\SUDEM L3 Unedited\x3323db_2015_L3a_GDAL_HS.tif"
SET ASPECTFILE="D:\Data\Development\Projects\PhD GeoInformatics\Data\CGA\SUDEM L3 Unedited\x3323db_2015_L3a_GDAL_Aspect.tif"
SET ROUGHNESSFILE="D:\Data\Development\Projects\PhD GeoInformatics\Data\CGA\SUDEM L3 Unedited\x3323db_2015_L3a_GDAL_Roughness.tif"
REM Make hillshade with SUN at midday and altitude corresp to Sept 22 - see http://keisan.casio.com/exec/system/1224682277
echo Hillshade---------------------------------------------------------
gdaldem hillshade %DEMFILE% %HSFILE% -az 0 -alt 56 -alg Horn -co "COMPRESS=JPEG"
gdaladdo -ro -r average --config COMPRESS_OVERVIEW JPEG %HSFILE% 2 4 8 16 32
echo Aspect---------------------------------------------------------
gdaldem aspect %DEMFILE% %ASPECTFILE% -alg Horn -zero_for_flat -co "COMPRESS=DEFLATE"
gdaladdo -ro -r average --config COMPRESS_OVERVIEW DEFLATE %ASPECTFILE% 2 4 8 16 32
echo Roughness---------------------------------------------------------
gdaldem roughness %DEMFILE% %ROUGHNESSFILE% -co "COMPRESS=DEFLATE"
gdaladdo -ro -r average --config COMPRESS_OVERVIEW DEFLATE %ASPECTFILE% 2 4 8 16 32
pause
