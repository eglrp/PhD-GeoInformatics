REM convert PCI DEM output to single band DEM hopefully compatible with qgisthreejs
echo off
SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

CALL %OSGEO4W_ROOT%bin\o4w_env.bat

REM SET PATH=%PATH%;%OTB%\applications
REM SET PYTHONPATH=%QGIS%\python;%OTB%\python;%PYTHONPATH%
REM setlocal enabledelayedexpansion

cd "D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\DEM Working"
D:

SET PCIDEMFILE="D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\DEM Working\3324c_2015_1004_05_0183_to_0186_DEM_NoWallisOutBlend.pix"
SET OUTDEMFILE="D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\DEM Working\3324c_2015_1004_05_0183_to_0186_DEM_NoWallisOutBlend.tif"

gdal_translate -of Gtiff -co "COMPRESS=DEFLATE" -co "BIGTIFF=YES" -b 1 -a_nodata -150 %PCIDEMFILE% %OUTDEMFILE%
gdaladdo -ro --config COMPRESS_OVERVIEW DEFLATE -r average %OUTDEMFILE% 2 4 8 16 32

SET PCIDEMFILE="D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\DEM Working\3324c_2015_1004_05_0183_to_0186_DEM_WallisOutHighest.pix"
SET OUTDEMFILE="D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\DEM Working\3324c_2015_1004_05_0183_to_0186_DEM_WallisOutHighest.tif"

gdal_translate -of Gtiff -co "COMPRESS=DEFLATE" -co "BIGTIFF=YES" -b 1 -a_nodata -150 %PCIDEMFILE% %OUTDEMFILE%
gdaladdo -ro --config COMPRESS_OVERVIEW DEFLATE -r average %OUTDEMFILE% 2 4 8 16 32

pause
