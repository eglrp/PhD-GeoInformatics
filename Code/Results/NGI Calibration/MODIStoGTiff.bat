echo off
SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

CALL %OSGEO4W_ROOT%bin\o4w_env.bat

REM SET PATH=%PATH%;%OTB%\applications
REM SET PYTHONPATH=%QGIS%\python;%OTB%\python;%PYTHONPATH%
setlocal enabledelayedexpansion
cd D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.A2015241.h20v12.005.2015260133130
D:

echo "1"
gdal_translate -of Gtiff -sds MCD43A4.A2015241.h20v12.005.2015260133130.hdf ModisBand.tif
echo "2"
gdal_merge -of Gtiff -separate -a_nodata 0 -o "MCD43A4.A2015241.h20v12.005.2015260133130.tif" ModisBand_1.tif ModisBand_2.tif ModisBand_3.tif ModisBand_4.tif ModisBand_5.tif ModisBand_6.tif ModisBand_7.tif
echo "3"
gdal_merge -of Gtiff -separate -a_nodata 0 -o "MCD43A4.A2015241.h20v12.005.2015260133130.NgiBandOrder.tif" ModisBand_1.tif ModisBand_4.tif ModisBand_3.tif ModisBand_2.tif
pause

REM ModisBand_3.tif ModisBand_4.tif ModisBand_1.tif ModisBand_2.tif  QB order 