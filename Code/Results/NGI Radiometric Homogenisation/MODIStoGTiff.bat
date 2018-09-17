echo off
SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

CALL %OSGEO4W_ROOT%bin\o4w_env.bat

REM SET PATH=%PATH%;%OTB%\applications
REM SET PYTHONPATH=%QGIS%\python;%OTB%\python;%PYTHONPATH%
setlocal enabledelayedexpansion
cd "D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.A2016298.h19v12.006.2016307063215\"
D:

gdal_merge.py -separate -o MCD43A4.A2016302.h19v12.006.2016311094217.tif HDF4_EOS:EOS_GRID:"MCD43A4.A2016302.h19v12.006.2016311094217.hdf":MOD_Grid_BRDF:Nadir_Reflectance_Band1 HDF4_EOS:EOS_GRID:"MCD43A4.A2016302.h19v12.006.2016311094217.hdf":MOD_Grid_BRDF:Nadir_Reflectance_Band2 HDF4_EOS:EOS_GRID:"MCD43A4.A2016302.h19v12.006.2016311094217.hdf":MOD_Grid_BRDF:Nadir_Reflectance_Band3 HDF4_EOS:EOS_GRID:"MCD43A4.A2016302.h19v12.006.2016311094217.hdf":MOD_Grid_BRDF:Nadir_Reflectance_Band4 HDF4_EOS:EOS_GRID:"MCD43A4.A2016302.h19v12.006.2016311094217.hdf":MOD_Grid_BRDF:Nadir_Reflectance_Band5 HDF4_EOS:EOS_GRID:"MCD43A4.A2016302.h19v12.006.2016311094217.hdf":MOD_Grid_BRDF:Nadir_Reflectance_Band6 HDF4_EOS:EOS_GRID:"MCD43A4.A2016302.h19v12.006.2016311094217.hdf":MOD_Grid_BRDF:Nadir_Reflectance_Band7

gdal_merge.py -separate -o MCD43A4.A2016302.h19v12.006.2016311094217.NgiBandOrder.tif HDF4_EOS:EOS_GRID:"MCD43A4.A2016302.h19v12.006.2016311094217.hdf":MOD_Grid_BRDF:Nadir_Reflectance_Band1 HDF4_EOS:EOS_GRID:"MCD43A4.A2016302.h19v12.006.2016311094217.hdf":MOD_Grid_BRDF:Nadir_Reflectance_Band4 HDF4_EOS:EOS_GRID:"MCD43A4.A2016302.h19v12.006.2016311094217.hdf":MOD_Grid_BRDF:Nadir_Reflectance_Band3
 
gdalinfo -proj4 "E:\Rectified\3318B_2016_1142\o3318B_2016_1142_01_0001_RGB.tif" 

REM '+proj=tmerc +lat_0=0 +lon_0=19 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs '

gdalwarp -r cubicspline -t_srs "+proj=tmerc +lat_0=0 +lon_0=19 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" MCD43A4.A2016302.h19v12.006.2016311094217.NgiBandOrder.tif MCD43A4.A2016302.h19v12.006.2016311094217.NgiBandOrder.TmLo19.tif
pause
exit


echo "1"
gdal_translate -of Gtiff -sds MCD43A4.A2016298.h19v12.006.2016307063215.hdf ModisBandSet1.tif ModisBandSet2.tif
echo "2"
gdal_merge -of Gtiff -separate -a_nodata 0 -o "MCD43A4.A2016298.h19v12.006.2016307063215.hdf" ModisBand_08.tif ModisBand_09.tif ModisBand_10.tif ModisBand_11.tif ModisBand_12.tif ModisBand_13.tif ModisBand_14.tif
echo "3"
gdal_merge -of Gtiff -separate -a_nodata 0 -o "MCD43A4.A2016298.h19v12.006.2016307063215.NgiBandOrder.tif" ModisBand_08.tif ModisBand_11.tif ModisBand_10.tif ModisBand_9.tif
pause

REM ModisBand_3.tif ModisBand_4.tif ModisBand_1.tif ModisBand_2.tif  QB order 