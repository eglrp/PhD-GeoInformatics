REM Xcalib just the NGI images we need to covering current sampling points in GEF study site 
REM This is a bit of a hack as we are using the MODIS im for 3324c dates for both 3323d and 3324c images 

@echo off
SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

CALL %OSGEO4W_ROOT%bin\o4w_env.bat
echo Paths....
setlocal EnableDelayedExpansion

SET PATH=%path%;C:\Data\Development\Toolboxes\OpenCV-2.4.8\build\x64\vc12\bin
REM SET XCALIBEXE="C:\Data\Development\Projects\PhD GeoInformatics\Code\Cross Calibration\Sungis08\CrossCalibration.exe"
SET XCALIBEXE="C:\Data\Development\Projects\PhD GeoInformatics\Code\Cross Calibration\x64\Release\CrossCalibration.exe"
REM SET MODISFILENAME="D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.A2015241.h20v12.005.2015260133130\MCD43A4.A2015241.h20v12.005.2015260133130.NgiBandOrder.Lo25Wgs84.tif"
SET MODISFILENAME="D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.3324C_Mosaic.NgiBandOrder.Lo25Wgs84.tif"
REM gdalinfo -proj4 %QBFILENAME%
SET GEFFILE1="V:\Data\NGI\Rectified\3323D_2015_1001\RGBN\o3323d_2015_1001_02_0077_RGBN.tif"
SET GEFFILE1_WARP="V:\Data\NGI\Rectified\3323D_2015_1001\RGBN\o3323d_2015_1001_02_0077_Lo25Wgs84_RGBN.tif"
SET GEFFILE2="V:\Data\NGI\Rectified\3324C_2015_1004\RGBN\o3324c_2015_1004_02_0044_RGBN.tif"

REM gdalwarp -overwrite -t_srs "+proj=tmerc +lat_0=0 +lon_0=25 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" -r cubic -multi -srcnodata 0 -dstnodata 0 -of GTiff  -co "COMPRESS=LZW" %GEFFILE1% %GEFFILE1_WARP%

REM gdaladdo -ro --config COMPRESS_OVERVIEW LZW -r average %GEFFILE1_WARP% 2 4 8 16 32


%XCALIBEXE% -o -w 2 2 %MODISFILENAME% %GEFFILE1_WARP%
%XCALIBEXE% -o -w 2 2 %MODISFILENAME% %GEFFILE2%

pause
exit
REM 3323D...77
REM 3323D...79
REM 3323D...81
REM 3323D...82
