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
REM SET MODISFILENAME="D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.3324C_Mosaic.NgiBandOrder.Lo25Wgs84.tif"
SET MODISFILENAME="D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.A2016298.h19v12.006.2016307063215\MCD43A4.A2016302.h19v12.006.2016311094217.NgiBandOrder.TmLo19.tif"


REM gdalwarp -overwrite -t_srs "+proj=tmerc +lat_0=0 +lon_0=25 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" -r cubic -multi -srcnodata 0 -dstnodata 0 -of GTiff  -co "COMPRESS=LZW" %GEFFILE1% %GEFFILE1_WARP%

REM gdaladdo -ro --config COMPRESS_OVERVIEW LZW -r average %GEFFILE1_WARP% 2 4 8 16 32
cd "E:\Rectified\3318B_2016_1142\"
E:

for %%i in (o*.tif) do (
set jj=%%i
%XCALIBEXE% -w 1 1 -p 1 %MODISFILENAME% %%i
REM  -co "TILED=YES" -co "BLOCKXSIZE=512" -co "BLOCKYSIZE=512" 
REM !jj:~1!
REM gdal_translate.exe -mo "BitsPerSample=12" -a_nodata 0 %%i ./PostProc/%%i
echo %%i
REM gdaladdo -clean %%i
)

REM post proc xcalib files to make overviews
for %i in (o*_RGB_XCALIB.tif) do (
set jj=%i
echo %i
gdaladdo -ro -r average --config COMPRESS_OVERVIEW DEFLATE %i 2 4 8 16 32
)

pause
exit
REM 3323D...77
REM 3323D...79
REM 3323D...81
REM 3323D...82

%XCALIBEXE% -o -w 1 1 -p 1 %MODISFILENAME% "E:\Rectified\3318B_2016_1142\o3318B_2016_1142_01_0001_RGB.tif"
%XCALIBEXE% -o -w 1 1 -p 1 %MODISFILENAME% "E:\Rectified\3318B_2016_1142\o3318B_2016_1142_01_0002_RGB.tif"
%XCALIBEXE% -o -w 1 1 -p 1 %MODISFILENAME% "E:\Rectified\3318B_2016_1142\o3318B_2016_1142_01_0003_RGB.tif"
%XCALIBEXE% -o -w 1 1 -p 1 %MODISFILENAME% "E:\Rectified\3318B_2016_1142\o3318B_2016_1142_01_0004_RGB.tif"
%XCALIBEXE% -o -w 1 1 -p 1 %MODISFILENAME% "E:\Rectified\3318B_2016_1142\o3318B_2016_1142_01_0005_RGB.tif"

REM gdaladdo -ro -r average --config COMPRESS_OVERVIEW DEFLATE "E:\Rectified\3318B_2016_1142\o3318B_2016_1142_01_0001_RGB.tif" 2 4 8 16 32
gdaladdo -ro -r average --config COMPRESS_OVERVIEW DEFLATE "E:\Rectified\3318B_2016_1142\o3318B_2016_1142_01_0002_RGB_XCALIB.tif" 2 4 8 16 32
gdaladdo -ro -r average --config COMPRESS_OVERVIEW DEFLATE "E:\Rectified\3318B_2016_1142\o3318B_2016_1142_01_0003_RGB_XCALIB.tif" 2 4 8 16 32
gdaladdo -ro -r average --config COMPRESS_OVERVIEW DEFLATE "E:\Rectified\3318B_2016_1142\o3318B_2016_1142_01_0004_RGB_XCALIB.tif" 2 4 8 16 32
gdaladdo -ro -r average --config COMPRESS_OVERVIEW DEFLATE "E:\Rectified\3318B_2016_1142\o3318B_2016_1142_01_0005_RGB_XCALIB.tif" 2 4 8 16 32
pause 
exit 
