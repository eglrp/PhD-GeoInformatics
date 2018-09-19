REM EG of how to make a VRT of tifs in a directory and make vrt coarse overviews (close-up overviews are provided by individual tif overviews)
REM The above is usable in QGIS although not fast
@echo off
SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

CALL %OSGEO4W_ROOT%bin\o4w_env.bat
echo Paths....
setlocal EnableDelayedExpansion

E:
cd "E:\Rectified\3318D_2016_1143"
dir /b *RGB.tif > dirlist.txt
gdalbuildvrt -input_file_list dirlist.txt 3318D_2016_1143_OrthoRect.vrt -srcnodata 0 -vrtnodata 0
gdaladdo -r average -ro --config COMPRESS_OVERVIEW DEFLATE 3318D_2016_1143_OrthoRect.vrt 32 64 128 256 512

cd "E:\Homogenised\3318D_2016_1143"
dir /b *XCALIB.tif > dirlist.txt
gdalbuildvrt -input_file_list dirlist.txt 3318D_2016_1143_Homog.vrt -srcnodata 0 -vrtnodata 0
gdaladdo -r average -ro --config COMPRESS_OVERVIEW DEFLATE 3318D_2016_1143_Homog.vrt 32 64 128 256 512
