@echo off
if -%2-==-- call :printhelp & exit /b

set tpix=500
if -%3-==-- (
set tpix=500
) else (
set tpix=%3) 
)

echo Input file spec: %1
echo Output file: %2
echo Pixel size: %tpix%
echo.

REM echo off
SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

CALL %OSGEO4W_ROOT%bin\o4w_env.bat

setlocal EnableDelayedExpansion

REM del mosaic.tif
REM echo "%~d1%~p1"
REM %~d1
REM cd "%~d1%~p1"
REM dir /b %1 > tiff_list.txt
REM pause 

REM create mosaic file
REM echo gdal_merge.py -o tmp.tif -ps %tpix% %tpix% -v -n 0 -a_nodata 0 -tap -createonly --optfile tiff_list.txt 
REM gdal_merge.py -o tmp.tif -ps %tpix% %tpix% -v -n 0 -a_nodata 0 -tap -createonly --optfile tiff_list.txt 
REM gdal_merge.py -o %2 -co "COMPRESS=DEFLATE" -ps %tpix% %tpix% -v -n 0 -a_nodata 0 -tap -createonly %1
REM pause 
REM override projection
REM gdal_translate -a_srs "+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" tmp.tif %2
REM  pause 

REM add tiffs to mosaic 
gdalwarp --config GDAL_CACHEMAX 3000 -wm 3000 -t_srs "+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" -tr %tpix% %tpix% -r cubicspline -tap -srcnodata 0 -dstnodata 0 %1 %2

gdaladdo -ro -r average --config COMPRESS_OVERVIEW DEFLATE %2 2 4 8 16 32 64

goto :eof

:printhelp
echo.
echo Usage: Mosiac [input file spec] [output file]
echo    [input file spec]: A wildcard pattern matching .tif files to be mosaiced, 
echo                             can include path eg C:/dirName/*_RGBN.tif
echo    [output file]: Name of mosaic file
echo    [tgt pixel size]: Target pixel size in m (optional - default=500)
goto :eof

