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

REM This 1st creates a VRT and then translates that into a tif.  It is faster than gdalwarp (and you have the VRT)

REM add tiffs to mosaic 
gdalbuildvrt %2.vrt %1 -a_srs "+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" -tr %tpix% %tpix% -tap -r cubicspline -srcnodata 0 -vrtnodata 0 
gdal_translate -of GTiff -co "COMPRESS=DEFLATE" -co "BIGTIFF=YES" -co "TILED=YES" -co "NBITS=12" %2.vrt %2
REM gdalwarp --config GDAL_CACHEMAX 3000 -wm 3000 -t_srs "+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" -tr %tpix% %tpix% -r cubicspline -tap -srcnodata 0 -dstnodata 0 %1 %2
gdaladdo -ro -r average --config COMPRESS_OVERVIEW DEFLATE %2.vrt 2 4 8 16 32 64
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

