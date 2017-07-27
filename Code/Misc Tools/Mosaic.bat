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

setlocal EnableDelayedExpansion

REM del mosaic.tif

gdalwarp -t_srs "+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" -tr %tpix% %tpix% -r cubicspline -tap -srcnodata 0 -dstnodata 0 %1 %2

goto :eof

:printhelp
echo.
echo Usage: Mosiac [input file spec] [output file]
echo    [input file spec]: A wildcard pattern matching .tif files to be mosaiced, 
echo                             can include path eg C:/dirName/*_RGBN.tif
echo    [output file]: Name of mosaic file
echo    [tgt pixel size]: Target pixel size in m (optional - default=500)
goto :eof

