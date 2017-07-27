@echo off
if -%2-==-- call :printhelp & exit /b
echo Input file: %1
type %1
echo Output file: %2
echo.

setlocal EnableDelayedExpansion
SET OSGEO4W_ROOT=C:\OSGeo4W64

python c:\osgeo4w64\bin\gdal_merge.py -o %2 -v -n 0 -a_nodata 0 -co "TILED=YES" -co "BLOCKXSIZE=512" -co "BLOCKYSIZE=512" -co "COMPRESS=LZW" --optfile %1

goto :eof

:printhelp
echo.
echo Usage: MasterMosiac [input file] [output file]
echo    [input file]: A text file containing the paths of tifs to mosaic/merge
echo    [output file]: Name of mosaic file
goto :eof

