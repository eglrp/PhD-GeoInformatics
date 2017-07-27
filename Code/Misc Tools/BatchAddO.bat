@echo off
if -%1-==-- call :printhelp & exit /b
echo Input file spec: %1
echo.

setlocal EnableDelayedExpansion

echo Adding overviews.... 
for %%i in (%1) do (
REM gdaladdo -r AVERAGE -ro --config JPEG_QUALITY_OVERVIEW 75 --config COMPRESS_OVERVIEW JPEG --config GDAL_TIFF_OVR_BLOCKSIZE 512 "%%i" 2 4 8 16 32
echo "%%i":
gdaladdo -r AVERAGE -ro --config COMPRESS_OVERVIEW DEFLATE --config GDAL_TIFF_OVR_BLOCKSIZE 512 "%%i" 3 6 12 24 48
)


REM "if %ERRORLEVEL%" doesn't actually work here as wtee overwrite the ERRORLEVEL with 0 but it is more important to have the log files 
REM if %ERRORLEVEL% equ 0 echo Binning.... & Bin.exe -c %~dp0\BinConfig.yaml %3000.bsg | wtee -a %3.log
REM if %ERRORLEVEL% equ 0 echo Exporting.... & VolExport.exe -f -z 99.9 -s %3@@ | wtee -a %3.log

goto :eof

:printhelp
echo.
echo Adds external jpeg overviews
echo Usage: BatchAddO [input file spec]
echo    [input file spec]: A wildcard pattern matching .tif files to be calibrated, 
echo                             can include path eg C:/dirName/*_RGBN.tif
goto :eof


REM gdaladdo -r AVERAGE -ro --config COMPRESS_OVERVIEW DEFLATE --config GDAL_TIFF_OVR_BLOCKSIZE 512  16 32 64 128