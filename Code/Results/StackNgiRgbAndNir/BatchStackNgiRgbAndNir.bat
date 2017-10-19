@echo off
if -%3-==-- call :printhelp & exit /b
echo RGB dir: %1
echo NIR dir: %2
echo Out dir: %3
echo.

SET OSGEO4W_ROOT=C:\OSGeo4W64\

CALL %OSGEO4W_ROOT%bin\o4w_env.bat

setlocal EnableDelayedExpansion

set srcDir=%~dp0
echo Source dir: %srcDir%
cd %~dp0

echo Stacking.... 
for %%i in (%2\*n.tif) do (
REM CrossCalibration.exe %1 "%%i"
echo "%%i"
set rgbFile=%%~ni
REM set rgbFile=!rgbFile:~~0=_0!
set rgbFileL=!rgbFile:~0,18!
set rgbFileR=!rgbFile:~19!
set rgbFile=!rgbFileL!_!rgbFileR!
set rgbFile=!rgbFile:_n=_RGB.tif!
REM set tmp=%%~ni
REM set tmp=!tmp:~0,-2!
REM set outFile=!tmp!_RGBN.tif
set outFile=!rgbFile:_RGB.tif=_RGBN.tif!
REM echo rgbFileL: !rgbFileL!
REM echo rgbFileR: !rgbFileR!
echo rgbFile: !rgbFile!
echo outFile: !outFile!
REM echo tmp: !tmp!
.\gdal_merge.py -o %3\!outFile! -co "COMPRESS=JPEG" -co "JPEG_QUALITY=85" -co "NBITS=12" -a_nodata 0 -separate "%1\!rgbFile!" "%%i" 
REM gdaladdo --config COMPRESS_OVERVIEW JPEG -ro -r average RGBNtest2.tif 2 4 8 16 32
)

pause
REM "if %ERRORLEVEL%" doesn't actually work here as wtee overwrite the ERRORLEVEL with 0 but it is more important to have the log files 
REM if %ERRORLEVEL% equ 0 echo Binning.... & Bin.exe -c %~dp0\BinConfig.yaml %3000.bsg | wtee -a %3.log
REM if %ERRORLEVEL% equ 0 echo Exporting.... & VolExport.exe -f -z 99.9 -s %3@@ | wtee -a %3.log


:printhelp
echo.
echo Usage: BatchStackNgiRgbAndNir [rgb dir] [nir dir] [out dir]
goto :eof

pause 

