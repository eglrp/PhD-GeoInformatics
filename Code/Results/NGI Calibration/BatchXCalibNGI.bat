@echo off
if -%2-==-- call :printhelp & exit /b
echo Ref file: %1
echo Input file spec: %2
echo.

echo off
SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

CALL %OSGEO4W_ROOT%bin\o4w_env.bat
echo Paths....
setlocal EnableDelayedExpansion
SET PATH=%path%;C:\Data\Development\Toolboxes\OpenCV-2.4.8\build\x64\vc12\bin
REM SET XCALIBEXE="C:\Data\Development\Projects\PhD GeoInformatics\Code\Cross Calibration\Sungis08\CrossCalibration.exe"
SET XCALIBEXE="C:\Data\Development\Projects\PhD GeoInformatics\Code\Cross Calibration\x64\Debug\CrossCalibration.exe"
SET MODISFILENAME="D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.A2015241.h20v12.005.2015260133130\MCD43A4.A2015241.h20v12.005.2015260133130.NgiBandOrder.Lo25Wgs84.tif"

REM gdalinfo -proj4 %QBFILENAME%

REM %XCALIBEXE% -w 1 1 -o %MODISFILENAME% %QBFILENAME%
echo Calibrating....
echo %XCALIBEXE%

echo Calibrating.... 
for %%i in (%2) do (
REM CrossCalibration.exe %1 "%%i"
echo %XCALIBEXE% -o -w 1 1 %MODISFILENAME% "%%i"
%XCALIBEXE% -o -w 1 1 %MODISFILENAME% "%%i"
)


REM "if %ERRORLEVEL%" doesn't actually work here as wtee overwrite the ERRORLEVEL with 0 but it is more important to have the log files 
REM if %ERRORLEVEL% equ 0 echo Binning.... & Bin.exe -c %~dp0\BinConfig.yaml %3000.bsg | wtee -a %3.log
REM if %ERRORLEVEL% equ 0 echo Exporting.... & VolExport.exe -f -z 99.9 -s %3@@ | wtee -a %3.log

goto :eof

:printhelp
echo.
echo Usage: BatchXCalib [ref file] [input file spec]
echo    [ref file]: Surf refl reference file in same projection and same band order as input files 
echo                (eg MODIS reference)
echo    [input file spec]: A wildcard pattern matching .tif files to be calibrated, 
echo                             can include path eg C:/dirName/*_RGBN.tif
goto :eof



pause
