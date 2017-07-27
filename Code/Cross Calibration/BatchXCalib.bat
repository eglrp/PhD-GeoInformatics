@echo off
if -%2-==-- call :printhelp & exit /b
echo Ref file: %1
echo Input file spec: %2
echo.

setlocal EnableDelayedExpansion

echo Calibrating.... 
for %%i in (%2) do (
CrossCalibration.exe %1 "%%i"
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
