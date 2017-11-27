echo off
SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

CALL %OSGEO4W_ROOT%bin\o4w_env.bat

setlocal EnableDelayedExpansion
REM echo on
REM CMD /V:ON /C
REM E:
REM cd "E:\Data\NGI\Rectified\3324C_2015_1004\RGBN"
V:
cd "V:\Data\NGI\Rectified\3323D_2015_1001\RGBN"

for %%i in (*_RGBN.tif) do (
set jj=%%i
REM  -co "TILED=YES" -co "BLOCKXSIZE=512" -co "BLOCKYSIZE=512" 
REM !jj:~1!
REM gdal_translate.exe -mo "BitsPerSample=12" -a_nodata 0 %%i ./PostProc/%%i
echo %%i
REM gdaladdo -clean %%i
gdaladdo -ro -r average --config COMPRESS_OVERVIEW PACKBITS %%i 2 4 8 16 32
)

pause
exit
quit
