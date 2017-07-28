echo off 
SET OSGEO4W_ROOT=C:\OSGeo4W64\
SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox
 
CALL %OSGEO4W_ROOT%bin\o4w_env.bat
 
SET PATH=%PATH%;%OTB%\applications
REM SET PYTHONPATH=%QGIS%\python;%OTB%\python;%PYTHONPATH%
setlocal enabledelayedexpansion

REM Mosaic the 1st 2 tiles 

set mulR1C1="D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\056549293010_01\056549293010_01_P001_MUL\03NOV18082012-M2AS_R1C1-056549293010_01_P001.TIF"
set mulR1C2="D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\056549293010_01\056549293010_01_P001_MUL\03NOV18082012-M2AS_R1C2-056549293010_01_P001.TIF"
set outFile="D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\056549293010_01\056549293010_01_P001_MUL\Assembled\03NOV18082012-M2AS_R1C12-056549293010_01_P001.TIF"
set rpcFile="D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\056549293010_01\056549293010_01_P001_MUL\03NOV18082012-M2AS-056549293010_01_P001.RPB"
REM echo %mulR1C1%
gdal_merge.py -o %outFile% %mulR1C1% %mulR1C2%
copy %rpcFile% "D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\056549293010_01\056549293010_01_P001_MUL\Assembled\03NOV18082012-M2AS-R1C12-056549293010_01_P001.RPB"
gdaladdo -ro --config COMPRESS_OVERVIEW DEFLATE %outFile% 4 8 16 32 64


set panR1C1="D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\056549293010_01\056549293010_01_P001_PAN\03NOV18082012-P2AS_R1C1-056549293010_01_P001.TIF"
set panR1C2="D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\056549293010_01\056549293010_01_P001_PAN\03NOV18082012-P2AS_R1C2-056549293010_01_P001.TIF"
set outFile="D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\056549293010_01\056549293010_01_P001_PAN\Assembled\03NOV18082012-P2AS_R1C12-056549293010_01_P001.TIF"
set rpcFile="D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\056549293010_01\056549293010_01_P001_PAN\03NOV18082012-P2AS-056549293010_01_P001.RPB"
gdal_merge.py -o %outFile% %panR1C1% %panR1C2%
copy %rpcFile% "D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\056549293010_01\056549293010_01_P001_PAN\Assembled\03NOV18082012-P2AS-R1C12-056549293010_01_P001.RPB"
gdaladdo -ro --config COMPRESS_OVERVIEW DEFLATE %outFile% 4 8 16 32 64

pause
exit 
