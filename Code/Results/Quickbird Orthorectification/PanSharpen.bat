SET OSGEO4W_ROOT=C:\OSGeo4W64\
SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox
 
CALL %OSGEO4W_ROOT%bin\o4w_env.bat
 
SET PATH=%PATH%;%OTB%\applications
REM SET PYTHONPATH=%QGIS%\python;%OTB%\python;%PYTHONPATH%
echo off 
setlocal enabledelayedexpansion

set data=%1
set search=%2
set replace=%3
set "data=!data:%search%=%replace%!"
echo %data%

Set Pattern="P2AS"
Set Replace="M2AS"

REM For %%a in (*.jpg) Do (
REM     Ren "%%a" "!File:%Pattern%=%Replace%!"
REM )

REM Pause&Exit
D:
cd "D:\Data\Development\Projects\MSc GeoInformatics\Data\Digital Globe\056549293010_01"

REM With GDAL 
for /r %%i in (.\056549293010_01_P001_PAN\*P001.tif) do (
	echo %%~ni
    set MsFile=%%~ni
	:: Note that !x! must be used inside loop where variables are changed and not %x% because of "delayed local expansion" 
	REM set str=!str:%original%=%replacement%!
	set MsFile=!MsFile:P2AS=M2AS!
	echo !MsFile!
	REM gdal_pansharpen.py -w 0.7 -w 0.2 -w 0.1 panchro.tif multispectral.tif pansharpened_out.tif
	REM weights from  http://desktop.arcgis.com/en/arcmap/10.3/manage-data/raster-and-images/fundamentals-of-panchromatic-sharpening.htm 
	gdal_pansharpen.py -w 0.85 -w 0.7 -w 0.35 -w 1.0 .\056549293010_01_P001_PAN\%%~ni.tif .\056549293010_01_P001_MUL\!MsFile!.tif .\PanSharpen\%%~ni_GdalPanSharp.tif 
	gdaladdo -ro --config COMPRESS_OVERVIEW DEFLATE .\PanSharpen\%%~ni_GdalPanSharp.tif 4 8 16 32 64
)


pause
exit 

REM OTB version 
D:
cd "D:\Data\Development\Projects\MSc GeoInformatics\Data\Digital Globe\056549293010_01"

REM With OTB
for /r %%i in (.\056549293010_01_P001_PAN\*P001.tif) do (
	echo %%~ni
	REM resample MS image to PAN grid
	otbcli_BundleToPerfectSensor -inp .\056549293010_01_P001_PAN\%%~ni.tif -inxs .\056549293010_01_P001_MUL\%%~ni.tif -out .\PanSharpen\temp.tif 
	pause 
	otbcli_Pansharpening -inp .\056549293010_01_P001_PAN\%%~ni.tif -inxs .\PanSharpen\temp.tif -out .\PanSharpen\%%~ni"PanSharp.tif" uint16
	pause 
)

pause 
exit

pause 

gdalwarp -rpc -to RPC_DEM="D:\Data\Development\Projects\MSc GeoInformatics\Data\CGA\SUDEM\x3324c_15_15_L2a.tif" "D:/Data/Development/Projects/MSc GeoInformatics/Data/Digital Globe/056549293010_01/056549293010_01_P001_PAN/03NOV18082012-P2AS_R1C1-056549293010_01_P001.TIF" 03NOV18082012-P2AS_R1C1-056549293010_01_P001_GdalRect.tif 

pause

cd "D:\Data\Development\Projects\MSc GeoInformatics\Data\DigitalGlobe\056549293010_01\056549293010_01_P001_PAN\"
otbcli_OrthoRectification -io.in 03NOV18082012-P2AS_R1C1-056549293010_01_P001.TIF -io.out 03NOV18082012-P2AS_R1C1-056549293010_01_P001_OtbRect.tif -elev.dem "D:\Data\Development\Projects\MSc GeoInformatics\Data\CGA\SUDEM\3324C" -interpolator "nn" -opt.ram 4096
