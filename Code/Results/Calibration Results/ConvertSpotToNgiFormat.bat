SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

REM CALL %OSGEO4W_ROOT%bin\o4w_env.bat

setlocal EnableDelayedExpansion
set PATH=%PATH%;"C:\Data\Development\Projects\PhD GeoInformatics\Code\Misc Tools"

echo on


REM SET SpotInFile="D:\Data\Development\Projects\PhD GeoInformatics\Data\SPOT\S131022114824832\Orthorectification\oATCORCorrected_METADATA_00812098_AutoGCP.tif"
SET SpotInFile="D:\Data\Development\Projects\PhD GeoInformatics\Data\SPOT\S131022114824832\ATCOR\ATCORCorrected_METADATA_14120009.pix"
SET SpotOutFile="D:\Data\Development\Projects\PhD GeoInformatics\Data\SPOT\S131022114824832\Orthorectification\oATCORCorrected_METADATA_00812098_AutoGCP_NgiFormat.tif"

del %SpotOutFile% tmp.tif 
REM gdal_translate ATCORCorrected_oS131022114824832b_139240261.tif ATCORCorrected_oS131022114824832b_139240261_Uint16.tif -ot UInt16 -scale 0 50 0 5000 -co "COMPRESS=LZW" -a_nodata 0
gdalwarp -srcnodata 0 -dstnodata 0 -r cubicspline -tr 10 10 -tap -t_srs "+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" %SpotInFile% tmp.tif
gdal_translate tmp.tif %SpotOutFile% -ot UInt16 -scale 0 50 0 5000 -co "COMPRESS=DEFLATE" -a_nodata 0 -projwin 45900.0 -3661740.0 112290.0 -3736410.0
gdaladdo -r average -ro --config COMPRESS_OVERVIEW DEFLATE %SpotOutFile% 2 5 8 16 32
del tmp.tif 
pause 
REM 