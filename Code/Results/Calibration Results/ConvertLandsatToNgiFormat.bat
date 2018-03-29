SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

REM CALL %OSGEO4W_ROOT%bin\o4w_env.bat

setlocal EnableDelayedExpansion
set PATH=%PATH%;"C:\Data\Development\Projects\PhD GeoInformatics\Code\Misc Tools"

echo on


REM SET SpotInFile="D:\Data\Development\Projects\PhD GeoInformatics\Data\SPOT\S131022114824832\Orthorectification\oATCORCorrected_METADATA_00812098_AutoGCP.tif"
SET LandsatInFile="D:\Data\Development\Projects\PhD GeoInformatics\Data\Landsat\LE71730832010034ASN00\LE07_L1TP_173083_20100203_20161217_01_T1_sr_MERGE.tif"
SET LandsatOutFile="D:\Data\Development\Projects\PhD GeoInformatics\Data\Landsat\LE71730832010034ASN00\LE07_L1TP_173083_20100203_20161217_01_T1_sr_MERGE_NgiFormat.tif"
SET LandsatOutFile2="D:\Data\Development\Projects\PhD GeoInformatics\Data\Landsat\LE71730832010034ASN00\LE07_L1TP_173083_20100203_20161217_01_T1_sr_MERGE_NgiFormat_SpotMask.tif"
SET NgiMosaicInFile="W:\PhD GeoInformatics\Data\NGI\Cross Calibration\XCalibMosaic30m.tif"
SET NgiMosaicOutFile="W:\PhD GeoInformatics\Data\NGI\Cross Calibration\XCalibMosaic30m_SpotMask.tif"

del tmp.tif 
goto thing
del %LandsatOutFile% 
REM gdal_translate ATCORCorrected_oS131022114824832b_139240261.tif ATCORCorrected_oS131022114824832b_139240261_Uint16.tif -ot UInt16 -scale 0 50 0 5000 -co "COMPRESS=LZW" -a_nodata 0
gdalwarp -srcnodata -9999 -dstnodata -9999 -r cubicspline -tr 30 30 -tap -t_srs "+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" %LandsatInFile% tmp.tif
REM RGBN band order 
gdal_translate tmp.tif %LandsatOutFile% -ot Int16 -co "COMPRESS=DEFLATE" -a_nodata -9999 -projwin 45690.0 -3652260.0 140730.0 -3765270.0 -b 3 -b 2 -b 1 -b 4
gdaladdo -r average -ro --config COMPRESS_OVERVIEW DEFLATE %LandsatOutFile% 2 5 8 16 32

del %LandsatOutFile2% 

gdalwarp -tr 30.0 -30.0 -tap -cutline "W:/PhD GeoInformatics/Data/Misc/spotExtent.shp" -crop_to_cutline -co "COMPRESS=DEFLATE" %LandsatOutFile% %LandsatOutFile2%
gdaladdo -r average -ro --config COMPRESS_OVERVIEW DEFLATE %LandsatOutFile2% 2 5 8 16 32
del tmp.tif 

REM make NGI mosaic of same extent etc 
del %NgiMosaicOutFile% 
gdalwarp -srcnodata 0 -dstnodata 0 -r cubicspline -tr 30 30 -tap  -cutline "W:/PhD GeoInformatics/Data/Misc/spotExtent.shp" -crop_to_cutline -co "COMPRESS=DEFLATE" %NgiMosaicInFile% %NgiMosaicOutFile%
gdaladdo -r average -ro --config COMPRESS_OVERVIEW DEFLATE %NgiMosaicOutFile% 2 5 8 16 32

:thing
del "W:/PhD GeoInformatics/Data/NGI/Cross Calibration/XCalibMosaic10m_SpotExtent.tif"
gdal_translate -projwin 45900.0 -3661740.0 112290.0 -3736410.0 -co "COMPRESS=DEFLATE" -b 4 -b 1 -b 2 -b 3 "W:/PhD GeoInformatics/Data/NGI/Cross Calibration/XCalibMosaic10m.tif" "W:/PhD GeoInformatics/Data/NGI/Cross Calibration/XCalibMosaic10m_SpotExtent.tif"
gdaladdo -r average -ro --config COMPRESS_OVERVIEW DEFLATE "W:/PhD GeoInformatics/Data/NGI/Cross Calibration/XCalibMosaic10m_SpotExtent.tif" 2 5 8 16 32

del tmp.tif 
pause 
REM -projwin 45690.0 -3652260.0 140730.0 -3765270.0

