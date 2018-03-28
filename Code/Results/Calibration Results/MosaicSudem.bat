SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox

CALL %OSGEO4W_ROOT%bin\o4w_env.bat
echo Paths....
setlocal EnableDelayedExpansion

W:
cd "W:\PhD GeoInformatics\Data\CGA\DEMs\SUDEM"

gdalwarp -co "COMPRESS=DEFLATE" x3321b_1_L2.tif x3322a_1_L2.tif x3321d_1_L2.tif x3322c_1_L2.tif SudemL2_LittleKaroo_Mosaic.tif
gdaladdo -r average -ro --config COMPRESS_OVERVIEW DEFLATE SudemL2_LittleKaroo_Mosaic.tif 2 4 8 16 32 64
pause