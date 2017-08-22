echo off 
SET OSGEO4W_ROOT=C:\OSGeo4W64\
SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox
 
CALL %OSGEO4W_ROOT%bin\o4w_env.bat
 
REM SET PATH=%PATH%;%OTB%\applications
REM SET PYTHONPATH=%QGIS%\python;%OTB%\python;%PYTHONPATH%
setlocal enabledelayedexpansion

REM -------------------------------------------------------------
REM transform from ellipsoidal to geoid heights using sageoid2010
REM first run SaGeoidToGtx.py 
REM -------------------------------------------------------------
C:
cd "C:\Data\Development\Projects\PhD GeoInformatics\Docs\Misc\Baviaanskloof"
ogr2ogr -t_srs "+proj=longlat +datum=WGS84 +geoidgrids=sageoid2010_25.gtx +vunits=m +no_defs" BaviiaansPeCorrectedGcpMay2017Combined_SaGeoid2010.shp BaviiaansPeCorrectedGcpMay2017Combined.shp
pause
ogrinfo BaviiaansPeCorrectedGcpMay2017Combined_SaGeoid2010.shp BaviiaansPeCorrectedGcpMay2017Combined_SaGeoid2010
pause
ogr2ogr -t_srs "+proj=longlat +datum=WGS84 +geoidgrids=egm96_15.gtx +vunits=m +no_defs" BaviiaansPeCorrectedGcpMay2017Combined_Egm96.shp BaviiaansPeCorrectedGcpMay2017Combined.shp
pause
ogrinfo BaviiaansPeCorrectedGcpMay2017Combined_Egm96.shp BaviiaansPeCorrectedGcpMay2017Combined_Egm96
pause
