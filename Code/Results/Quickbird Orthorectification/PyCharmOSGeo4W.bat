@echo off
echo "OSGEO4W ENV SETUP..."
call "C:\OSGEO4W64\bin\o4w_env.bat"
@echo off
path %OSGEO4W_ROOT%\apps\qgis\bin;%PATH%
set QGIS_PREFIX_PATH=%OSGEO4W_ROOT:\=/%/apps/qgis
set GDAL_FILENAME_IS_UTF8=YES
rem Set VSI cache to be used as buffer, see #6448
set VSI_CACHE=TRUE
set VSI_CACHE_SIZE=1000000
set QT_PLUGIN_PATH=%OSGEO4W_ROOT%\apps\qgis\qtplugins;%OSGEO4W_ROOT%\apps\qt4\plugins;%QT_PLUGIN_PATH%
set PYTHONPATH=%OSGEO4W_ROOT%\apps\qgis\python;%PYTHONPATH%
REM set IPYTHON=1
REM set PYTHONHOME=C:\ProgramData\Anaconda3\envs\py27
echo "DONE"
SET PYCHARM="C:\Program Files\PyCharm Community Edition 2017.1.3\bin\pycharm64.exe"
 
start "PyCharm aware of QGIS & OSGeo4W" /B %PYCHARM% %*

REM SET OSGEO4W_ROOT=C:\OSGeo4W64\
REM SET QGISNAME=qgis
REM SET QGIS=%OSGEO4W_ROOT%apps\%QGISNAME%
REM SET QGIS_PREFIX_PATH=%QGIS%
REM SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox
 
REM CALL %OSGEO4W_ROOT%bin\o4w_env.bat
 
REM SET PATH=%PATH%;%QGIS%\bin;%OTB%\applications
REM SET PYTHONPATH=%QGIS%\python;%OTB%\python;%PYTHONPATH%

REM SET PYCHARM="C:\Program Files\PyCharm Community Edition 2017.1.3\bin\pycharm64.exe"
 
REM start "PyCharm aware of QGIS & OSGeo4W" /B %PYCHARM% %*
