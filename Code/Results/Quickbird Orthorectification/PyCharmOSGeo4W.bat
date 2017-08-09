SET OSGEO4W_ROOT=C:\OSGeo4W64\
SET QGISNAME=qgis
SET QGIS=%OSGEO4W_ROOT%apps\%QGISNAME%
SET QGIS_PREFIX_PATH=%QGIS%
SET OTB=%OSGEO4W_ROOT%apps\orfeotoolbox
SET PYCHARM="C:\Program Files\PyCharm Community Edition 2017.1.3\bin\pycharm64.exe"
 
CALL %OSGEO4W_ROOT%bin\o4w_env.bat
 
SET PATH=%PATH%;%QGIS%\bin;%OTB%\applications
SET PYTHONPATH=%QGIS%\python;%OTB%\python;%PYTHONPATH%
 
start "PyCharm aware of QGIS & OSGeo4W" /B %PYCHARM% %*

