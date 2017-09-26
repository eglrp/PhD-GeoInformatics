SET XCALIBEXE="C:\Data\Development\Projects\PhD GeoInformatics\Code\Cross Calibration\x64\Release\CrossCalibration.exe"
%XCALIBEXE%
SET MODISFILENAME="D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.A2003313.h20v12.005.2008026015724\MCD43A4.A2003313.h20v12.005.2008026015724.QbBandOrderAndCoords.tif"
SET QBFILENAME="D:\Data\Development\Projects\PhD GeoInformatics\Data\Digital Globe\056844553010_01\PCI Output\ATCOR1\ATCORCorrected_056844553010_01_P001_OrthoPanSharpen_05644032.tif"

gdalinfo -proj4 %QBFILENAME%

%XCALIBEXE% -w 1 1 -o %MODISFILENAME% %QBFILENAME%

pause
