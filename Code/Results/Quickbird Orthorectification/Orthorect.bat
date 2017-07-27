echo on

gdalwarp -rpc -to RPC_DEM="D:\Data\Development\Projects\MSc GeoInformatics\Data\CGA\SUDEM\x3324c_15_15_L2a.tif" "D:/Data/Development/Projects/MSc GeoInformatics/Data/Digital Globe/056549293010_01/056549293010_01_P001_PAN/03NOV18082012-P2AS_R1C1-056549293010_01_P001.TIF" 03NOV18082012-P2AS_R1C1-056549293010_01_P001_GdalRect.tif 

pause

cd "D:\Data\Development\Projects\MSc GeoInformatics\Data\DigitalGlobe\056549293010_01\056549293010_01_P001_PAN\"
otbcli_OrthoRectification -io.in 03NOV18082012-P2AS_R1C1-056549293010_01_P001.TIF -io.out 03NOV18082012-P2AS_R1C1-056549293010_01_P001_OtbRect.tif -elev.dem "D:\Data\Development\Projects\MSc GeoInformatics\Data\CGA\SUDEM\3324C" -interpolator "nn" -opt.ram 4096