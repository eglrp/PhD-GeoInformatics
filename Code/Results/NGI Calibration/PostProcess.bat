echo on
REM CMD /V:ON /C
setlocal EnableDelayedExpansion

for %%i in (*_RGBN_DS_2_5.tif) do (
set jj=%%i
REM  -co "TILED=YES" -co "BLOCKXSIZE=512" -co "BLOCKYSIZE=512" 
REM !jj:~1!
gdal_translate.exe -mo "BitsPerSample=12" -a_nodata 0 %%i ./PostProc/%%i
)

pause
exit
quit
