gdal_translate -of Gtiff -sds MCD43A4.A2003313.h20v12.005.2008026015724.hdf ModisBand.tif
gdal_merge -of Gtiff -separate -a_nodata 0 -o "MCD43A4.A2003313.h20v12.005.2008026015724.tif" ModisBand_1.tif ModisBand_2.tif ModisBand_3.tif ModisBand_4.tif ModisBand_5.tif ModisBand_6.tif ModisBand_7.tif
gdal_merge -of Gtiff -separate -a_nodata 0 -o "MCD43A4.A2003313.h20v12.005.2008026015724.QbBandOrder.tif" ModisBand_3.tif ModisBand_4.tif ModisBand_1.tif ModisBand_2.tif

REM QB co-ord sys
REM AUTHORITY["EPSG","32735"]]
REM PROJ.4 string is:
REM '+proj=utm +zone=35 +south +datum=WGS84 +units=m +no_defs '

gdalwarp -r cubicspline -t_srs '+proj=utm +zone=35 +south +datum=WGS84 +units=m +no_defs' "MCD43A4.A2003313.h20v12.005.2008026015724.QbBandOrder.tif" "MCD43A4.A2003313.h20v12.005.2008026015724.QbBandOrderAndCoOrds.tif"

pause