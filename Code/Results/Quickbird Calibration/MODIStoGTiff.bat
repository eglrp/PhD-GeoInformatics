D:
cd "D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.A2003313.h20v12.005.2008026015724"
gdal_translate -of Gtiff -sds MCD43A4.A2003313.h20v12.005.2008026015724.hdf ModisBand.tif
gdal_merge -of Gtiff -separate -a_nodata 0 -o "MCD43A4.A2003313.h20v12.005.2008026015724.tif" ModisBand_1.tif ModisBand_2.tif ModisBand_3.tif ModisBand_4.tif ModisBand_5.tif ModisBand_6.tif ModisBand_7.tif
gdal_merge -of Gtiff -separate -a_nodata 0 -o "MCD43A4.A2003313.h20v12.005.2008026015724.QbBandOrder.tif" ModisBand_3.tif ModisBand_4.tif ModisBand_1.tif ModisBand_2.tif