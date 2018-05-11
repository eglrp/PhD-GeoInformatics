import subprocess
import os
import glob
import shutil


########################################################################################################
spotFileNameNew = "D:\Data\Development\Projects\PhD GeoInformatics\Data\SPOT\S131022114824832\Orthorectification\oATCORCorrected_METADATA_00812098_AutoGCP_NgiFormat.tif"  #improved orthorect
ngiCalibFileName = "D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\StudyAreaXCalibMosaicSpotMask.tif"
# ngiCalibFileName = "W:\PhD GeoInformatics\Data\NGI\Cross Calibration\XCalibMosaic10m_SpotExtent.tif"   #seamline fix mosaic
outFileName = "D:\Data\Development\Projects\PhD GeoInformatics\Data\SPOT\S131022114824832\Orthorectification\AbsDiffWithNgi.tif"


# gdal_calc.py -A input.tif --outfile=result.tif --calc="A*(A>0)" --NoDataValue=0
bandFnList=[]
for b in [0,1,2]:
    (ofn, oext) = os.path.splitext(outFileName)
    ofnb = "%s%d%s" % (ofn, b, oext)
    bandFnList.append(ofnb)
    # subprocess.call('gdal_calc -A "{0}" --A_band={3} -B "{1}" --B_band={3} --overwrite --outfile="{2}" --calc="abs(A-B)"'.format(
    #     spotFileNameNew, ngiCalibFileName, ofnb, b+1), shell=True)

subprocess.call('gdal_merge -o "{3}" -separate "{0}" "{1}" "{2}"'.format(
    bandFnList[0], bandFnList[1], bandFnList[2], outFileName), shell=True)

subprocess.call('gdaladdo -ro -r cubicspline --config COMPRESS_OVERVIEW DEFLATE "{0}" 2 4 8 16 32 64'.format(outFileName),
                    shell=True)


