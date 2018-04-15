import subprocess
import os
import glob
import shutil
####################################################################################################
# NOTE
# CrossCalibration or any code linked with osgeo4w does cannot be called from anaconda
# - it needs to be called from osgeo4w python.  Perhaps the anaconda gdal/env interferes with
# the osgeo linking (?)  You can however call gdal utilities from anaconda ...
# Ah - I see now that Anaconda also comes with the gdal utilites and libs - so perhaps we should
# compile with these.  I'm pretty sure this is the issue now.
# We could dump osgeo completely and use gdal from anaconda and download qgis straight
# No, actually not a good a idea, anaconda does not update gdal as often as osgeo4w and
# when i link with the anaconda gdal version, i get a weird bug when applying params.  so lets
# just stick to osgeo4w

# run XCalib on files for specific parameters
# move output files to labelled folder
# mosaic output files to SPOT 10m
# (get SPOT in same proj and tap as xcalib files)
# read in mosaic and compare to spot as previously


########################################################################################################
# Run xcalib
inputDir = 'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\XCalib Experiments\Source2\\'
# inputFiles = ['3322a_320_09_0295_rgbn_CMP.tif', '3322a_320_09_0296_rgbn_CMP.tif', '3322a_320_09_0297_rgbn_CMP.tif', '3322a_320_09_0298_rgbn_CMP.tif',
#               '3322a_320_10_0359_rgbn_CMP.tif', '3322a_320_10_0360_rgbn_CMP.tif', '3322a_320_10_0361_rgbn_CMP.tif', '3322a_320_10_0362_rgbn_CMP.tif',
#               '3322a_320_11_0373_rgbn_CMP.tif', '3322a_320_11_0374_rgbn_CMP.tif', '3322a_320_11_0375_rgbn_CMP.tif', '3322a_320_11_0376_rgbn_CMP.tif']
# refFile = 'D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.A2010025.h19v12.005.2010043064233.Lo23.RGBN.tif'
refFile = "D:\Data\Development\Projects\PhD GeoInformatics\Data\MODIS\MCD43A4.A2010025.h19v12.005.2010043064233.Lo21.RGBN.tif"
xcalibExe = 'C:/Data/Development/Projects/PhD GeoInformatics/Code/Cross Calibration//x64/Release/CrossCalibration'
calibRootDir = 'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\XCalib Experiments\Calibrated2\\'

paramList = ['-w 1 1 -p 1', '-w 2 2 -p 1', '-w 3 3 -p 1', '-w 4 4 -p 1', '-w 5 5 -p 1', '-w 7 7 -p 1'] #, '-w 1 1 -p 4', '-w 3 1 -p 4', '-w 3 3 -p 4', '-w 5 5 -p 4']
paramList = ['-w 3 1 -p 1', '-w 1 1 -p 4', '-w 3 1 -p 4', '-w 3 3 -p 4', '-w 5 5 -p 4', '-w 7 7 -p 4'] #, '-w 1 1 -p 4', '-w 3 1 -p 4', '-w 3 3 -p 4', '-w 5 5 -p 4']
paramList = ['-w 3 3 -p 2', '-w 5 5 -p 2', '-w 7 7 -p 2'] #, '-w 1 1 -p 4', '-w 3 1 -p 4', '-w 3 3 -p 4', '-w 5 5 -p 4']
# paramList = ['-w 1 1 -p 4', '-w 2 2 -p 4', '-w 3 3 -p 4', '-w 4 4 -p 4', '-w 5 5 -p 4', '-w 7 7 -p 4'] #, '-w 1 1 -p 4', '-w 3 1 -p 4', '-w 3 3 -p 4', '-w 5 5 -p 4']
# paramList = ['-w 1 1 -p 4', '-w 3 1 -p 4', '-w 3 3 -p 4', '-w 5 5 -p 4']
# paramList = ['-w 3 3 -p 2', '-w 5 5 -p 2', '-w 7 7 -p 2']
# paramList = ['-w 7 7 -p 1', '-w 7 7 -p 4']
# paramList = ['-w 2 2 -p 1', '-w 2 2 -p 2', '-w 2 2 -p 4', '-w 4 4 -p 1', '-w 4 4 -p 2', '-w 4 4 -p 4'] #, '-w 1 1 -p 4', '-w 3 1 -p 4', '-w 3 3 -p 4', '-w 5 5 -p 4']
for param in  paramList:
    # inputFiles = os.listdir()
    for inputFile in glob.glob(os.path.join(inputDir, '*_CMP.tif')):
        print 'Processing {0:%s}' % (inputFile)
        subprocess.call('"{0}" -o {1} "{2}" "{3}"'.format(xcalibExe, param, refFile, inputFile), shell=True)
        # subprocess.call([xcalibExe, '-o', 'w 1 1', '-p 1', '"%s"'%refFile, '"%s"'%(inputDir+inputFile)], shell=True,
        #                 stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # move the xcalib files to a subdir
    subDir = param.replace(' ', '')
    subDir = subDir.replace('-', '')
    subDir = os.path.join(calibRootDir, subDir)
    if not os.path.exists(subDir):
        os.mkdir(subDir)
    for file in glob.glob(os.path.join(inputDir, '*_CMP_*.tif')):
        shutil.move(file, subDir)


########################################################################################################
# Add overviews and make mosaic
import os



# add overviews
# (root, dirs, files) = os.walk(calibRootDir)
for subDir in os.listdir(calibRootDir):
    if os.path.isdir(os.path.join(calibRootDir, subDir)):
        print subDir
        procList = []
        for file in os.listdir(os.path.join(calibRootDir, subDir)):
            print file
            if os.path.isfile(os.path.join(calibRootDir, subDir, file)):
                fn = os.path.join(calibRootDir, subDir, file)
                if fn.endswith('_CMP_XCALIB.tif') and not os.path.exists(fn+'.ovr'):
                    print 'gdaladdo ' + file
                    p = subprocess.call('gdaladdo -ro -r cubicspline --config COMPRESS_OVERVIEW DEFLATE "{0}" 2 4 8 16 32 64'.format(fn),
                                 shell=True)
                    #procList.append(p)
        # for p in procList:
        #     p.wait()
# make mosaics
for subDir in os.listdir(calibRootDir):
    if os.path.isdir(os.path.join(calibRootDir, subDir)):
        print subDir
        outFile = os.path.join(calibRootDir, subDir, subDir + '_Mosaic10m.tif')
        if not os.path.exists(outFile):
            print 'Mosaic ' + subDir
            subprocess.call('mosaic.bat "{0}" "{1}" 10'.format(os.path.join(calibRootDir, subDir, '*_CMP_XCALIB.tif'), outFile),
                            shell=True)

