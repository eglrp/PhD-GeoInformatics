import otbApplication
import os.path as path
import glob
import os

#This uses Orfeo - you must be using the osgeo4w python interpretor and pycharm must be started
# with"C:/OSGeo4W64/PyCharmOSGeo4W.bat"

# run "C:/OSGeo4W64/PyCharmOSGeo4W.bat" to use OSGeo4W stuff

panFiles = glob.glob("D:/Data/Development/Projects/MSc GeoInformatics/Data/Digital Globe/056549293010_01/"
                     "056549293010_01_P001_PAN/*P001.tif")

# ["D:/Data/Development/Projects/MSc GeoInformatics/Data/Digital Globe/056549293010_01/056549293010_01_P001_PAN/03NOV18082012-P2AS_R1C1-056549293010_01_P001.TIF"]
msFiles = glob.glob("D:/Data/Development/Projects/MSc GeoInformatics/Data/Digital Globe/056549293010_01/"
                     "056549293010_01_P001_MUL/*P001.tif")

outDir = "D:/Data/Development/Projects/MSc GeoInformatics/Data/Digital Globe/056549293010_01/PanSharpen/"
# The following line creates an instance of the Pansharpening application


for i in range(3, 4): #panFiles.__len__()):
    BundleToPerfectSensor = otbApplication.Registry.CreateApplication("BundleToPerfectSensor")
    Pansharpening = otbApplication.Registry.CreateApplication("Pansharpening")
    tempFile = outDir + "temp.tif"
    print "Resampling " + msFiles[i]
    os.remove(tempFile)
    BundleToPerfectSensor.SetParameterString("inp", panFiles[i])
    BundleToPerfectSensor.SetParameterString("inxs", msFiles[i])
    BundleToPerfectSensor.SetParameterString("out", tempFile)
    BundleToPerfectSensor.SetParameterString("ram", "4096")
    BundleToPerfectSensor.ExecuteAndWriteOutput()
    BundleToPerfectSensor = 0

    # The following lines set all the application parameters:
    print "Pan sharpening " + panFiles[i]
    Pansharpening.SetParameterString("inp", panFiles[i])
    Pansharpening.SetParameterString("inxs", tempFile)
    Pansharpening.SetParameterString("out", outDir + path.basename(panFiles[i]))
    Pansharpening.SetParameterOutputImagePixelType("out", 3)
    Pansharpening.SetParameterString("ram", "4096")
    # Pansharpening.SetParameterString("method", "bayes")
    # Does not work unless the input pan and ms files have same extent and pixel size...
    Pansharpening.ExecuteAndWriteOutput()
    Pansharpening = 0


