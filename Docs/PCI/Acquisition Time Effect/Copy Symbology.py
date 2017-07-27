myMXD = arcpy.mapping.MapDocument("Current")
lstDataFrames=arcpy.mapping.ListDataFrames(myMXD)
lstLayers=arcpy.mapping.ListLayers(myMXD)

myOrthoLayer=arcpy.mapping.ListLayers(myMXD,'Later')

allFrames=arcpy.mapping.ListDataFrames(myMXD)

for dataFrame in allFrames:
     myMXD.activeView=dataFrame
     for groupLayer in myOrthoLayer:
		for imLayer in groupLayer:
			if imLayer.isRasterLayer:
				print "RASTA: " + str(imLayer.name)
				arcpy.ApplySymbologyFromLayer_management(imLayer, r"D:\Data\Development\Projects\MSc GeoInformatics\Docs\PCI\Acquisition Time Effect\o3321C_2010_318_01_0002_RGB.pix.lyr")