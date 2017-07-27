myMXD = arcpy.mapping.MapDocument("Current")
gtLayer=arcpy.mapping.ListLayers(myMXD,'Ground Truth')

for groupLayer in gtLayer:
	print "Group: " + str(groupLayer.name)	
	for layer in groupLayer:
		if layer.isFeatureLayer:
			print "Feature: " + str(layer.name)	
			try:
				arcpy.ApplySymbologyFromLayer_management(layer, r"D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\Groenfontein.lyr")
			except:	
				print "Error"