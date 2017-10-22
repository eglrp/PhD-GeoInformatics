# ----------------------------------------------------------------------------------------------------------------------
# Read in Allometric model parameters

allometryFileName = "C:/Data/Development/Projects/PhD GeoInformatics/Data/GEF Field Trial/AllometricModels.xlsx"
woodyFileName = "C:/Data/Development/Projects/PhD GeoInformatics/Data/GEF Field Trial/GEF_Woody canopy_2017.10.16_Mdoda.xlsx"

from openpyxl import load_workbook

wb = load_workbook(allometryFileName)
ws = wb["Allometric Models"]

first_row = ws[0 + 1]
header = []
for c in first_row:
    header.append(c.value)

allometry_models = {}
for r in ws[2:ws.max_row]:  #how to find num rows?
    if r[0].value is None:
        break
    species = str(r[0].value[0]) + '. ' + str(r[1].value)
    model = {}
    model['vars'] = r[3].value
    model['ay'] = r[6].value
    model['by'] = r[7].value
    if str(r[14].value) == 'x':
        model['useWdRatio'] = False
    else:
        model['useWdRatio'] = True
    allometry_models[species] = model
wb = None
ws = None

# ----------------------------------------------------------------------------------------------------------------------
# Read in trial woody measurements

wb = load_workbook(woodyFileName)

plots = {}
for ws in wb:
    print ws.title, ' rows: ', ws.max_row
    first_row = ws[0 + 5]
    header = []
    for c in first_row:
        header.append(c.value)

    plot = []
    for r in ws[6:ws.max_row]:
        if r[1].value is None:
            break
        record = {}
        record['species'] = str(r[1].value)
        record['canopyWidth'] = r[2].value
        record['canopyLength'] = r[3].value
        record['height'] = r[4].value
        record['bsd'] = r[5].value
        plot.append(record)
    plots[ws.title] = plot
wb = None