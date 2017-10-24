# ----------------------------------------------------------------------------------------------------------------------
# Read in Allometric model parameters

allometryFileName = "C:/Data/Development/Projects/PhD GeoInformatics/Data/GEF Field Trial/AllometricModels.xlsx"
woodyFileName = "C:/Data/Development/Projects/PhD GeoInformatics/Data/GEF Field Trial/GEF_Woody canopy_2017.10.16_Mdoda.xlsx"

from openpyxl import load_workbook
import numpy as np
import pylab
from scipy.stats import gaussian_kde

def EvalRecordCs(allometricModels, record):
    # vars = [model['vars'] for model in allometricModels.values()]
    if not allometricModels.__contains__(record['species']):
        print record['species'], " not found"
        return 0.
    model = allometricModels[record['species']]
    x = 0.
    if model['vars'] == 'CA.H':
        x = np.pi*record['canopyLength']*record['canopyWidth']*record['height']/4.  # should this be elliptical or rectangular area
    elif model['vars'] == 'CA.SL':
        x = np.pi*record['canopyLength'] * record['canopyWidth'] * record['height']/4.   # assume cos measured these aloe speciosa as "hgt"
    elif model['vars'] == 'CD':
        x = np.mean([record['canopyLength'], record['canopyWidth']])
    elif model['vars'] == 'CD.H':
        x = np.mean([record['canopyLength'], record['canopyWidth']]) * record['height']
    elif model['vars'] == 'Hgt':
        x = record['height']
    else:
        print model['vars'], " unknown variable"
        return 0.
    yn = model['ay']*x**model['by']   # "naive"
    Yc = yn*model['duanz']
    return Yc

def EvalPlotCs(allometricModels, plot):
    plotYc = []
    for record in plot:
        Yc = EvalRecordCs(allometricModels, record)
        if Yc > 0.:
            plotYc.append(Yc)
        print record['species'], " - ", str(Yc)
    return plotYc

# ---------------------------------------------------------------------------------------------------------------------


wb = load_workbook(allometryFileName)
ws = wb["Allometric Models"]

first_row = ws[0 + 1]
header = []
for c in first_row:
    header.append(c.value)

allometricModels = {}
for r in ws[2:ws.max_row]:  #how to find num rows?
    if r[0].value is None:
        break
    species = str.strip(str(r[0].value[0]) + '. ' + str(r[1].value))
    model = {}
    model['vars'] = r[3].value
    model['ay'] = r[6].value
    model['by'] = r[7].value
    model['duanz'] = r[12].value
    if str(r[14].value) == 'x':
        model['useWdRatio'] = False
    else:
        model['useWdRatio'] = True
    allometricModels[species] = model
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
        species = str(r[1].value).strip()
        # hack for unknown plants
        if species == 'A. ferox' or species.__contains__('Aloe'):
            species = 'A. speciosa'
        # if species.__contains__('unk'):
        if species.__contains__('Asparagus') or species.__contains__('hairy'):
            species = 'A. capensis'
        if species.__contains__('Crassula') or species.__contains__('horns'):
            species = 'C. ovata'

        record['species'] = species
        record['canopyWidth'] = r[2].value
        record['canopyLength'] = r[3].value
        record['height'] = r[4].value
        record['bsd'] = r[5].value
        yc = EvalRecordCs(allometricModels, record)
        record['yc'] = yc
        plot.append(record)
    plots[ws.title] = plot
wb = None

# vars = [model['vars'] for model in allometricModels.values()]
# print np.unique(vars)
pylab.close('all')
i = 1
ycTtl = 0.
pylab.figure()
for plotKey, plot in plots.iteritems():
    yc = np.array([record['yc'] for record in plot])
    height = np.float64([record['height'] for record in plot])
    ycTtl += yc.sum()

    kde = gaussian_kde(height)  #, bw_method=bandwidth / height.std(ddof=1))
    heightGrid = np.linspace(0, 300, 100)
    heightKde = kde.evaluate(heightGrid)
    pylab.subplot(2, 3, i)
    pylab.plot(heightGrid, heightKde)
    pylab.grid('on')
    pylab.xlabel('Height (cm)')
    pylab.ylabel('KDE(Height)')
    pylab.title(plotKey)
    i += 1

i = 1
pylab.figure()
for plotKey, plot in plots.iteritems():
    yc = np.array([record['yc'] for record in plot])
    height = np.float64([record['height'] for record in plot])
    idx = np.argsort(height)

    ycCumSum = np.cumsum(-yc[idx])
    pylab.subplot(2, 3, i)
    pylab.plot(height[idx], ycCumSum)
    pylab.grid('on')
    pylab.xlabel('Height (cm)')
    pylab.ylabel('Cdf(CS)')
    pylab.xlim([0, 350])
    pylab.title(plotKey)
    i += 1

    idx = height > 50
    print str(plotKey), ": "
    print "Ttl: ", str(yc.sum())
    print "Hgt>50: ", str(yc[idx].sum())
    print "Hgt>50/Plot Ttl: ", str(yc[idx].sum()/yc.sum())
    print "Hgt>50/Ttl: ", str(yc[idx].sum()/ycTtl)
