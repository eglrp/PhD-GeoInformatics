
deadwoodFileName = "C:/Data/Development/Projects/PhD GeoInformatics/Data/GEF Field Trial/GEF_Deadwood_2017.10.16_Mdoda.xlsx"

from openpyxl import load_workbook
import numpy as np
import pylab
from scipy.stats import gaussian_kde
import collections



def EvalRecordCs(record):
    # dist in m !!
    # thumb suck params
    cfTree = 0.6
    Dj = 0.6
    L = 100   # length of transects
    densityReduction = [1., 0.8, 0.45]
    # formula from CDM doc: A/R Methodological Tool
    yc = (10.**2) * (44./12) * cfTree * Dj * (np.pi*np.pi/(8*L)) * ((record['diameter']/1000.)**2) * densityReduction[record['densityState']-1]
    return yc

def EvalPlotCs(allometricModels, plot):
    plotYc = []
    for record in plot:
        yc = EvalRecordCs(allometricModels, record)
        if yc > 0.:
            plotYc.append(yc)
    return plotYc


fontSize = 16

wb = load_workbook(deadwoodFileName)

sheetNames = ['INT-Lying Deadwood', 'INT-Standing Deadwood', 'TCH-Lying Deadwood', 'TCH-Standing Deadwood']
plots = collections.OrderedDict()
for sheetName in sheetNames:
    ws = wb[sheetName]
    for colStart in np.arange(0, ws.max_column, 2):
        if ws[7][colStart].value is None:
            break
        plotName = str(ws[7][colStart].value)
        plot = []
        print sheetName, ' - ', plotName
        for r in ws[9:ws.max_row]:
            if r[colStart].value is None:
                break
            record = {}
            record['diameter'] = r[colStart].value
            record['densityState'] = r[colStart + 1].value
            if record['densityState'] == 'D':
                record['densityState'] = 3L
            elif record['densityState'] == 'C':
                record['densityState'] = 2L
            elif record['densityState'] == 'B':
                record['densityState'] = 1L
            elif record['densityState'] == 'A':
                record['densityState'] = 1L
            elif type(record['densityState']) is not long:
                print 'WARNING: unknown density state: ', record['densityState']
            record['yc'] = EvalRecordCs(record)
            plot.append(record)
            print '.',
        if plots.__contains__(plotName):
            plots[plotName] = plots[plotName] + plot  # concat lying and standing deadwood
        else:
            plots[plotName] = plot
        print ''
wb = None

# vars = [model['vars'] for model in allometricModels.values()]
# print np.unique(vars)
pylab.close('all')
i = 1
ycTtl = 0.
pylab.figure()
for plotKey, plot in plots.iteritems():
    height = np.float64([record['diameter'] for record in plot])

    kde = gaussian_kde(height)  #, bw_method=bandwidth / height.std(ddof=1))
    heightGrid = np.linspace(0, 300, 100)
    heightKde = kde.evaluate(heightGrid)
    pylab.subplot(2, 3, i)
    pylab.plot(heightGrid, heightKde)
    axLim = pylab.axis()
    # h = pylab.plot([50, 50], [0, heightKde.max()], 'r')
    pylab.grid('on')
    pylab.xlabel('Plant Diam. (mm)', fontdict={'size':fontSize})
    pylab.ylabel('Prob.(diam.)', fontdict={'size':fontSize})
    pylab.title('Diam. distribution for plot %s' % (plotKey), fontdict={'size':fontSize})
    # if i >= plots.__len__():
    #     pylab.legend(h, ['50mm threshold'], loc='upper left', bbox_to_anchor=(1.2, 1), prop={'size':fontSize})
    pylab.axis([axLim[0], axLim[1], 0, heightKde.max()])
    i += 1

i = 1
pylab.figure()
for plotKey, plot in plots.iteritems():
    yc = np.array([record['yc'] for record in plot])
    height = np.float64([record['diameter'] for record in plot])
    idx = np.argsort(height)

    ycCumSum = np.cumsum(yc[idx])
    pylab.subplot(2, 3, i)
    pylab.plot(height[idx], ycCumSum)
    axLim = pylab.axis()
    # h = pylab.plot([50, 50], [axLim[2], axLim[3]], 'r')
    pylab.axis(axLim)
    pylab.grid('on')
    pylab.xlabel('Plant Diam. (mm)', fontdict={'size':fontSize})
    pylab.ylabel('Cum. Distr.(C. stock) (kg)', fontdict={'size':fontSize})
    # pylab.xlim([0, 350])
    pylab.title('Diam. / C. stock relation for plot %s' % (plotKey), fontdict={'size':fontSize})
    # if i >= plots.__len__():
    #     pylab.legend(h, ['50cm threshold'], loc='upper left', bbox_to_anchor=(1.2, 1), prop={'size':fontSize})
    i += 1

    idx = height > 50
    print str(plotKey), ": "
    print "Ttl: ", str(yc.sum())
    print "Hgt>50: ", str(yc[idx].sum())
    print "Hgt>50/Plot Ttl: ", str(yc[idx].sum()/yc.sum())
    # print "Hgt>50/Ttl: ", str(yc[idx].sum()/ycTtl)
