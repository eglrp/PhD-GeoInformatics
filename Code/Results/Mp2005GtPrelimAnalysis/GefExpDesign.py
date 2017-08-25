import gdal
import ogr
import numpy as np
import osr
import pylab
from scipy import stats as stats
from matplotlib import patches
import matplotlib.pyplot as plt

# Following "Good Practice in Designing a Forest Inventory", calculate num samples etc for GEF experiments based on MP's 2005 data

csGtFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/Misc/BMR Carbon Stocks/abf_agc_191_plots.shp"

# read in cs gt
ds = gdal.OpenEx(csGtFile, gdal.OF_VECTOR)
if ds is None:
    print "Open failed./n"

lyr = ds.GetLayerByIndex(0)
lyr.ResetReading()
cs_gt_spatial_ref = lyr.GetSpatialRef()

# gcpList = []
cs_gt_dict = {}
for (i, feat) in enumerate(lyr):
    if i > 190:
        break
    print '.',
    feat_defn = lyr.GetLayerDefn()
    f = {}
    for i in range(feat_defn.GetFieldCount()):
        field_defn = feat_defn.GetFieldDefn(i)
        f[field_defn.GetName()] = feat.GetField(i)
    geom = feat.GetGeometryRef()
    if geom is not None and (geom.GetGeometryType() == ogr.wkbPoint or geom.GetGeometryType() == ogr.wkbPoint25D):
        print "%s %.6f, %.6f" % (f['PLOT'], geom.GetX(), geom.GetY())
        f['geom'] = geom
        f['X'] = geom.GetX()
        f['Y'] = geom.GetY()
    else:
        print "no point geometry/n"
    # gcpList.append(f)
    cs_gt_dict[f['PLOT']] = f
print ' '

ds = None

# for k in cs_gt_dict.keys():
#     # fixing some quirks in the keys
#     gpsK = k
#     gpsK = str(gpsK).replace('KADSS', 'KASS')
#     gpsK = str(gpsK).replace('GHOL1_', 'RHOL1_')
#     if gpsK == 'GHOL3_5':
#         gpsK = 'GHOL3_'

#

classi = np.zeros((cs_gt_dict.__len__()))

for i, plot in enumerate(cs_gt_dict.values()):
    if 'DST' in plot['PLOT']:
        classi[i] = 1
    elif 'ST' in plot['PLOT']:
        classi[i] = 2
    elif 'OL' in plot['PLOT']:
        classi[i] = 0
    else:
        classi[i] = 2
#

class_labels = ['OL','DST','ST', 'All']
tagc = np.log10(np.array([plot['TAGC'] for plot in cs_gt_dict.values()]))
cvi = np.zeros(3, float)
si = np.zeros(3, float)

pylab.figure()
for i, class_label in enumerate(class_labels[:-1]):
    class_idx = classi == i
    tagci = tagc[class_idx]
    #h = np.histogram(tagci, np.linspace(tagc.min(), tagc.max(), 20))
    h,hb = np.histogram(tagci, bins='auto', density=True)
    # pylab.subplot(1, 3, 1+i)
    pylab.plot(hb[:-1] + np.diff(hb)/2, h)
    cvi[i] = tagci.std()/np.abs(tagci.mean())
    si[i] = tagci.std()

h, hb = np.histogram(tagc, bins='auto', density=True)
# pylab.subplot(1, 3, 1+i)
pylab.plot(hb[:-1] + np.diff(hb) / 2, h, 'k')
pylab.xlabel('log(TAGC)')
pylab.legend(class_labels)

cv = tagc.std() / np.abs(tagc.mean())

# -----------------------------------------------
# Plot area effect on cv
# cv2**2 = (cv1**2)*np.sqrt(a1/a2)

print "%s CVi: %s" % (str(class_labels[:-1]), str(cvi))
print "5.x5. CV: %f" % (100*cv)

for sidelen in range(10, 30, 5):
    cv2 = np.sqrt(((100*cv)**2)*np.sqrt((5.*5.)/(sidelen**2)))
    print "%dx%d CVall: %f" % (sidelen, sidelen, cv2)

sidelens = np.array([25., 5., 5.])
for inc in range(5, 30, 5):
    sidelens_inc = sidelens + inc
    cvi2 = np.sqrt(((100*cvi)**2)*np.sqrt((sidelens**2)/(sidelens_inc**2)))
    print "Plot size: %s CVi: %s" % (str(sidelens_inc), str(cvi2))

# note the cv decreases non-lin and slower than sidelen i.e. diminishing returns
# also note that CV is markedly worse for OL and DST which suggests their plot sizes should be increased!


# -----------------------------------------------
# Initial sampling size

# NOTE: TAGC must not be log and must be in same area units as rest of qty's ??
L = 1.   # num strata
tst = 2. # t student value for 95% CI
E = 10.  # allowable error %
sidelen = 10.  # plot size
cvmax = np.max([cv, cvi.max()])
CV = np.sqrt(((100*cvi)**2)*np.sqrt((5.*5.)/(sidelen**2))) # cv*100 #np.max([cv, cvi.max()])*100 # highest CV in %
Ni = np.array([10000, 10000, 10000])  # max possible sampling units per stratum i.e. if plot is 1m2 and startum is 100m2 then Ni=100
N = np.sum(Ni) #

#from doc
n = ((tst**2)*(CV**2))/(E**2 + (tst**2 + CV**2)/N)
ni = n*Ni/N
print "Initial sampling size: %d" % (n)
print "Initial stratum sampling sizes: %s" % (str(ni))

# recalculation and realloation of sampling size (from doc)
# NOTE: this version looks fishy to me - how do they decide how to incorporate cost?  Prefer the above version
# in this case ni do not sum to n but can be less


Si = tagc.std()  #  si  # std dev of stratum i
Wi = Ni/N
Ci = np.array([1])  # cost per plot in stratum i

n = ((tst/E)**2) * (Wi*Si*np.sqrt(Ci)).sum() * (Wi*Si/np.sqrt(Ci)).sum()
ni = n * (Wi*Si/np.sqrt(Ci)) / ((Wi*Si/np.sqrt(Ci)).sum())
