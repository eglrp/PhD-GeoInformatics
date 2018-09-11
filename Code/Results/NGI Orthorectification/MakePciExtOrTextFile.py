# convert geospace exterior orientation file to PCI format

import csv
import numpy as np

gpsInsFileNames = ["E:\Flight Data\\3318D_2016_1143\\3318D_2016_1143_GNSS_INS_DATA.txt",
               "E:\Flight Data\\3318B_2016_1142\\3318B_2016_1142_GNSS_INS_DATA.txt"]

imDirs = ["E:\Raw\\3318B_2016_1142", "E:\Raw\\3318D_2016_1143"]

imFileNameFormatStrings = ["3318B_2016_1142_{0:02d}_{1:04d}_RGB", "3318D_2016_1143_{0:02d}_{1:04d}_RGB"]

pciExtOrFileNames = ["E:\Raw\\3318B_2016_1142\\3318D_2016_1143_PCI_EXTOR.txt",
               "E:\Raw\\3318D_2016_1143\\3318B_2016_1142_PCI_EXTOR.txt"]

for i in range(0, gpsInsFileNames.__len__()):
    with open(gpsInsFileNames[i]) as gpsInsfile:
        gpsInsLines = gpsInsfile.readlines()
        headerEndLine = gpsInsLines.index('[End of Header]\n')
        gpsInsArray = np.genfromtxt(gpsInsLines[headerEndLine+1:])
        imNums = gpsInsArray[:,1]
        imFileNames = [imFileNameFormatStrings[i].format(int(row[0]), int(row[1])) for row in gpsInsArray]

    with open(pciExtOrFileNames[i], 'w') as pciExtOrFile:
        for imFileName, gpsInsRow in zip(imFileNames, gpsInsArray):
            print '{:s} {d[3]:f} {d[4]:f} {d[5]:f} {d[6]:f} {d[7]:f} {d[8]:f}'.format(imFileName, d=gpsInsRow)
            pciExtOrFile.write('{:s} {d[3]:f} {d[4]:f} {d[5]:f} {d[6]:f} {d[7]:f} {d[8]:f}\n'.format(imFileName, d=gpsInsRow))
