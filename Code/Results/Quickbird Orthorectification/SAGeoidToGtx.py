#!/usr/bin/env python
#******************************************************************************
# 
#  Project:  Vertical Datum
#  Purpose:  Translate EGM2008 binary file into gtx format.
#  Author:   Frank Warmerdam, warmerdam@pobox.com
# 
#******************************************************************************
#  Copyright (c) 2010, Frank Warmerdam
# 
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
#******************************************************************************

import sys
import string
from osgeo import osr
from osgeo import gdal
import numpy as np
import pylab

# =============================================================================

# my adaptation to make a gtx file for sa geoid
# if __name__ == '__main__':
if True:
    # argv = gdal.GeneralCmdLineProcessor( sys.argv )
    # if argv is None:
    #     sys.exit( 0 )
    #
    # if len(argv) != 1:
    #     Usage()

    geoidFile = "D:/Data/Development/Projects/PhD GeoInformatics/Data/NGI/SA Geoid/SAGEOID2010.dat"
    dat = np.loadtxt(geoidFile, dtype=np.float32) #, delimiter = ',')
    ncols = 409
    nrows = 128017/ncols
    heights = np.reshape(dat[:, 2], (nrows, ncols))
    pylab.imshow(heights)
    lat_delta = np.median(np.diff(dat[:ncols-1, 1]))
    long_delta = np.median(np.diff(dat[::ncols, 0]))

#     vrt = """<VRTDataset rasterXSize="8640" rasterYSize="4321">
#   <VRTRasterBand dataType="Float32" band="1" subClass="VRTRawRasterBand">
#     <SourceFilename>Und_min2.5x2.5_egm2008_isw=82_WGS84_TideFree</SourceFilename>
#     <ImageOffset>4</ImageOffset>
#     <PixelOffset>4</PixelOffset>
#     <LineOffset>34568</LineOffset>
#     <ByteOrder>MSB</ByteOrder>
#   </VRTRasterBand>
# </VRTDataset>"""
#
#     ds_vrt = gdal.Open(vrt)
#     egm_m = ds_vrt.GetRasterBand(1).ReadAsArray()
#     ds_vrt = None

    gtx_file = 'D:/Data/Development/Projects/PhD GeoInformatics/Data/NGI/SA Geoid/sageoid2010_25.gtx'
    drv = gdal.GetDriverByName('GTX')
    ds_gtx = drv.Create(gtx_file,
                        ncols, nrows, 1, gdal.GDT_Float32 )
    ps_25 = 2.5 / 60.0
    # ds_gtx.SetGeoTransform( (-180 - ps_25,ps_25,0,90 + ps_25,0,-1 * ps_25) )

    # the below needs checking/fixing
    # ds_gtx.SetGeoTransform((np.min(dat[:,1]) - lat_delta, lat_delta, 0, np.max(dat[:, 0]) + long_delta, 0, -1 * long_delta))
    ds_gtx.SetGeoTransform((np.min(dat[:,1]) - ps_25/2, ps_25, 0, np.max(dat[:, 0]) + ps_25/2, 0, -1 * ps_25))

    # east = egm_m[:,:4320] * 1.0
    # west = egm_m[:,4320:] * 1.0

    ds_gtx.GetRasterBand(1).WriteArray(np.flipud(heights)) #( east, 4320, 0 )
    #ds_gtx.GetRasterBand(1).WriteArray( west, 0, 0 )
    
    ds_gtx = None
    
    import shutil
    shutil.copy2(gtx_file, 'C:/OSGeo4W64/share/proj/')