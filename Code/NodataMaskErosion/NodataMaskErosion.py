# %gui wx
# %matplotlib wx

import gdal
import cv2
import numpy as np
import matplotlib as mpl
import pylab

fileName = "F:/MSc GeoInformatics/Data/NGI/My Rectified/3322A_2010_320/SeamLineExperiment/3322a_320_02_0045_rgbn_CMP.tif"

imageDs = gdal.Open(fileName)
try:
    band = imageDs.GetRasterBand(1)
    band.GetNoDataValue()
    maskBand = band.GetMaskBand()
    mask = maskBand.ReadAsArray()
    pylab.figure()
    pylab.subplot(121)
    pylab.imshow(mask)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (500, 500))
    maskErode = cv2.erode(mask, kernel, iterations = 1)
    pylab.subplot(122)
    pylab.imshow(maskErode)

finally:
    imageDs = None


#
# class bob:
#     etc = 0
#
# @decorator(param=1)
# def f(x):
#     """ Syntax Highlighting Demo
#         @param x Parameter"""
#     s = ("Test", 2+3, {'a': 'b'}, x)   # Comment
#     print s[0].lower()
#
# class Foo:
#     def __init__(self):
#         byte_string = 'newline:\n also newline:\x0a'
#         text_string = u"Cyrillic Я is \u042f. Oops: \u042g"
#         self.makeSense(whatever=1)
#
#     def makeSense(self, whatever):
#         self.sense = whatever
#
# x = len('abc')
# print(f.__doc__)
