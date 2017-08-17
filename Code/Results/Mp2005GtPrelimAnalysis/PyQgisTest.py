# import sip
# # switch on QString in Python3
# sip.setapi('QString', 1)

from qgis.gui import *
from qgis.core import *
from PyQt4.QtGui import QAction, QMainWindow, QMdiSubWindow
from PyQt4.QtCore import SIGNAL, Qt #, QString

def QString(s):
    return str(s)

class MyWnd(QMdiSubWindow):
  def __init__(self, layer):
    # QMainWindow.__init__(self)
    QMdiSubWindow.__init__(self)
    self.canvas = QgsMapCanvas()
    self.canvas.setCanvasColor(Qt.white)

    # self.canvas.setExtent(layer.extent())
    # self.canvas.setLayerSet([QgsMapCanvasLayer(layer)])

    # self.setCentralWidget(self.canvas)
    self.setWidget(self.canvas)

    actionZoomIn = QAction(QString("Zoom in"), self)
    actionZoomOut = QAction(QString("Zoom out"), self)
    actionPan = QAction(QString("Pan"), self)

    actionZoomIn.setCheckable(True)
    actionZoomOut.setCheckable(True)
    actionPan.setCheckable(True)

    self.connect(actionZoomIn, SIGNAL("triggered()"), self.zoomIn)
    self.connect(actionZoomOut, SIGNAL("triggered()"), self.zoomOut)
    self.connect(actionPan, SIGNAL("triggered()"), self.pan)

    # self.toolbar = self.addToolBar("Canvas actions")
    # self.toolbar.addAction(actionZoomIn)
    # self.toolbar.addAction(actionZoomOut)
    # self.toolbar.addAction(actionPan)
    #
    # # create the map tools
    self.toolPan = QgsMapToolPan(self.canvas)
    self.toolPan.setAction(actionPan)
    self.toolZoomIn = QgsMapToolZoom(self.canvas, False) # false = in
    self.toolZoomIn.setAction(actionZoomIn)
    self.toolZoomOut = QgsMapToolZoom(self.canvas, True) # true = out
    self.toolZoomOut.setAction(actionZoomOut)

    self.pan()

  def zoomIn(self):
    self.canvas.setMapTool(self.toolZoomIn)

  def zoomOut(self):
    self.canvas.setMapTool(self.toolZoomOut)

  def pan(self):
    self.canvas.setMapTool(self.toolPan)


import sys

if __name__ == "__main__":
    qgs = QgsApplication([], True)

    # load providers
    qgs.initQgis()
    # qgs.exec_()

    m = MyWnd(None)
    m.show()
    sys.exit(qgs.exec_())