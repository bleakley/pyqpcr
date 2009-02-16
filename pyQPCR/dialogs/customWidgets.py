# -*- coding: utf-8 -*-
#
# pyQPCR, an application to analyse qPCR raw data
# Copyright (C) 2008 Thomas Gastine
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

from PyQt4.QtGui import *
from PyQt4.QtCore import QSize
from pyQPCR.wellGeneSample import *
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4 import NavigationToolbar2QT
from matplotlib.figure import Figure
import re

class GeneEchComboBox(QComboBox):

    def __init__(self, parent=None):
        QComboBox.__init__(self, parent)

    def addItem(self, obj):
        item = obj.name
        QComboBox.addItem(self, item)

    def addItems(self, listObj):
        self.listObj = listObj
        pix = QPixmap(32, 32)
        pix.fill(Qt.blue)
        for ind, obj in enumerate(self.listObj):
            self.addItem(obj)
            if obj.isRef == 2:
                self.setItemIcon(ind, QIcon(pix))

    def currentObj(self):
        index = QComboBox.currentIndex(self)
        return self.listObj[index]


class MatplotlibWidget(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        # We don't want the axes cleared every time plot() is called
        #self.axes.hold(True)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def sizeHint(self):
        w, h = self.get_width_height()
        return QSize(w, h)

    def minimumSizeHint(self):
        return QSize(10, 10)

    #def compute_initial_figure(self):
        #self.axes.plot([])
