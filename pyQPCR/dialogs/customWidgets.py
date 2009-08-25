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
import pyQPCR.qrc_resources

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"


class HeaderModel(QStandardItemModel):

    def __init__(self, parent=None):
        QStandardItemModel.__init__(self, parent)

    def flags(self, index):
        flags =QStandardItemModel.flags(self, index)
        if not index.isValid():
            return flags
        value = index.data(Qt.UserRole).toString()
        if value == "header":
            flags &= Qt.ItemIsSelectable
            flags &= Qt.ItemIsEnabled
        return flags

class ComboDelegate(QItemDelegate):

    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)

    def paint(self, painter, option, index):
        text = index.model().data(index, Qt.DisplayRole).toString()
        if not (index.model().flags(index) and Qt.ItemIsEnabled):
            fontBold = QFont()
            fontBold.setBold(True)
            fontBold.setItalic(True)
            painter.setFont(fontBold)
            palette = QApplication.palette()
            painter.drawText(option.rect, Qt.AlignLeft|Qt.TextSingleLine,
                             text)
        else:
            QItemDelegate.paint(self, painter, option, index)

class GeneEchComboBox(QComboBox):

    def __init__(self, parent=None):
        QComboBox.__init__(self, parent)
        self.setModel(HeaderModel(self))
        self.setItemDelegate(ComboDelegate(self))

    def addItem(self, obj, *args):
        if hasattr(obj, "name"):
            item = obj.name
        else:
            item =  obj
        QComboBox.addItem(self, item, *args)

    def addItems(self, listObj):
        self.listObj = listObj
        for ind, obj in enumerate(self.listObj):
            self.addItem(obj)
            if hasattr(obj, "isRef") and obj.isRef == Qt.Checked:
                self.setItemIcon(ind+1, QIcon(":/reference.png"))

    def currentObj(self):
        index = QComboBox.currentIndex(self)
        nheader = 0
        for ind in range(QComboBox.count(self)):
            value = QComboBox.itemData(self, ind, Qt.UserRole).toString()
            if value == "header":
                nheader += 1
        return self.listObj[index-nheader]


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

