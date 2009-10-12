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
from PyQt4.QtCore import *
from pyQPCR.wellGeneSample import *
import pyQPCR.qrc_resources
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4 import NavigationToolbar2QT
from matplotlib.figure import Figure
import re
import os
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

class NavToolBar(NavigationToolbar2QT):

    def __init__(self, canvas, parent, coordinates=True):
        NavigationToolbar2QT.__init__(self, canvas, parent, coordinates)

    def _init_toolbar(self):
        """
        modification of the toolbar definition in order to get
        the oxygen icons in the matplotib toolbar
        """
        a = self.addAction(QIcon(':/home.png'), 'Home', self.home)
        a.setToolTip('Reset original view')
        a = self.addAction(QIcon(':/undo.png'), 'Back', self.back)
        a.setToolTip('Back to previous view')
        a = self.addAction(QIcon(':/redo.png'), 'Forward', self.forward)
        a.setToolTip('Forward to next view')
        self.addSeparator()
        a = self.addAction(QIcon(':/move.png'), 'Pan', self.pan)
        a.setToolTip('Pan axes with left mouse, zoom with right')
        a = self.addAction(QIcon(':/zoom.png'), 'Zoom', self.zoom)
        a.setToolTip('Zoom to rectangle')
        self.addSeparator()
        a = self.addAction(QIcon(':/settings.png'), 'Subplots',
                self.configure_subplots)
        a.setToolTip('Configure subplots')
        a = self.addAction(QIcon(':/filesave.png'), 'Save',
                self.save_figure)
        a.setToolTip('Save the figure')

        self.buttons = {}

        if self.coordinates:
            self.locLabel = QLabel( "", self )
            self.locLabel.setAlignment(Qt.AlignRight|Qt.AlignTop )
            self.locLabel.setSizePolicy( QSizePolicy(QSizePolicy.Expanding,
                                  QSizePolicy.Ignored))
            labelAction = self.addWidget(self.locLabel)
            labelAction.setVisible(True)

        # reference holder for subplots_adjust window
        self.adj_window = None


