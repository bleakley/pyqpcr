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
        flags = QStandardItemModel.flags(self, index)
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
            try:
                item =  QString("%.2f" % obj)
            except TypeError:
                item = QString(obj)
        QComboBox.addItem(self, item, *args)

    def addItems(self, hashObj, editDialog = False):
        self.hashObj = hashObj
        if editDialog is False:
            for key in self.hashObj.keys():
                if key != '':
                    obj = self.hashObj[key]
                    self.addItem(obj)
                    if hasattr(obj, "isRef") and obj.isRef == Qt.Checked:
                        self.setItemIcon(self.hashObj.index(key), 
                                         QIcon(":/reference.png"))
        else:
            for key in self.hashObj.keys():
                obj = self.hashObj[key]
                self.addItem(obj)
                if hasattr(obj, "isRef") and obj.isRef == Qt.Checked:
                    self.setItemIcon(self.hashObj.index(key), 
                                     QIcon(":/reference.png"))

    def currentObj(self):
        objName = QComboBox.currentText(self)
        if objName != "header":
            return self.hashObj[objName]


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
        self.setIconSize(QSize(16, 16))

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

class PlateWidget(QTableWidget):

    def __init__(self, parent=None):
        QTableWidget.__init__(self, parent)
        self.tableLabels = ["A", "B", "C", "D", "E", "F", "G", "H"]
        self.setRowCount(8)
        self.setColumnCount(12)
        self.setVerticalHeaderLabels(self.tableLabels)
        for i in range(12):
            self.horizontalHeader().setResizeMode(i, QHeaderView.Stretch)
        for j in range(8):
            self.verticalHeader().setResizeMode(j, QHeaderView.Stretch)
        self.setEditTriggers(QTableWidget.NoEditTriggers)

#        self.connect(self, SIGNAL("cellDoubleClicked(int,int)"),
#                     self.editWell)

    def editWell(self):
        print "coucou"

    def clear(self):
        QTableWidget.clear(self)
        self.setVerticalHeaderLabels(self.tableLabels)

    def populateTable(self, plate):
        for well in plate.listePuits:
            if well.type in (QString('unknown'), QString('negative')):
                name = "%s\n%s" % (well.ech, well.gene.name)
            elif well.type == QString('standard'):
                try:
                    name = "%.2f\n%s" % (well.amount, well.gene)
                except TypeError:
                    name = "%s\n%s" % (well.amount, well.gene)
            tipname = "ct=%s\namount=%s" % (str(well.ct), str(well.amount))
            it = self.createItem(name, tip=tipname, status=well.name,
                                 back=well.type, icon=well.enabled)
            if well.warning == True and well.enabled == True:
                # if there is a warning and the well is enabled, then
                # we put the warning icon
                it.setIcon(QIcon(":/warning"))
            self.setItem(well.xpos, well.ypos, it)

    def createItem(self, text, tip=None, status=None, back=Qt.white,
                   fore=Qt.black, icon=None):
        """
        This method highly simplifies the creation of QTableWidgetItem
        """
        item = QTableWidgetItem(text)
        item.setForeground(fore)
        if tip is not None:
            item.setToolTip(tip)
        if status is not None:
            item.setStatusTip(status)
        if icon is not None:
            if icon == False:
                item.setIcon(QIcon(":/disable"))
        if back == QString('unknown'):
            item.setBackground(QColor(116, 167, 227))
        elif back == QString('standard'):
            item.setBackground(QColor(233, 0, 0))
        elif back == QString('negative'):
            item.setBackground(QColor(255, 250, 80))
        else:
            item.setBackground(Qt.white)
        item.setTextAlignment(Qt.AlignCenter|Qt.AlignVCenter)
        return item

class ResultWidget(QTableWidget):

    def __init__(self, parent=None):
        QTableWidget.__init__(self, parent)
        self.resultLabels=["Well", "Target", "Ct", "<Ct>", "E(Ct)", "Amount",
                "Sample", "Eff", "Type", "NRQ"]
        self.setRowCount(96)
        self.setColumnCount(10)
        self.setHorizontalHeaderLabels(self.resultLabels)
        for i in range(len(self.resultLabels)):
            self.horizontalHeader().setResizeMode(i, QHeaderView.Stretch)
        self.setEditTriggers(QTableWidget.NoEditTriggers)
        self.setAlternatingRowColors(True)
        self.setSizePolicy(QSizePolicy(QSizePolicy.Maximum,
                                       QSizePolicy.Maximum))

    def clear(self):
        QTableWidget.clear(self)
        self.setVerticalHeaderLabels(self.resultLabels)

    def populateResult(self, plaque):
        for ind, well in enumerate(plaque.listePuits):
            if well.enabled == True:
                item = QTableWidgetItem("")
                #item.setFont(QFont("Sans Serif", 16))
                item.setIcon(QIcon(":/enable"))
                self.setVerticalHeaderItem(ind, item)
            else:
                item = QTableWidgetItem("")
                #item.setFont(QFont("Sans Serif", 16))
                item.setIcon(QIcon(":/disable"))
                self.setVerticalHeaderItem(ind, item)
            if well.warning == True and well.enabled == True:
                item.setIcon(QIcon(":/warning"))
            itWell = QTableWidgetItem(well.name)
            itWell.setFont(QFont("Sans Serif", 16))
            itGene = QTableWidgetItem(well.gene.name)
            try:
                itCt = QTableWidgetItem("%.2f" % well.ct)
            except TypeError:
                itCt = QTableWidgetItem("%s" % well.ct)
            try:
                itCtmean = QTableWidgetItem("%.2f" % well.ctmean)
            except TypeError:
                itCtmean = QTableWidgetItem(well.ctmean)
            try:
                itCtdev = QTableWidgetItem("%.2f" % well.ctdev)
            except TypeError:
                itCtdev = QTableWidgetItem(well.ctdev)
            try:
                itAmount = QTableWidgetItem("%.2f" % well.amount)
            except TypeError:
                itAmount = QTableWidgetItem(well.amount)
            itEch = QTableWidgetItem(well.ech.name)
            itEff = QTableWidgetItem("%.2f%%%s%.2f" % (well.gene.eff,
                                     unichr(177), well.gene.pm))
            itType = QTableWidgetItem(well.type)
            try:
                itNRQ = QTableWidgetItem("%.2f%s%.2f" % (well.NRQ,
                                         unichr(177), well.NRQerror))
            except TypeError:
                itNRQ = QTableWidgetItem("%s%s" % (str(well.NRQ),
                                         str(well.NRQerror)))
            self.setItem(ind, 0, itWell)
            self.setItem(ind, 1, itGene)
            self.setItem(ind, 2, itCt)
            self.setItem(ind, 3, itCtmean)
            self.setItem(ind, 4, itCtdev)
            self.setItem(ind, 5, itAmount)
            self.setItem(ind, 6, itEch)
            self.setItem(ind, 7, itEff)
            self.setItem(ind, 8, itType)
            self.setItem(ind, 9, itNRQ)


