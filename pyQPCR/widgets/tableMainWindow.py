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
import pyQPCR.qrc_resources

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"


class PlateWidget(QTableWidget):
    """
    This class allows to construct and populate a Q-PCR plate (A-H lines
    and 1-12 lines). The different table elements are filled depending on 
    their type.

    @ivar tableLabels: the list of wells names (i.e. A-H)
    @type tableLabels: list
    """

    def __init__(self, parent=None):
        """
        Constructor of the Plate Widget. 

        @param parent: the QWidget parent
        @type parent: PyQt4.QtGui.QWidget
        """
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

    def clear(self):
        """
        Overload the clear method by puting the tableLabels back after
        clearing.
        """
        QTableWidget.clear(self)
        self.setVerticalHeaderLabels(self.tableLabels)

    def populateTable(self, plate):
        """
        This method is used to fill the PCR plate thanks to the wells data.

        @param plate: the input plate of the table
        @type plate: pyQPCR.plate.Plaque
        """
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

        @param text: the texte of the cell
        @type text: string
        @param tip: the tool tip of the cell
        @type tip: string
        @param status: the status tip of the cell
        @type status: string
        @param back: the background color of the cell
        @type back: PyQt4.QtGui.QColor
        @param fore: the foreground color of the cell
        @type fore: PyQt4.QtGui.QColor
        @param icon: the possible icon of the cell
        @type icon: string
        @return: the cell of the table
        @rtype: PyQt4.QtGui.QTableWidgetItem
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
    """
    This class allows to construct and populate the result table of a
    Q-PCR experiment.

    @ivar resultLabels: the list of the names of the different columns used
                        in the result table
    @type resultLabels: list
    """

    def __init__(self, parent=None, typeCalc='Relative quantification'):
        """
        Constructor of the Result Widget. 

        @param parent: the QWidget parent
        @type parent: PyQt4.QtGui.QWidget
        @param typeCalc: the type of calculation
        @type typeCalc: PyQt4.QtCore.QString
        """
        QTableWidget.__init__(self, parent)
        if typeCalc == 'Relative quantification':
            self.resultLabels = ["Well", "Target", "Ct", "<Ct>", "E(Ct)", "Amount",
                    "Sample", "Eff", "Type", "NRQ"]
        elif typeCalc == 'Absolute quantification':
            self.resultLabels = ["Well", "Target", "Ct", "<Ct>", "E(Ct)", "Amount",
                    "Sample", "Eff", "Type", "Qabs"]
        self.setRowCount(96)
        self.setColumnCount(10)
        self.setHorizontalHeaderLabels(self.resultLabels)
        for i in range(len(self.resultLabels)):
            self.horizontalHeader().setResizeMode(i, QHeaderView.Stretch)
        self.setEditTriggers(QTableWidget.NoEditTriggers)
        self.setAlternatingRowColors(True)
        self.setSizePolicy(QSizePolicy(QSizePolicy.Maximum,
                                       QSizePolicy.Maximum))

        self.copyAction = QAction("Copy",  self)
        self.copyAction.setStatusTip('Copy data to the clipboard')
        self.copyAction.setToolTip('Copy data to the clipboard')
        self.copyAction.setIcon(QIcon(":/copy.png"))
        self.copyAction.setShortcut(QKeySequence.Copy)
        self.addAction(self.copyAction)

        self.setContextMenuPolicy(Qt.ActionsContextMenu)
        self.addAction(self.copyAction)

        self.connect(self.copyAction, SIGNAL("triggered()"), self.copyCells)

    def clear(self):
        """
        Overload the clear method by puting the headers of the table back after
        clearing.
        """
        QTableWidget.clear(self)
        self.setVerticalHeaderLabels(self.resultLabels)

    def populateResult(self, plaque, typeCalc):
        """
        This method is used to fill the PCR result thanks to the computed data?

        @param plaque: the input plate of the table
        @type plaque: pyQPCR.plate.Plaque
        """
        if typeCalc == 'Relative quantification':
            self.resultLabels = ["Well", "Target", "Ct", "<Ct>", "E(Ct)", "Amount",
                    "Sample", "Eff", "Type", "NRQ"]
        elif typeCalc == 'Absolute quantification':
            self.resultLabels = ["Well", "Target", "Ct", "<Ct>", "E(Ct)", "Amount",
                    "Sample", "Eff", "Type", "Qabs"]
        self.setHorizontalHeaderLabels(self.resultLabels)

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

    def copyCells(self):
        """
        A method to copy selected cells to the clipboard. The different
        columns are separated by tabulations by default.
        """
        selRange  = self.selectedRanges()[0]#just take the first range
        topRow = selRange.topRow()
        bottomRow = selRange.bottomRow()
        rightColumn = selRange.rightColumn()
        leftColumn = selRange.leftColumn()
        #item = self.tableWidget.item(topRow, leftColumn)
        clipStr = QString()
        for row in xrange(topRow, bottomRow+1):
            for col in xrange(leftColumn, rightColumn+1):
                cell = self.item(row, col)
                if cell:
                    clipStr.append(cell.text())
                else:
                    clipStr.append(QString(""))
                clipStr.append(QString("\t"))
            clipStr.chop(1)
            clipStr.append(QString("\r\n"))
        
        cb = QApplication.clipboard()
        cb.setText(clipStr)
