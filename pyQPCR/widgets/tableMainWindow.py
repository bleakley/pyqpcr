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
from numpy import nan

__author__ = "$Author: tgastine $"
__date__ = "$Date: 2010-10-14 21:02:08 +0200 (jeu. 14 oct. 2010) $"
__version__ = "$Rev: 324 $"


class MyQTableWidget(QTableWidget):
    """
    This class is used to reimplement the resizeEvent function
    of classical QTableWidget. It allows to have always a 12x8 table
    (even for 384-wells plate)
    """

    def resizeEvent(self, event):
        """
        This method allows to resize the plate in a 12x8 wells shape.
        """
        if self.columnCount() == 8 and self.rowCount() == 9:
            width = event.size().width()/8.
            height = event.size().height()/9.
        elif self.columnCount() == 10 and self.rowCount() == 10:
            width = event.size().width()/10.
            height = event.size().height()/10.
        else:
            width = event.size().width()/12.
            height = event.size().height()/8.
        for i in range(self.columnCount()):
            self.setColumnWidth(i, width)
        for j in range(self.rowCount()):
            self.setRowHeight(j, height)



class PlateWidget(MyQTableWidget):
    """
    This class allows to construct and populate a Q-PCR plate (A-H lines
    and 1-12 lines). The different table elements are filled depending on 
    their type. It is used in the main window of pyQPCR
    """

    def __init__(self, parent=None, plateType='96'):
        """
        Constructor of the Plate Widget. 

        :param parent: the QWidget parent
        :type parent: PyQt4.QtGui.QWidget
        :param plateType: the type of plates (96 wells, 384 wells, 72 or 100)
        :type plateType: PyQt4.QtCore.QString
        """
        QTableWidget.__init__(self, parent)

        if plateType == '96':
            self.tableLabels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
            self.setRowCount(8)
            self.setColumnCount(12)
        elif plateType == '384':
            self.tableLabels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
                                'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
            self.setRowCount(16)
            self.setColumnCount(24)
        elif plateType == '16':
            self.tableLabels = ['A']
            self.setRowCount(1)
            self.setColumnCount(16)
        elif plateType == '72':
            self.tableLabels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
                                'I']
            self.setRowCount(9)
            self.setColumnCount(8)
        elif plateType == '100':
            self.tableLabels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
                                'I', 'J']
            self.setRowCount(10)
            self.setColumnCount(10)

        self.setVerticalHeaderLabels(self.tableLabels)
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

        :param plate: the input plate of the table
        :type plate: pyQPCR.plate.Plaque
        """
        for well in plate.listePuits:
            if well.type in (QString('unknown'), QString('negative')):
                name = "%s\n%s" % (well.ech, well.gene.name)
            elif well.type == QString('standard'):
                try:
                    if well.amount >= 1e-2 and well.amount <= 1e3:
                        name = "%.2f\n%s" % (well.amount, well.gene)
                    else:
                        name = "%.1e\n%s" % (well.amount, well.gene)
                except TypeError:
                    name = "%s\n%s" % (well.amount, well.gene)
            tipname = "ct=%s\namount=%s" % (str(well.ct), str(well.amount))
            it = self.createItem(name, tip=tipname, status=well.name,
                                 back=well.type, icon=well.enabled)
            if well.warning == True and well.enabled == True:
                # if there is a warning and the well is enabled, then
                # we put the warning icon
                it.setIcon(QIcon(":/warning.png"))
            self.setItem(well.xpos, well.ypos, it)

    def createItem(self, text, tip=None, status=None, back=Qt.white,
                   fore=Qt.black, icon=None):
        """
        This method highly simplifies the creation of QTableWidgetItem

        :param text: the texte of the cell
        :type text: string
        :param tip: the tool tip of the cell
        :type tip: string
        :param status: the status tip of the cell
        :type status: string
        :param back: the background color of the cell
        :type back: PyQt4.QtGui.QColor
        :param fore: the foreground color of the cell
        :type fore: PyQt4.QtGui.QColor
        :param icon: the possible icon of the cell
        :type icon: string
        :return: the cell of the table
        :rtype: PyQt4.QtGui.QTableWidgetItem
        """
        item = QTableWidgetItem(text)
        item.setForeground(fore)
        if tip is not None:
            item.setToolTip(tip)
        if status is not None:
            item.setStatusTip(status)
        if icon is not None:
            if icon == False:
                item.setIcon(QIcon(":/disable.png"))
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
    Q-PCR experiment. It is used in the main window of pyQPCR
    """

    def __init__(self, parent=None, typeCalc='Relative quantification'):
        """
        Constructor of the Result Widget. 

        :param parent: the QWidget parent
        :type parent: PyQt4.QtGui.QWidget
        :param typeCalc: the type of calculation
        :type typeCalc: PyQt4.QtCore.QString
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

        :param plaque: the input plate of the table
        :type plaque: pyQPCR.plate.Plaque
        :param typeCalc: the type of calculation
        :type typeCalc: PyQt4.QtCore.QString
        """
        #self.setSortingEnabled(False)

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
                item.setIcon(QIcon(":/enable.png"))
                self.setVerticalHeaderItem(ind, item)
            else:
                item = QTableWidgetItem("")
                #item.setFont(QFont("Sans Serif", 16))
                item.setIcon(QIcon(":/disable.png"))
                self.setVerticalHeaderItem(ind, item)
            if well.warning == True and well.enabled == True:
                item.setIcon(QIcon(":/warning.png"))
            itWell = QTableWidgetItem(well.name)
            itWell.setFont(QFont("Sans Serif", 16))
            itGene = QTableWidgetItem(well.gene.name)
            try:
                itCt = QTableWidgetItem("%.2f" % well.ct)
            except TypeError:
                itCt = QTableWidgetItem("%s" % well.ct)
            if well.ctmean is nan:
                itCtmean = QTableWidgetItem('')
            else:
                try:
                    itCtmean = QTableWidgetItem("%.2f" % well.ctmean)
                except TypeError:
                    itCtmean = QTableWidgetItem(well.ctmean)
            if well.ctdev is nan:
                itCtdev = QTableWidgetItem('')
            else:
                try:
                    itCtdev = QTableWidgetItem("%.2f" % well.ctdev)
                except TypeError:
                    itCtdev = QTableWidgetItem(well.ctdev)
            try:
                if well.amount >= 1e-2 and well.amount <= 1e3:
                    itAmount = QTableWidgetItem("%.2f" % well.amount)
                else:
                    itAmount = QTableWidgetItem("%.2e" % well.amount)
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

        #self.setSortingEnabled(True)

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
