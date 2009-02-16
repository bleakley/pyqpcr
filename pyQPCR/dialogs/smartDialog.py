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

import copy
import pyQPCR.qrc_resources
from PyQt4.QtGui import *
from PyQt4.QtCore import *

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"

def fromPosition(x, y):
    lettres = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    dict = {}
    for xpos, l in enumerate(lettres):
        dict[xpos] = l
    if dict.has_key(x):
        return dict[x]+str(y)

class SmartDialog(QDialog):

    def __init__(self, parent=None, plaque=None):
        self.parent = parent
        QDialog.__init__(self, parent)
# Widgets
        lab1 = QLabel("&Wells per replicate:")
        self.trip = QSpinBox()
        self.trip.setRange(1, 6)
        lab1.setBuddy(self.trip)
        lab2 = QLabel("&Direction:")
        self.cboxDir = QComboBox()
        self.cboxDir.addItems(["Horizontal", "Vertical"])
        self.cboxDir.setItemIcon(0, QIcon(":/arrow-right.png"))
        self.cboxDir.setItemIcon(1, QIcon(":/arrow-down.png"))
        lab2.setBuddy(self.cboxDir)
        lab3 = QLabel("&Initial amount:")
        self.editAmount = QLineEdit()
        lab3.setBuddy(self.editAmount)
        lab4 = QLabel("&Dilution factor:")
        self.editDil = QLineEdit()
        lab4.setBuddy(self.editDil)
#
        layout = QGridLayout()
        layout.addWidget(lab1, 0, 0)
        layout.addWidget(self.trip, 0, 1)
        layout.addWidget(lab2, 1, 0)
        layout.addWidget(self.cboxDir, 1, 1)
        layout.addWidget(lab3, 2, 0)
        layout.addWidget(self.editAmount, 2, 1)
        layout.addWidget(lab4, 3, 0)
        layout.addWidget(self.editDil, 3, 1)
#
        self.but = QPushButton("&Preview")
        self.but.setCheckable(True)
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|
                                     QDialogButtonBox.Cancel)
        buttonBox.addButton(self.but,QDialogButtonBox.ActionRole)
#
        
        self.table = QTableWidget()
        tableLabels=["A", "B", "C", "D", "E", "F", "G", "H"]
        self.table.setRowCount(8)
        self.table.setColumnCount(12)
        self.table.setVerticalHeaderLabels(tableLabels)
        self.table.resizeRowsToContents()
        self.table.resizeColumnsToContents()
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionMode(QTableWidget.NoSelection)
#
        self.table.hide()
#
        mainLayout = QVBoxLayout()
        mainLayout.setSizeConstraint(QLayout.SetFixedSize)
        mainLayout.addLayout(layout)
        mainLayout.addWidget(self.table)
        mainLayout.addWidget(buttonBox)
        self.setLayout(mainLayout)
#
        self.connect(buttonBox, SIGNAL("accepted()"), self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"), self, SLOT("reject()"))
        self.connect(self.but, SIGNAL("toggled(bool)"), self.preview)
        self.setWindowTitle("Smart")

    def preview(self):
        self.table.setVisible(self.but.isChecked())
        if self.parent is not None:
            #trip = []
            #k = 1
            setX = set() ; setY = set()
            for it in self.parent.table.selectedItems():
                nom = str(it.statusTip())
                well = getattr(self.plaque, nom)
                setX.add(it.row())
                setY.add(it.column())

            xmin, xmax = min(list(setX)), min(list(setX))
            ymin, ymax = min(list(setY)), min(list(setY))
            if self.cboxDir.currentText() == "Horizontal":
                for xind in range(xmin, xmax+1):
                    for yind in range(ymin, ymax+1, int(self.trip.value())):
                        name = fromPosition(xind, yind)
                        well = getattr(self.plaque, nom)
                        if str(well.type) == 'standard':
                            well.setAmount(float(self.editAmount.text())/coeff)
            elif self.cboxDir.currentText() == "Vertical":
                for yind in range(ymin, ymax+1):
                    for xind in range(xmin, xmax+1, int(self.trip.value())):
                        for indice in range(int(self.trip.value())):
                            name = fromPosition(xind+indice, yind)
                            well = getattr(self.plaque, nom)
                            if str(well.type) == 'standard':
                                well.setAmount(float(self.editAmount.text())/coeff)

                """
                if str(well.type) == 'standard':
                    trip.append(well)
                    if len(trip) == int(self.trip.value()):
                        coeff = k*float(self.editDil.text())
                        for well in trip:
                            well.setAmount(float(self.editAmount.text())/coeff)
                        trip = []
                        k += 1
                    well.setAmount
                """

        print "Not implemented yet"


if __name__=="__main__":
    import sys
    from plaque import *
    app = QApplication(sys.argv)
    pl = Plaque('sortiesrealplex/test.txt')
    f = SmartDialog(plaque=pl)
    f.show()
    app.exec_()
