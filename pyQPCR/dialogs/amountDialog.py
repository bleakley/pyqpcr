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

import re
import copy
from pyQPCR.wellGeneSample import *
from PyQt4.QtGui import *
from PyQt4.QtCore import *
import pyQPCR.qrc_resources

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"

class AmountDialog(QDialog):
    
    def __init__(self, parent=None, plaque=None):
        self.parent = parent
        QDialog.__init__(self, parent)

        self.listWidget = QListWidget()
        self.listWidget.setAlternatingRowColors(True)
        if plaque is not None:
            self.plaque = copy.deepcopy(plaque)
            self.populateList()

        self.listWidget.setCurrentRow(-1)
        buttonAdd = QPushButton("&Add")
        buttonEdit = QPushButton("&Edit")
        buttonRemove = QPushButton("&Remove")
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|
                                     QDialogButtonBox.Cancel)

        vlayout = QVBoxLayout()
        vlayout.addWidget(buttonAdd)
        vlayout.addWidget(buttonEdit)
        vlayout.addWidget(buttonRemove)
        vlayout.addStretch()
        hlayout = QHBoxLayout()
        hlayout.addWidget(self.listWidget)
        hlayout.addLayout(vlayout)
        layout = QVBoxLayout()
        layout.addLayout(hlayout)
        layout.addWidget(buttonBox)
        self.setLayout(layout)

        self.connect(buttonBox, SIGNAL("accepted()"), self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"), self, SLOT("reject()"))
        self.connect(buttonAdd, SIGNAL("clicked()"), self.add)
        self.connect(buttonEdit, SIGNAL("clicked()"), self.edit)
        self.connect(buttonRemove, SIGNAL("clicked()"), self.remove)
        self.setWindowTitle("New amount")

    def populateList(self):
        self.listWidget.clear()
        for ind, it in enumerate(self.plaque.listAmount[1:]):
            item = QListWidgetItem(it)
            self.listWidget.addItem(item)

    def add(self):
        dialog = AddAmDialog(self)
        if dialog.exec_():
            am, status = dialog.am.text().toFloat()
            if not self.plaque.adresseAmount.has_key(str(am)):
                self.plaque.listAmount.append(str(am))
                self.plaque.adresseAmount[str(am)] = len(self.plaque.listAmount)-1
                self.populateList()

    def edit(self):
        row = self.listWidget.currentRow()
        am_before = self.plaque.listAmount[row+1]
        dialog = AddAmDialog(self, am=am_before)
        if dialog.exec_():
            am, status = dialog.am.text().toFloat()
            am = str(am)
            self.plaque.listAmount[row+1] = am
            self.populateList()
        if self.plaque.dicoAmount.has_key(str(am_before)):
                self.plaque.dicoAmount[am] = \
                     self.plaque.dicoAmount[str(am_before)]
                self.plaque.adresseAmount[am] = \
                     self.plaque.adresseAmount[str(am_before)]
                for well in self.plaque.dicoAmount[str(am_before)]:
                    well.setAmount(float(am))
                self.plaque.dicoAmount.__delitem__(str(am_before))
                self.plaque.adresseAmount.__delitem__(str(am_before))

    def remove(self):
        row = self.listWidget.currentRow()
        am = self.plaque.listAmount[row+1]
        item = self.listWidget.item(row)
        if item is None:
            return
        reply = QMessageBox.question(self, "Remove",
                        "Remove %s ?" % am,
                        QMessageBox.Yes|QMessageBox.No)
        if reply == QMessageBox.Yes:
            self.plaque.listAmount.__delitem__(row+1)
            self.populateList()
            if self.plaque.dicoAmount.has_key(str(am)):
                for well in self.plaque.dicoAmount[str(am)]:
                    well.setAmount('')

class AddAmDialog(QDialog):
    
    def __init__(self, parent=None, am=None):
        self.parent = parent
        QDialog.__init__(self, parent)
        lab = QLabel("Amount:")
        if am is not None:
            self.am = QLineEdit(am)
        else:
            self.am = QLineEdit()
        self.am.setValidator(QDoubleValidator(self))
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|
                                     QDialogButtonBox.Cancel)

        layout = QGridLayout()
        layout.addWidget(lab, 0, 0)
        layout.addWidget(self.am, 0, 1)
        layout.addWidget(buttonBox, 1, 0, 1, 2)
        self.setLayout(layout)

        self.connect(buttonBox, SIGNAL("accepted()"), self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"), self, SLOT("reject()"))
        self.setWindowTitle("New sample")

if __name__=="__main__":
    app = QApplication(sys.argv)
    pl = Plaque('toto.csv')
    f = AmountDialog(plaque=pl)
    f.show()
    app.exec_()
