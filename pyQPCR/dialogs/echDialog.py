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

class EchDialog(QDialog):
    
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
        self.setWindowTitle("New sample")

    def populateList(self):
        self.listWidget.clear()
        for ind, it in enumerate(self.plaque.listEch[1:]):
            item = QListWidgetItem(it.name)
            if it.isRef == Qt.Checked:
                item.setIcon(QIcon(":/reference.png"))
            self.listWidget.addItem(item)

    def add(self):
        dialog = AddEchDialog(self)
        if dialog.exec_():
            nomech = dialog.ech.text()
            state = dialog.ref.checkState()
            ech = Ech(nomech, state)
            ech.setColor(QColor(Qt.black))
            if state == Qt.Checked:
                self.plaque.echRef = ech
                for sample in self.plaque.listEch:
                    if sample.isRef == Qt.Checked:
                        sample.setRef(Qt.Unchecked)
            if not self.plaque.adresseEch.has_key(nomech):
                self.plaque.listEch.append(ech)
                self.plaque.adresseEch[nomech] = len(self.plaque.listEch)-1
                self.populateList()

    def edit(self):
        row = self.listWidget.currentRow()
        ech = self.plaque.listEch[row+1]
        ech_before = ech.name
        dialog = AddEchDialog(self, ech=ech)
        if dialog.exec_():
            name = dialog.ech.text()
            state = dialog.ref.checkState()
            ech.setName(name)
            ech.setRef(state)
            if state == Qt.Checked:
                self.plaque.echRef = ech
                for ind, sample in enumerate(self.plaque.listEch):
                    if sample.isRef == Qt.Checked and ind != row+1:
                        sample.setRef(Qt.Unchecked)
            self.populateList()
            if self.plaque.dicoEch.has_key(ech_before):
                self.plaque.dicoEch[name] = \
                     self.plaque.dicoEch[ech_before]
                self.plaque.adresseEch[name] = \
                     self.plaque.adresseEch[ech_before]
                for well in self.plaque.dicoEch[ech_before]:
                    well.setEch(ech)
                self.plaque.dicoEch.__delitem__(ech_before)
                self.plaque.adresseEch.__delitem__(ech_before)

    def remove(self):
        row = self.listWidget.currentRow()
        ech = self.plaque.listEch[row+1]
        item = self.listWidget.item(row)
        if item is None:
            return
        reply = QMessageBox.question(self, "Remove",
                        "Remove %s ?" % ech,
                        QMessageBox.Yes|QMessageBox.No)
        if reply == QMessageBox.Yes:
            self.plaque.listEch.__delitem__(row+1)
            self.populateList()
# On ajuste les puits de la plaque concernes sur l'echantillon vide
# a condition que le puit soit concerne
            if self.plaque.dicoEch.has_key(QString(ech.name)):
                for well in self.plaque.dicoEch[QString(ech.name)]:
                    well.setEch(Ech(''))

class AddEchDialog(QDialog):
    
    def __init__(self, parent=None, ech=None):
        self.parent = parent
        QDialog.__init__(self, parent)
        lab = QLabel("Sample:")
        if ech is not None:
            self.ech = QLineEdit(ech.name)
        else:
            self.ech = QLineEdit()
        labRef = QLabel("Reference:")
        self.ref = QCheckBox()
        if ech is not None:
            self.ref.setCheckState(ech.isRef)
        else:
            self.ref.setCheckState(Qt.Unchecked)
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|
                                     QDialogButtonBox.Cancel)

        layout = QGridLayout()
        layout.addWidget(lab, 0, 0)
        layout.addWidget(self.ech, 0, 1)
        layout.addWidget(labRef, 1, 0)
        layout.addWidget(self.ref, 1, 1)
        layout.addWidget(buttonBox, 2, 0, 1, 2)
        self.setLayout(layout)

        self.connect(buttonBox, SIGNAL("accepted()"), self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"), self, SLOT("reject()"))
        self.setWindowTitle("New sample")

if __name__=="__main__":
    import sys
    from plaque import *
    app = QApplication(sys.argv)
    pl = Plaque('toto.csv')
    f = EchDialog(plaque=pl)
    f.show()
    app.exec_()
