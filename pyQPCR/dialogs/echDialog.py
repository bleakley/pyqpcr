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
    
    def __init__(self, parent=None, project=None):
        self.parent = parent
        QDialog.__init__(self, parent)

        self.listWidget = QListWidget()
        self.listWidget.setAlternatingRowColors(True)
        self.listWidget.setSelectionMode(QAbstractItemView.ExtendedSelection)

        if project is not None:
            self.project = copy.deepcopy(project)
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
        for it in self.project.hashEch.values()[1:]:
            item = QListWidgetItem(it.name)
            if it.isRef == Qt.Checked:
                item.setIcon(QIcon(":/reference.png"))
            item.setStatusTip(it.name)
            self.listWidget.addItem(item)

    def add(self):
        dialog = AddEchDialog(self, listPlates=self.project.dicoPlates)
        if dialog.exec_():
            nomech = dialog.ech.text()
            state = dialog.ref.checkState()
            e = Ech(nomech, state)

            if state == Qt.Checked:
                if dialog.whichPlates.currentText() == QString('All Plates'):
                    for pl in self.project.dicoPlates.values():
                        pl.echRef = nomech
                    for ech in self.project.hashEch.values():
                        ech.setRef(Qt.Unchecked)

                else:
                    currentPlate = dialog.whichPlates.currentText()
                    pl = self.project.dicoPlates[currentPlate]
                    pl.echRef = nomech
                    for echName in pl.dicoEch.keys():
                        ech = self.project.hashEch[echName]
                        if ech.isRef == Qt.Checked and ech.name != nomech:
                            ech.setRef(Qt.Unchecked)

            if not self.project.hashEch.has_key(nomech):
                self.project.hashEch[nomech] = e
                self.populateList()
            else:
                QMessageBox.warning(self, "Already exist",
                        "The sample %s is already defined !" % e)

    def edit(self):
        if len(self.listWidget.selectedItems()) == 0:
            return
        ech_before = self.listWidget.currentItem().statusTip()
        ech = self.project.hashEch[ech_before]
        dialog = AddEchDialog(self, ech=ech, listPlates=self.project.dicoPlates)
        if dialog.exec_():
            name = dialog.ech.text()
            state = dialog.ref.checkState()
# Si l'echantillon etait ech de reference et qu'il est desactive
# alors la project n'a plus d'echantillon de reference
            if ech.isRef == Qt.Checked and state == Qt.Unchecked:
                delattr(self.project, "echRef")

            ech.setName(name)
            ech.setRef(state)

            if state == Qt.Checked:
                if dialog.whichPlates.currentText() == QString('All Plates'):
                    for pl in self.project.dicoPlates.values():
                        pl.echRef = name
                    for e in self.project.hashEch.values():
                        if e.isRef == Qt.Checked and e.name != name:
                            e.setRef(Qt.Unchecked)

                else:
                    currentPlate = dialog.whichPlates.currentText()
                    pl = self.project.dicoPlates[currentPlate]
                    pl.echRef = name
                    for echName in pl.dicoEch.keys():
                        e = self.project.hashEch[echName]
                        if e.isRef == Qt.Checked and e.name != name:
                            e.setRef(Qt.Unchecked)
# dico
            ind = None
            for pl in self.project.dicoPlates.values():
                if pl.dicoEch.has_key(ech_before) and ech_before != name:
                    ind = pl.dicoEch.index(ech_before)
                    pl.dicoEch.insert(ind, name, pl.dicoEch[ech_before])
                    ind = self.project.hashEch.index(ech_before)

                    for well in pl.dicoEch[name]:
                        well.setEch(ech)
                    pl.dicoEch.__delitem__(ech_before)

            if ind is not None:
                self.project.hashEch.insert(ind, name,
                                            self.project.hashEch[ech_before])
                self.project.hashEch.__delitem__(ech_before)
                self.project.unsaved = True
            self.populateList()

    def remove(self):
        echs = []
        if len(self.listWidget.selectedItems()) == 0:
            return
        for it in self.listWidget.selectedItems():
            ech = self.project.hashEch[it.statusTip()]
            echs.append(ech)

        reply = QMessageBox.question(self, "Remove",
                        "Remove %s ?" % echs,
                        QMessageBox.Yes|QMessageBox.No)
        if reply == QMessageBox.Yes:
            for ech in echs:
                for pl in self.project.dicoPlates.values():
                    if pl.dicoEch.has_key(ech.name):
                        for well in pl.dicoEch[ech.name]:
                            well.setEch(Ech(''))
                        pl.setDicoEch()

                    if pl.echRef == ech.name:
                        pl.echRef = ''
                self.project.hashEch.__delitem__(ech.name)
            self.project.unsaved = True
            self.populateList()

class AddEchDialog(QDialog):
    
    def __init__(self, parent=None, ech=None, listPlates=None):
        self.parent = parent
        QDialog.__init__(self, parent)
        lab = QLabel("Sample:")
        if ech is not None:
            self.ech = QLineEdit(ech.name)
        else:
            self.ech = QLineEdit()
        labRef = QLabel("Reference:")
        self.ref = QCheckBox()
        self.whichPlates = QComboBox()
        self.whichPlates.addItem("All Plates")
        if listPlates.keys() is not None:
            self.whichPlates.addItems(listPlates.keys())

        if ech is not None:
            self.ref.setCheckState(ech.isRef)
        else:
            self.ref.setCheckState(Qt.Unchecked)
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|
                                     QDialogButtonBox.Cancel)

        liste = []
        if ech is not None:
            for pl in listPlates.keys():
                if ech.name == listPlates[pl].echRef:
                    liste.append(pl)
        if len(liste) != 1:
            self.whichPlates.setCurrentIndex(0)
        else:
            self.whichPlates.setCurrentIndex(listPlates.index(liste[0])+1)


        layout = QGridLayout()
        layout.addWidget(lab, 0, 0)
        layout.addWidget(self.ech, 0, 1)
        hLay = QHBoxLayout()
        hLay.addWidget(labRef)
        hLay.addWidget(self.ref)
        hLay.addWidget(self.whichPlates)
        layout.addLayout(hLay, 1, 0, 1, 2)
        layout.addWidget(buttonBox, 2, 0, 1, 2)
        self.setLayout(layout)

        self.connect(buttonBox, SIGNAL("accepted()"), self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"), self, SLOT("reject()"))
        self.setWindowTitle("New sample")

if __name__ == "__main__":
    import sys
    from project import *
    app = QApplication(sys.argv)
    pl = project('toto.csv')
    f = EchDialog(plaque=pl)
    f.show()
    app.exec_()
