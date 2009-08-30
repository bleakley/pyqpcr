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
from pyQPCR.wellGeneSample import *
from pyQPCR.dialogs.customWidgets import *
from pyQPCR.dialogs.echDialog import *
from pyQPCR.dialogs.geneDialog import *
from pyQPCR.dialogs.amountDialog import *
from PyQt4.QtGui import *
from PyQt4.QtCore import *

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"

class EditDialog(QDialog):

    def __init__(self, parent=None, plaque=None, selected=None):
        self.parent = parent
        QDialog.__init__(self, parent)
# Widgets
        lab1 = QLabel("T&ype:")
        self.cboxType = QComboBox()
        lab1.setBuddy(self.cboxType)
        self.cboxType.addItems(['unknown', 'standard', 'negative', ''])
        dico = {}
        dico[QString('unknown')] = 0
        dico[QString('standard')] = 1
        dico[QString('negative')] = 2
        pix = QPixmap(32, 32)
        bleu = QColor(116, 167, 227) ; pix.fill(bleu)
        self.cboxType.setItemIcon(0, QIcon(pix))
        rouge = QColor(233, 0, 0) ; pix.fill(rouge)
        self.cboxType.setItemIcon(1, QIcon(pix))
        jaune = QColor(255, 250, 80) ; pix.fill(jaune)
        self.cboxType.setItemIcon(2, QIcon(pix))
        lab2 = QLabel("&Target:")
        self.cboxGene = GeneEchComboBox()
        lab2.setBuddy(self.cboxGene)
        ic = QIcon(":/addgene")
        btnAddEch = QToolButton()
        btnAddEch.setIcon(ic)
        btnAddGene = QToolButton()
        btnAddGene.setIcon(ic)
        btnAddAm = QToolButton()
        btnAddAm.setIcon(ic)
#
        self.stackedWidget = QStackedWidget()
#
        sampleWidget = QWidget()
        sampleLayout = QHBoxLayout()
        self.labSample = QLabel("&Sample:")
        self.cboxSample = GeneEchComboBox()
        self.labSample.setBuddy(self.cboxSample)
        sampleLayout.addWidget(self.labSample)
        hLay = QHBoxLayout()
        hLay.addWidget(self.cboxSample)
        hLay.addWidget(btnAddEch)
        sampleLayout.addLayout(hLay)
        sampleWidget.setLayout(sampleLayout)
        self.stackedWidget.addWidget(sampleWidget)
#
        amountWidget = QWidget()
        amountLayout = QHBoxLayout(amountWidget)
        self.labAmount = QLabel("&Amount:")
        self.cboxAm = QComboBox()
        self.labAmount.setBuddy(self.cboxAm)
        hLay = QHBoxLayout()
        amountLayout.addWidget(self.labAmount)
        hLay.addWidget(self.cboxAm)
        hLay.addWidget(btnAddAm)
        amountLayout.addLayout(hLay)
        amountWidget.setLayout(amountLayout)
        self.stackedWidget.addWidget(amountWidget)
#
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|
                                     QDialogButtonBox.Cancel)

# Remplissage de la comboBox avec les genes
        if plaque is not None:
            self.plaque = copy.deepcopy(plaque)
            self.cboxGene.addItems(self.plaque.listGene)
            self.cboxSample.addItems(self.plaque.listEch)
            self.cboxAm.addItems(self.plaque.listAmount)
        if selected is not None:
            self.selected = selected
            self.populateEch()
            self.populateGene()
            self.populateAm()

            nType = list(self.selected[0])
# Determination de l'item courant pour le type
            if len(nType) == 1:
                self.cboxType.setCurrentIndex(dico[nType[0]])
            else:
                self.cboxType.setCurrentIndex(3)

        topLayout = QGridLayout()
        topLayout.addWidget(lab1, 0, 0)
        topLayout.addWidget(self.cboxType, 0, 1)
        topLayout.addWidget(lab2, 1, 0)
        hLay = QHBoxLayout()
        hLay.addWidget(self.cboxGene)
        hLay.addWidget(btnAddGene)
        topLayout.addLayout(hLay, 1, 1)
# Layout
        layout = QVBoxLayout()
        layout.addLayout(topLayout)
        layout.addWidget(self.stackedWidget)
        layout.addWidget(buttonBox)
        self.setLayout(layout)

        self.modifDialog()
# Connections
        self.connect(buttonBox, SIGNAL("accepted()"), self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"), self, SLOT("reject()"))
        self.connect(self.cboxType, SIGNAL("activated(int)"), self.modifDialog)
        self.connect(btnAddEch, SIGNAL("clicked()"), self.addEch)
        self.connect(btnAddAm, SIGNAL("clicked()"), self.addAm)
        self.connect(btnAddGene, SIGNAL("clicked()"), self.addGene)
        self.setWindowTitle("Edit")
        self.setMinimumSize(280,100)

    def modifDialog(self):
        if self.cboxType.currentText() in (QString('unknown'), QString('negative')):
            self.stackedWidget.setCurrentIndex(0)
        if self.cboxType.currentText() == QString('standard'):
            self.stackedWidget.setCurrentIndex(1)

    def populateAm(self):
        self.cboxAm.clear()
        self.cboxAm.addItems(self.plaque.listAmount)
        nAm = list(self.selected[3])
# Determination de l'item courant pour l'amount
        if len(nAm) == 1:
            try:
                ind = self.plaque.adresseAmount[str(nAm[0])]
                self.cboxAm.setCurrentIndex(ind)
                self.curAmIndex = ind
            except KeyError:
                self.cboxAm.setCurrentIndex(self.curAmIndex)
        else:
            self.cboxAm.setCurrentIndex(0)

    def populateEch(self):
        self.cboxSample.clear()
        self.cboxSample.addItems(self.plaque.listEch)
        nEch = list(self.selected[2])
# Determination de l'item courant pour l'echantillon
        if len(nEch) == 1:
            try:
                ind = self.plaque.adresseEch[nEch[0]]
                self.cboxSample.setCurrentIndex(ind)
                self.curEchIndex = ind
            except KeyError:
                self.cboxSample.setCurrentIndex(self.curEchIndex)
        else:
            self.cboxSample.setCurrentIndex(0)

    def populateGene(self):
        self.cboxGene.clear()
        self.cboxGene.addItems(self.plaque.listGene)
        nGene = list(self.selected[1])
# Determination de l'item courant pour le gene
        if len(nGene) == 1:
            try:
                ind = self.plaque.adresseGene[nGene[0]]
                self.cboxGene.setCurrentIndex(ind)
                self.curGeneIndex = ind
            except KeyError:
                self.cboxGene.setCurrentIndex(self.curGeneIndex)
        else:
            self.cboxGene.setCurrentIndex(0)

    def addEch(self):
        dialog = EchDialog(self, plaque=self.plaque)
        if dialog.exec_():
            self.plaque = dialog.plaque
            self.plaque.setDicoEch()
            self.populateEch()

    def addAm(self):
        dialog = AmountDialog(self, plaque=self.plaque)
        if dialog.exec_():
            self.plaque = dialog.plaque
            self.plaque.setDicoAm()
            self.populateAm()

    def addGene(self):
        dialog = GeneDialog(self, plaque=self.plaque)
        if dialog.exec_():
            self.plaque = dialog.plaque
            self.plaque.setDicoGene()
            self.populateGene()

    def accept(self):
        ge = self.cboxGene.currentObj()
        if self.cboxType.currentText() == QString('unknown'):
            ech = self.cboxSample.currentObj()
            for it in self.parent.table.selectedItems():
                nom = it.statusTip()
                well = getattr(self.plaque, str(nom))
                well.setType(QString('unknown'))
                if ge.name != "": well.setGene(ge)
                if ech.name != "": well.setEch(ech)

        if self.cboxType.currentText() == QString('standard'):
            ind = self.cboxAm.currentIndex()
            am = self.plaque.listAmount[ind]
            for it in self.parent.table.selectedItems():
                nom = it.statusTip()
                well = getattr(self.plaque, str(nom))
                well.setType(QString('standard'))
                if ge.name != '':
                    well.setGene(ge)
                if am != '':
                    well.setAmount(float(am))
        QDialog.accept(self)


if __name__=="__main__":
    import sys
    from plaque import *
    app = QApplication(sys.argv)
    pl = Plaque('sortiesrealplex/test_2.txt')
    f = EditDialog(plaque=pl)
    f.show()
    app.exec_()
