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
        self.cboxType.addItems(["unknown", "standard", "negative", ""])
        dico = {}
        dico['unknown'] = 0
        dico['standard'] = 1
        dico['negative'] = 2
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

#
        self.stackedWidget = QStackedWidget()
#
        sampleWidget = QWidget()
        sampleLayout = QHBoxLayout()
        self.labSample = QLabel("&Sample:")
        self.cboxSample = GeneEchComboBox()
        self.labSample.setBuddy(self.cboxSample)
        sampleLayout.addWidget(self.labSample)
        sampleLayout.addWidget(self.cboxSample)
        sampleLayout.setMargin(0)
        sampleWidget.setLayout(sampleLayout)
        self.stackedWidget.addWidget(sampleWidget)
#
        amountWidget = QWidget()
        amountLayout = QHBoxLayout(amountWidget)
        self.labAmount = QLabel("&Amount:")
        self.editAmount = QLineEdit(amountWidget)
        self.labAmount.setBuddy(self.editAmount)
        amountLayout.addWidget(self.labAmount)
        sizePolicy = QSizePolicy(QSizePolicy.Preferred,
            QSizePolicy.Fixed)
        sizePolicy2 = QSizePolicy(QSizePolicy.Expanding,
            QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy2.setHorizontalStretch(1)
        self.labAmount.setSizePolicy(sizePolicy)
        self.editAmount.setSizePolicy(sizePolicy2)
        amountLayout.addWidget(self.editAmount)
        amountLayout.setMargin(0)
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
        if selected is not None:
            nType = list(selected[0])
            nGene = list(selected[1])
            nEch = list(selected[2])
# Determination de l'item courant pour le type
            if len(nType) == 1:
                self.cboxType.setCurrentIndex(dico[str(nType[0])])
            else:
                self.cboxType.setCurrentIndex(3)
# Determination de l'item courant pour le gene
            if len(nGene) == 1:
                ind = self.plaque.adresseGene[str(nGene[0])]
                self.cboxGene.setCurrentIndex(ind)
            else:
                self.cboxGene.setCurrentIndex(0)
# Determination de l'item courant pour l'echantillon
            if len(nEch) == 1:
                ind = self.plaque.adresseEch[str(nEch[0])]
                self.cboxSample.setCurrentIndex(ind)
            else:
                self.cboxSample.setCurrentIndex(0)


        topLayout = QGridLayout()
        topLayout.addWidget(lab1, 0, 0)
        topLayout.addWidget(self.cboxType, 0, 1)
        topLayout.addWidget(lab2, 1, 0)
        topLayout.addWidget(self.cboxGene, 1, 1)
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
        self.setWindowTitle("Edit")
        self.setMinimumSize(280,100)

    def modifDialog(self):
        if self.cboxType.currentText() in ("unknown", "negative"):
            self.stackedWidget.setCurrentIndex(0)
        if self.cboxType.currentText() == "standard":
            self.stackedWidget.setCurrentIndex(1)

    def accept(self):
        ge = self.cboxGene.currentObj()
        if self.cboxType.currentText() == "unknown":
            ech = self.cboxSample.currentObj()
            for it in self.parent.table.selectedItems():
                nom = str(it.statusTip())
                well = getattr(self.plaque, nom)
                well.setType(QString("unknown"))
                if ge.name != "": well.setGene(ge)
                if ech.name != "": well.setEch(ech)

        if self.cboxType.currentText() == "standard":
            try:
                am = float(self.editAmount.text())
            except ValueError, e:
                QMessageBox.warning(self, "Error", str(e))
                return
            for it in self.parent.table.selectedItems():
                nom = str(it.statusTip())
                well = getattr(self.plaque, nom)
                well.setType(QString("standard"))
                if ge.name != '': well.setGene(ge)
                well.setAmount(am)
        QDialog.accept(self)


if __name__=="__main__":
    import sys
    from plaque import *
    app = QApplication(sys.argv)
    pl = Plaque('sortiesrealplex/test_2.txt')
    f = EditDialog(plaque=pl)
    f.show()
    app.exec_()
