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
from PyQt4.QtGui import *
from PyQt4.QtCore import *
import pyQPCR.qrc_resources

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"

class GeneDialog(QDialog):
    
    def __init__(self, parent=None, plaque=None):
        self.parent = parent
        QDialog.__init__(self, parent)

        self.listWidget = QListWidget()
        self.listWidget.setAlternatingRowColors(True)
        self.listWidget.setSelectionMode(3)
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
        self.setWindowTitle("New target")

    def populateList(self):
        self.listWidget.clear()
        for ind, it in enumerate(self.plaque.listGene[1:]):
            name = "%s (%.2f%%%s%.2f)" % (it.name, it.eff, unichr(177), it.pm)
            item = QListWidgetItem(name)
            if it.isRef == Qt.Checked:
                item.setIcon(QIcon(":/reference.png"))
            item.setStatusTip(it.name)
            self.listWidget.addItem(item)

    def add(self):
        dialog = AddGeneDialog(self)
        if dialog.exec_():
            gene = dialog.gene.text()
            eff = dialog.eff.value()
            pm = dialog.pmerror.value()
            state = dialog.ref.checkState()
            g = Gene(gene, eff, pm)
            g.setRef(state)
# Si le gene ajoute est un gene de reference alors l'autre gene de reference
# repasse en isRef = Qt.Unchecked
            if state == Qt.Checked:
                self.plaque.geneRef = g
                for gene in self.plaque.listGene:
                    if gene.isRef == Qt.Checked:
                        gene.setRef(Qt.Unchecked)
            g.setColor(QColor(Qt.black))
            if not self.plaque.adresseGene.has_key(gene):
                self.plaque.listGene.append(g)
                self.plaque.adresseGene[gene] = len(self.plaque.listGene)-1
                self.populateList()

    def edit(self):
        row = self.listWidget.currentRow()
        gene = self.plaque.listGene[row+1]
        gene_before = gene.name
        dialog = AddGeneDialog(self, ge=gene)
        if dialog.exec_():
            name = dialog.gene.text()
            eff = dialog.eff.value()
            pm = dialog.pmerror.value()
            state = dialog.ref.checkState()
            gene.setRef(state)
            gene.setEff(eff)
            gene.setPm(pm)
            gene.setName(name)
# Si le gene ajoute est un gene de reference alors l'autre gene de reference
# repasse en isRef = Qt.Unchecked
            if state == Qt.Checked:
                self.plaque.geneRef = gene
                for ind, g in enumerate(self.plaque.listGene):
                    if g.isRef == Qt.Checked and ind != row+1:
                        g.setRef(Qt.Unchecked)
# dico
            self.populateList()
            if self.plaque.dicoGene.has_key(gene_before):
                self.plaque.dicoGene[name] = \
                     self.plaque.dicoGene[gene_before]
                self.plaque.adresseGene[name] = \
                     self.plaque.adresseGene[gene_before]
                for well in self.plaque.dicoGene[gene_before]:
                    well.setGene(gene)
                self.plaque.dicoGene.__delitem__(gene_before)
                self.plaque.adresseGene.__delitem__(gene_before)
                self.plaque.unsaved = True

    def remove(self):
        genes = []
        if len(self.listWidget.selectedItems()) == 0:
            return
        for it in self.listWidget.selectedItems():
            row = self.plaque.adresseGene[it.statusTip()]
            gene = self.plaque.listGene[row]
            genes.append(gene)

        reply = QMessageBox.question(self, "Remove",
                        "Remove %s ?" % genes,
                        QMessageBox.Yes|QMessageBox.No)
        if reply == QMessageBox.Yes:
            for gene in genes:
                if self.plaque.dicoGene.has_key(gene.name):
                    for well in self.plaque.dicoGene[gene.name]:
                        well.setGene(Gene(''))
                    self.plaque.setDicoGene()
                    self.plaque.listGene.__delitem__( \
                                  self.plaque.adresseGene[gene.name])
                self.populateList()
                self.plaque.unsaved = True

class AddGeneDialog(QDialog):
    
    def __init__(self, parent=None, ge=None):
        self.parent = parent
        QDialog.__init__(self, parent)
        lab = QLabel("Target:")
        if ge is not None:
            g = copy.deepcopy(ge)
            self.gene = QLineEdit(g.name)
        else:
            self.gene = QLineEdit()
        lab2 = QLabel("Efficiency:")
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|
                                     QDialogButtonBox.Cancel)
        hlayout = QHBoxLayout()
        self.eff = QDoubleSpinBox()
# Pour changer les , par des . on force la locale
        self.eff.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.eff.setRange(0.0, 120.0)
        self.eff.setSuffix(" %")
        if ge is not None:
            self.eff.setValue(g.eff)
        else:
            self.eff.setValue(100.0)
        self.pmerror = QDoubleSpinBox()
# Pour changer les , par des . on force la locale
        self.pmerror.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.pmerror.setRange(0.0, 100.0)
        self.pmerror.setSuffix(" %")
        if ge is not None:
            self.pmerror.setValue(g.pm)
        else:
            self.pmerror.setValue(0.)
        hlayout.addWidget(self.eff)
        hlayout.addWidget(QLabel(unichr(177)))
        hlayout.addWidget(self.pmerror)
        labRef = QLabel("Reference:")
        self.ref = QCheckBox()
        if ge is not None:
            self.ref.setCheckState(g.isRef)
        else:
            self.ref.setCheckState(Qt.Unchecked)

        layout = QGridLayout()
        layout.addWidget(lab, 0, 0)
        layout.addWidget(self.gene, 0, 1)
        layout.addWidget(lab2, 1, 0)
        layout.addLayout(hlayout, 1, 1)
        layout.addWidget(labRef, 2, 0)
        layout.addWidget(self.ref, 2, 1)
        layout.addWidget(buttonBox, 3, 0, 1, 2)
        self.setLayout(layout)

        self.connect(buttonBox, SIGNAL("accepted()"), self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"), self, SLOT("reject()"))
        self.setWindowTitle("New target")


if __name__=="__main__":
    import sys
    from plaque import *
    app = QApplication(sys.argv)
    pl = Plaque('toto.csv')
    f = GeneDialog(plaque=pl)
    f.show()
    app.exec_()
