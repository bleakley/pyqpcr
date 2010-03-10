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
        self.setWindowTitle("New target")

    def populateList(self):
        self.listWidget.clear()
        for it in self.project.hashGene.values()[1:]:
            name = "%s (%.2f%%%s%.2f)" % (it.name, it.eff, unichr(177), it.pm)
            item = QListWidgetItem(name)
            if it.isRef == Qt.Checked:
                item.setIcon(QIcon(":/reference.png"))
            item.setStatusTip(it.name)
            self.listWidget.addItem(item)

    def add(self):
        dialog = AddGeneDialog(self, listPlates=self.project.dicoPlates)
        if dialog.exec_():
            nomgene = dialog.gene.text()
            eff = dialog.eff.value()
            pm = dialog.pmerror.value()
            state = dialog.ref.checkState()
            g = Gene(nomgene, eff, pm)
            g.setRef(state)

            if state == Qt.Checked:
# si le gene de ref est pour toutes les plaques, toutes les autres passent a zero
                if dialog.whichPlates.currentText() == QString('All Plates'):
                    for pl in self.project.dicoPlates.values():
                        pl.geneRef = nomgene
                    for gene in self.project.hashGene.values():
                        gene.setRef(Qt.Unchecked)
# sinon, tous les autres genes de cette plaque passent a non ref
                else:
                    currentPlate = dialog.whichPlates.currentText()
                    pl = self.project.dicoPlates[currentPlate]
                    pl.geneRef = nomgene
                    for geneName in pl.dicoGene.keys():
                        gene = self.project.hashGene[geneName]
                        if gene.isRef == Qt.Checked and gene.name != nomgene:
                            gene.setRef(Qt.Unchecked)

            if not self.project.hashGene.has_key(nomgene):
                self.project.hashGene[nomgene] = g
                self.populateList()
            else:
                QMessageBox.warning(self, "Already exist",
                        "The gene %s is already defined !" % nomgene)

    def edit(self):
        if len(self.listWidget.selectedItems()) == 0:
            return
        gene_before = self.listWidget.currentItem().statusTip()
        gene = self.project.hashGene[gene_before]
        dialog = AddGeneDialog(self, ge=gene, listPlates=self.project.dicoPlates)
        if dialog.exec_():
            name = dialog.gene.text()
            eff = dialog.eff.value()
            pm = dialog.pmerror.value()
            state = dialog.ref.checkState()

            gene.setRef(state)
            gene.setEff(eff)
            gene.setPm(pm)
            gene.setName(name)

            if state == Qt.Checked:
# si le gene de ref est pour toutes les plaques, toutes les autres passent a zero
                if dialog.whichPlates.currentText() == QString('All Plates'):
                    for pl in self.project.dicoPlates.values():
                        pl.geneRef = name
                    for g in self.project.hashGene.values():
                        if g.isRef == Qt.Checked and g.name != name:
                            g.setRef(Qt.Unchecked)
# sinon, tous les autres genes de cette plaque passent a non ref
                else:
                    currentPlate = dialog.whichPlates.currentText()
                    pl = self.project.dicoPlates[currentPlate]
                    pl.geneRef = name
                    for geneName in pl.dicoGene.keys():
                        g = self.project.hashGene[geneName]
                        if g.isRef == Qt.Checked and g.name != name:
                            g.setRef(Qt.Unchecked)

# dico
            ind = None
            for pl in self.project.dicoPlates.values():
                if pl.dicoGene.has_key(gene_before) and gene_before != name:
                    ind = pl.dicoGene.index(gene_before)
                    pl.dicoGene.insert(ind, name, pl.dicoGene[gene_before])
                    ind = self.project.hashGene.index(gene_before)

                    for well in pl.dicoGene[name]:
                        well.setGene(gene)
                    pl.dicoGene.__delitem__(gene_before)

            if ind is not None:
                self.project.hashGene.insert(ind, name,
                                            self.project.hashGene[gene_before])
                self.project.hashGene.__delitem__(gene_before)
                self.project.unsaved = True
            self.populateList()

    def remove(self):
        genes = []
        if len(self.listWidget.selectedItems()) == 0:
            return
# Liste des genes a supprimer
        for it in self.listWidget.selectedItems():
            gene = self.project.hashGene[it.statusTip()]
            genes.append(gene)

        reply = QMessageBox.question(self, "Remove",
                        "Remove %s ?" % genes,
                        QMessageBox.Yes|QMessageBox.No)

# Si la suppression est confirmee, on remet tous les puits correspondant
# au gene nul et on supprime l'entree dans hashGene
        if reply == QMessageBox.Yes:
            for gene in genes:
                for pl in self.project.dicoPlates.values():
                    if pl.dicoGene.has_key(gene.name):
                        for well in pl.dicoGene[gene.name]:
                            well.setGene(Gene(''))
                            pl.setDicoGene()

                    if pl.geneRef == gene.name:
                        pl.geneRef = ''
                self.project.hashGene.__delitem__(gene.name)
            self.populateList()
            self.project.unsaved = True

class AddGeneDialog(QDialog):
    
    def __init__(self, parent=None, ge=None, listPlates=None):
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
        self.whichPlates = QComboBox()
        self.whichPlates.addItem("All Plates")
        if listPlates.keys() is not None:
            self.whichPlates.addItems(listPlates.keys())
        if ge is not None:
            self.ref.setCheckState(g.isRef)
        else:
            self.ref.setCheckState(Qt.Unchecked)

        liste = []
        if ge is not None:
            for pl in listPlates.keys():
                if ge.name == listPlates[pl].geneRef:
                    liste.append(pl)
        if len(liste) != 1:
            self.whichPlates.setCurrentIndex(0)
        else:
            self.whichPlates.setCurrentIndex(listPlates.index(liste[0])+1)

        layout = QGridLayout()
        layout.addWidget(lab, 0, 0)
        layout.addWidget(self.gene, 0, 1)
        layout.addWidget(lab2, 1, 0)
        layout.addLayout(hlayout, 1, 1)
        hLay = QHBoxLayout()
        hLay.addWidget(labRef)
        hLay.addWidget(self.ref)
        hLay.addWidget(self.whichPlates)
        layout.addLayout(hLay, 2, 0, 1, 2)
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
