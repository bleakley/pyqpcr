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
        """
        Constructor
        """
        self.parent = parent
        QDialog.__init__(self, parent)

        self.listWidget = QListWidget()
        self.listWidget.setAlternatingRowColors(True)
        self.listWidget.setSelectionMode(QAbstractItemView.ExtendedSelection)
        if project is not None:
            self.project = copy.deepcopy(project)
            self.populateList()

        self.colors = [QColor(Qt.blue), QColor(Qt.red), QColor(Qt.green), 
                  QColor(Qt.yellow), QColor(Qt.magenta),
                  QColor(Qt.cyan), QColor(Qt.gray),
                  QColor(Qt.darkBlue), QColor(Qt.darkRed),
                  QColor(Qt.darkGreen), QColor(Qt.darkYellow),
                  QColor(Qt.darkMagenta), QColor(Qt.darkCyan),
                  QColor(Qt.darkGray), QColor(Qt.lightGray), 
                  QColor(Qt.black)]
        self.indexColor = 0

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
        """
        This methods allows to populate the QListWidget which contains the different
        targets of an experiment.
        """
        self.listWidget.clear()
        for it in self.project.hashGene.values()[1:]:
            name = "%s (%.2f%%%s%.2f)" % (it.name, it.eff, unichr(177), it.pm)
            item = QListWidgetItem(name)
            if it.isRef == Qt.Checked:
                item.setIcon(QIcon(":/reference.png"))
            item.setStatusTip(it.name)
            self.listWidget.addItem(item)

    def add(self):
        """
        This method allows to add a new Gene to your experiment. You cannot enter
        a name which already exists
        """
        dialog = AddGeneDialog(self, listPlates=self.project.dicoPlates)
        if dialog.exec_():
            nomgene = dialog.gene.text()
            eff = dialog.eff.value()
            pm = dialog.pmerror.value()
            state = dialog.ref.checkState()
            g = Gene(nomgene, eff, pm)
            g.setRef(state)
            # color attributions
            g.setColor(self.colors[self.indexColor])
            if self.indexColor < len(self.colors)-1:
                self.indexColor += 1
            else:
                self.indexColor = 0

            if not self.project.hashGene.has_key(nomgene):
                self.project.hashGene[nomgene] = g
                for plaque in self.project.dicoPlates.keys():
                    pl = self.project.dicoPlates[plaque]
                    pl.dicoGene[nomgene] = []
            else:
                QMessageBox.warning(self, "Already exist",
                        "The gene <b>%s</b> is already defined !" % nomgene)
                return

            if state == Qt.Checked:
# si le gene de ref est pour toutes les plaques, toutes les autres passent a zero
                for pl in dialog.refPlates:
                    plaque = self.project.dicoPlates[pl]
                    plaque.geneRef.append(nomgene)
            self.populateList()

    def edit(self):
        """
        This method allows to change the properties of a gene. You cannot change 
        its name in an existing one.
        """
        if len(self.listWidget.selectedItems()) == 0:
            return
        gene_before = self.listWidget.currentItem().statusTip()
        gene = copy.deepcopy(self.project.hashGene[gene_before])
        dialog = AddGeneDialog(self, ge=gene, listPlates=self.project.dicoPlates)
        if dialog.exec_():
            name = dialog.gene.text()
            eff = dialog.eff.value()
            pm = dialog.pmerror.value()
            state = dialog.ref.checkState()
            if self.project.hashGene.has_key(name) and gene_before != name:
                QMessageBox.warning(self, "Already exist",
                        "The gene <b>%s</b> is already defined !" % name)
                return
# Si le gene etait gene de reference et qu'il est desactive
# alors la plaque a un gene de ref de moins
            if gene.isRef == Qt.Checked and state == Qt.Unchecked:
                for pl in self.project.dicoPlates.values():
                    for geneName in pl.geneRef:
                        if geneName == gene.name or geneName == gene_before:
                            pl.geneRef.remove(gene.name)

            gene.setRef(state)
            gene.setEff(eff)
            gene.setPm(pm)
            gene.setName(name)

            if state == Qt.Checked:
# si le gene de ref est pour toutes les plaques, toutes les autres passent a zero
                for pl in dialog.refPlates:
                    plaque = self.project.dicoPlates[pl]
                    if not plaque.geneRef.__contains__(name):
                        if gene_before != name and plaque.geneRef.__contains__(gene_before):
                            plaque.geneRef.remove(gene_before)
                        plaque.geneRef.append(name)
# dico
            ind = None
            for plaque in self.project.dicoPlates.keys():
                pl = self.project.dicoPlates[plaque]
                if plaque not in dialog.refPlates:
                    for gg in pl.geneRef:
                        if gg == name or gg == gene_before:
                            pl.geneRef.remove(gg)

                if pl.dicoGene.has_key(gene_before) and \
                               gene != self.project.hashGene[gene_before]:
                    ind = pl.dicoGene.index(gene_before)
                    pl.dicoGene.insert(ind, name, pl.dicoGene[gene_before])
                    ind = self.project.hashGene.index(gene_before)

                    for well in pl.dicoGene[name]:
                        well.setGene(gene)
                    if gene_before != gene.name:
                        pl.dicoGene.__delitem__(gene_before)

            if ind is not None:
                self.project.hashGene.insert(ind, name, gene)
                if gene_before != gene.name:
                    self.project.hashGene.__delitem__(gene_before)
                self.project.unsaved = True
            self.populateList()

    def remove(self):
        """
        This method deletes an existing gene from your experiment. 
        Every well containing this gene are now on the empty '' gene.
        """
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

                    for geneName in pl.geneRef:
                        if geneName == gene.name:
                            pl.geneRef.remove(gene.name)
                self.project.hashGene.__delitem__(gene.name)
            self.populateList()
            self.project.unsaved = True

class AddGeneDialog(QDialog):
    
    def __init__(self, parent=None, ge=None, listPlates=None):
        self.parent = parent
        self.listPlates = listPlates
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

        self.widList = QListWidget()
        if len(self.listPlates.keys()) > 1:
            self.widList.setVisible(self.ref.isChecked())
        else:
            self.widList.setVisible(False)
        self.widList.setAlternatingRowColors(True)
        self.widList.setSelectionMode(QAbstractItemView.MultiSelection)
        self.populateList(ge)

        layout = QGridLayout()
        layout.addWidget(lab, 0, 0)
        layout.addWidget(self.gene, 0, 1)
        layout.addWidget(lab2, 1, 0)
        layout.addLayout(hlayout, 1, 1)
        hLay = QHBoxLayout()
        hLay.addWidget(labRef)
        hLay.addWidget(self.ref)
        layout.setSizeConstraint(QLayout.SetFixedSize)
        layout.addLayout(hLay, 2, 0, 1, 2)
        layout.addWidget(self.widList, 3, 0, 1, 2)
        layout.addWidget(buttonBox, 4, 0, 1, 2)
        self.setLayout(layout)

        self.connect(buttonBox, SIGNAL("accepted()"), self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"), self, SLOT("reject()"))
        self.connect(self.ref, SIGNAL("stateChanged(int)"), self.unHide)
        self.setWindowTitle("New target")

    def populateList(self, ge):
        self.widList.clear()
        for it in self.listPlates.keys():
            item = QListWidgetItem(it)
            item.setStatusTip(it)
            self.widList.addItem(item)
            if ge is not None:
                if ge.name in self.listPlates[it].geneRef:
                    self.widList.setItemSelected(item, True)

    def unHide(self):
        if len(self.listPlates.keys()) == 1:
            self.widList.clear()
            item = QListWidgetItem(self.listPlates.keys()[0])
            self.widList.addItem(item)
            self.widList.setItemSelected(item, True)
        else:
            self.widList.setVisible(self.ref.isChecked())

    def accept(self):
        self.refPlates = []
        if len(self.widList.selectedItems()) != 0:
            for it in self.widList.selectedItems():
                self.refPlates.append(it.text())
        else:
            self.ref.setCheckState(Qt.Unchecked)

        QDialog.accept(self)        


if __name__=="__main__":
    import sys
    from project import *
    app = QApplication(sys.argv)
    pl = Project('../../samples/2plates.xml')
    f = GeneDialog(project=pl)
    f.show()
    app.exec_()
