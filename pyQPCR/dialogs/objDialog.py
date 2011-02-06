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
from PyQt4.QtGui import *
from PyQt4.QtCore import *

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"

class PropDialog(QDialog):
    
    def __init__(self, parent=None, hashGene=None, hashEch=None):
        self.parent = parent
        QDialog.__init__(self, parent)

        self.geneListWidget = QListWidget()
        self.geneListWidget.setAlternatingRowColors(True)
        pix = QPixmap(32, 32)
        if hashGene is not None:
            self.hashGene = copy.deepcopy(hashGene)
            self.populateListGene()

        self.echListWidget = QListWidget()
        self.echListWidget.setAlternatingRowColors(True)
        pix = QPixmap(32, 32)
        if hashEch is not None:
            self.hashEch = copy.deepcopy(hashEch)
            self.populateListEch()

        buttonUp = QPushButton("&Up")
        buttonDown = QPushButton("&Down")
        buttonColor = QPushButton("&Set color...")
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|
                                     QDialogButtonBox.Cancel)

        vlayout = QVBoxLayout()
        vlayout.addWidget(buttonUp)
        vlayout.addWidget(buttonDown)
        vlayout.addWidget(buttonColor)
        vlayout.addStretch()

        hlayout = QHBoxLayout()
        hlayout.addWidget(self.geneListWidget)
        hlayout.addWidget(self.echListWidget)
        hlayout.addLayout(vlayout)
        layout = QVBoxLayout()
        layout.addLayout(hlayout)
        layout.addWidget(buttonBox)
        self.setLayout(layout)
        self.whois = None

        self.connect(buttonBox, SIGNAL("accepted()"), self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"), self, SLOT("reject()"))
        self.connect(buttonUp, SIGNAL("clicked()"), self.up)
        self.connect(buttonDown, SIGNAL("clicked()"), self.down)
        self.connect(buttonColor, SIGNAL("clicked()"), self.color)
        self.connect(self.geneListWidget, 
                     SIGNAL("itemChanged(QListWidgetItem *)"),
                     self.unselectGene)
        self.connect(self.geneListWidget, SIGNAL("itemSelectionChanged()"),
                     self.unselectEch)
        self.connect(self.geneListWidget, 
                     SIGNAL("itemClicked(QListWidgetItem *)"),
                     self.chooseGeneList)
        self.connect(self.echListWidget, 
                     SIGNAL("itemChanged(QListWidgetItem *)"),
                     self.unselectEch)
        self.connect(self.echListWidget, SIGNAL("itemSelectionChanged()"),
                     self.unselectGene)
        self.connect(self.echListWidget, 
                     SIGNAL("itemClicked(QListWidgetItem *)"),
                     self.chooseEchList)
        self.setWindowTitle("Plot properties")

    def chooseGeneList(self):
        self.whois = "geneList"

    def chooseEchList(self):
        self.whois = "echList"

    def accept(self):
        for ind, it in enumerate(self.hashGene.values()[1:]):
            if self.geneListWidget.item(ind).checkState() == Qt.Unchecked:
                it.setEnabled(False)
            else:
                it.setEnabled(True)
        for ind, it in enumerate(self.hashEch.values()[1:]):
            if self.echListWidget.item(ind).checkState() == Qt.Unchecked:
                it.setEnabled(False)
            else:
                it.setEnabled(True)
        QDialog.accept(self)

    def populateListGene(self):
        pix = QPixmap(32, 32)
        for it in self.hashGene.values()[1:]:
            item = QListWidgetItem(it.name)
            if it.enabled:
                item.setCheckState(Qt.Checked)
            else:
                item.setCheckState(Qt.Unchecked)
            pix.fill(it.color)
            item.setIcon(QIcon(pix))
            self.geneListWidget.addItem(item)

    def populateListEch(self):
        pix = QPixmap(32, 32)
        for it in self.hashEch.values()[1:]:
            item = QListWidgetItem(it.name)
            if it.enabled:
                item.setCheckState(Qt.Checked)
            else:
                item.setCheckState(Qt.Unchecked)
            pix.fill(it.color)
            item.setIcon(QIcon(pix))
            self.echListWidget.addItem(item)

    def unselectEch(self):
        for item in self.echListWidget.selectedItems():
            self.echListWidget.setItemSelected(item, False)
            if item.checkState() == Qt.Unchecked:
                self.hashEch[item.text()].setEnabled(False)
            if item.checkState() == Qt.Checked:
                self.hashEch[item.text()].setEnabled(True)

    def unselectGene(self):
        for item in self.geneListWidget.selectedItems():
            self.geneListWidget.setItemSelected(item, False)
            if item.checkState() == Qt.Unchecked:
                self.hashGene[item.text()].setEnabled(False)
            if item.checkState() == Qt.Checked:
                self.hashGene[item.text()].setEnabled(True)

    def up(self):
        activeGenes = self.geneListWidget.selectedItems()
        activeEchs = self.echListWidget.selectedItems()
        if len(activeGenes) != 0:
            row = self.geneListWidget.currentRow()
            if row >= 1:
                key, popkey = self.hashGene.popitem(row)
                self.hashGene.insert(row+1, key, popkey)
                self.geneListWidget.clear()
                self.populateListGene()
                self.geneListWidget.setCurrentRow(row-1)
        elif len(activeEchs) != 0:
            row = self.echListWidget.currentRow()
            if row >= 1:
                key, popkey = self.hashEch.popitem(row)
                self.hashEch.insert(row+1, key, popkey)
                self.echListWidget.clear()
                self.populateListEch()
                self.echListWidget.setCurrentRow(row-1)

    def down(self):
        activeGenes = self.geneListWidget.selectedItems()
        activeEchs = self.echListWidget.selectedItems()
        if len(activeGenes) != 0:
            row  = self.geneListWidget.currentRow()
            if row < self.geneListWidget.count() -1:
                key, popkey = self.hashGene.popitem(row+2)
                self.hashGene.insert(row+1, key, popkey)
                self.geneListWidget.clear()
                self.populateListGene()
                self.geneListWidget.setCurrentRow(row+1)
        elif len(activeEchs) != 0:
            row  = self.echListWidget.currentRow()
            if row < self.echListWidget.count() -1:
                key, popkey = self.hashEch.popitem(row+2)
                self.hashEch.insert(row+1, key, popkey)
                self.echListWidget.clear()
                self.populateListEch()
                self.echListWidget.setCurrentRow(row+1)

    def color(self):
        pix = QPixmap(32, 32)
        if self.whois == "geneList":
            item = self.geneListWidget.currentItem()
            col = self.hashGene[item.text()].color
            dialog = QColorDialog(self)
            dialog.setCurrentColor(col)
            if dialog.exec_():
                color = dialog.currentColor()
                pix.fill(color)
                item.setIcon(QIcon(pix))
                self.hashGene[item.text()].setColor(color)
        elif self.whois == "echList":
            item = self.echListWidget.currentItem()
            col = self.hashEch[item.text()].color
            dialog = QColorDialog(self)
            dialog.setCurrentColor(col)
            if dialog.exec_():
                color = dialog.currentColor()
                pix.fill(color)
                item.setIcon(QIcon(pix))
                self.hashEch[item.text()].setColor(color)


if __name__=="__main__":
    import sys
    from plate import *
    pl = Plaque("../../treated_data.txt")
    for g in pl.listGene:
        setattr(g, 'color', QColor('#000000'))
    pl.listGene[2].setEnabled(False)
    app = QApplication(sys.argv)
    f = PropDialog(listGene=pl.listGene[1:])
    f.show()
    app.exec_()
