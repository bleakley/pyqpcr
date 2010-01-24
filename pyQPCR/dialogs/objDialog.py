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
from pyQPCR.widgets.customCbox import *
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
            self.listGene = copy.deepcopy(hashGene.values()[1:])
            self.populateListGene()

        self.echListWidget = QListWidget()
        self.echListWidget.setAlternatingRowColors(True)
        pix = QPixmap(32, 32)
        if hashEch is not None:
            self.listEch = copy.deepcopy(hashEch.values()[1:])
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
        self.connect(self.geneListWidget, SIGNAL("itemSelectionChanged()"),
                     self.unselectEch)
        self.connect(self.geneListWidget, 
                     SIGNAL("itemActivated(QListWidgetItem *)"),
                     self.chooseGeneList)
        self.connect(self.echListWidget, SIGNAL("itemSelectionChanged()"),
                     self.unselectGene)
        self.connect(self.echListWidget, 
                     SIGNAL("itemActivated(QListWidgetItem *)"),
                     self.chooseEchList)
        self.setWindowTitle("Plot properties")

    def chooseGeneList(self):
        self.whois = "geneList"

    def chooseEchList(self):
        self.whois = "echList"

    def accept(self):
        for ind, it in enumerate(self.listGene):
            it.setEnabled(self.geneListWidget.item(ind).checkState())
        for ind, it in enumerate(self.listEch):
            it.setEnabled(self.echListWidget.item(ind).checkState())
        QDialog.accept(self)

    def populateListGene(self):
        pix = QPixmap(32, 32)
        for ind, it in enumerate(self.listGene):
            item = QListWidgetItem(it.name)
            item.setCheckState(it.enabled)
            pix.fill(it.color)
            item.setIcon(QIcon(pix))
            self.geneListWidget.addItem(item)

    def populateListEch(self):
        pix = QPixmap(32, 32)
        for ind, it in enumerate(self.listEch):
            item = QListWidgetItem(it.name)
            item.setCheckState(it.enabled)
            pix.fill(it.color)
            item.setIcon(QIcon(pix))
            self.echListWidget.addItem(item)

    def unselectEch(self):
        for item in self.echListWidget.selectedItems():
            self.echListWidget.setItemSelected(item, False)

    def unselectGene(self):
        for item in self.geneListWidget.selectedItems():
            self.geneListWidget.setItemSelected(item, False)

    def up(self):
        activeGenes = self.geneListWidget.selectedItems()
        activeEchs = self.echListWidget.selectedItems()
        if len(activeGenes) != 0:
            row = self.geneListWidget.currentRow()
            if row >= 1:
                thisObj = self.listGene[row]
                prevObj = self.listGene[row-1]
                self.listGene[row] = prevObj
                self.listGene[row-1] = thisObj
                self.geneListWidget.clear()
                self.populateListGene()
                self.geneListWidget.setCurrentRow(row-1)
        elif len(activeEchs) != 0:
            row = self.echListWidget.currentRow()
            if row >= 1:
                thisObj = self.listEch[row]
                prevObj = self.listEch[row-1]
                self.listEch[row] = prevObj
                self.listEch[row-1] = thisObj
                self.echListWidget.clear()
                self.populateListEch()
                self.echListWidget.setCurrentRow(row-1)

    def down(self):
        activeGenes = self.geneListWidget.selectedItems()
        activeEchs = self.echListWidget.selectedItems()
        if len(activeGenes) != 0:
            row  = self.geneListWidget.currentRow()
            if row < self.geneListWidget.count() -1:
                thisObj = self.listGene[row]
                nextObj = self.listGene[row+1]
                self.listGene[row] = nextObj
                self.listGene[row+1] = thisObj
                self.geneListWidget.clear()
                self.populateListGene()
                self.geneListWidget.setCurrentRow(row+1)
        elif len(activeEchs) != 0:
            row  = self.echListWidget.currentRow()
            if row < self.echListWidget.count() -1:
                thisObj = self.listEch[row]
                nextObj = self.listEch[row+1]
                self.listEch[row] = nextObj
                self.listEch[row+1] = thisObj
                self.echListWidget.clear()
                self.populateListEch()
                self.echListWidget.setCurrentRow(row+1)

    def color(self):
        pix = QPixmap(32, 32)
        if self.whois == "geneList":
            item = self.geneListWidget.currentItem()
            row = self.geneListWidget.currentRow()
            col = self.listGene[row].color
            dialog = QColorDialog(self)
            dialog.setCurrentColor(col)
            if dialog.exec_():
                color = dialog.currentColor()
                pix.fill(color)
                item.setIcon(QIcon(pix))
                self.listGene[row].setColor(color)
        elif self.whois == "echList":
            item = self.echListWidget.currentItem()
            row = self.echListWidget.currentRow()
            col = self.listEch[row].color
            dialog = QColorDialog(self)
            dialog.setCurrentColor(col)
            if dialog.exec_():
                color = dialog.currentColor()
                pix.fill(color)
                item.setIcon(QIcon(pix))
                self.listEch[row].setColor(color)


if __name__=="__main__":
    import sys
    from plate import *
    pl = Plaque("../../treated_data.txt")
    for g in pl.listGene:
        setattr(g, 'color', QColor('#000000'))
    pl.listGene[2].setEnabled(Qt.Unchecked)
    app = QApplication(sys.argv)
    f = PropDialog(listGene=pl.listGene[1:])
    f.show()
    app.exec_()
