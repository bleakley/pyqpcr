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
from pyQPCR.dialogs.customWidgets import *
from PyQt4.QtGui import *
from PyQt4.QtCore import *

class PropDialog(QDialog):
    
    def __init__(self, parent=None, listObj=None):
        self.parent = parent
        QDialog.__init__(self, parent)

        self.listWidget = QListWidget()
        self.listWidget.setAlternatingRowColors(True)
        pix = QPixmap(32, 32)
        if listObj is not None:
            self.listObj = copy.deepcopy(listObj)
            self.populateList()

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
        hlayout.addWidget(self.listWidget)
        hlayout.addLayout(vlayout)
        layout = QVBoxLayout()
        layout.addLayout(hlayout)
        layout.addWidget(buttonBox)
        self.setLayout(layout)

        self.connect(buttonBox, SIGNAL("accepted()"), self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"), self, SLOT("reject()"))
        self.connect(buttonUp, SIGNAL("clicked()"), self.up)
        self.connect(buttonDown, SIGNAL("clicked()"), self.down)
        self.connect(buttonColor, SIGNAL("clicked()"), self.color)
        self.setWindowTitle("Plot properties")

    def accept(self):
        for ind,it in enumerate(self.listObj):
            it.setEnabled(self.listWidget.item(ind).checkState())
        QDialog.accept(self)

    def populateList(self):
        pix = QPixmap(32, 32)
        for ind, it in enumerate(self.listObj):
            item = QListWidgetItem(it.name)
            item.setCheckState(it.enabled)
            pix.fill(it.color)
            item.setIcon(QIcon(pix))
            self.listWidget.addItem(item)

    def up(self):
        row = self.listWidget.currentRow()
        if row >= 1:
            thisObj = self.listObj[row]
            prevObj = self.listObj[row-1]
            self.listObj[row] = prevObj
            self.listObj[row-1] = thisObj
            self.listWidget.clear()
            self.populateList()
            self.listWidget.setCurrentRow(row-1)

    def down(self):
        row  = self.listWidget.currentRow()
        if row < self.listWidget.count() -1:
            thisObj = self.listObj[row]
            nextObj = self.listObj[row+1]
            self.listObj[row] = nextObj
            self.listObj[row+1] = thisObj
            self.listWidget.clear()
            self.populateList()
            self.listWidget.setCurrentRow(row+1)

    def color(self):
        pix = QPixmap(32, 32)
        item = self.listWidget.currentItem()
        row = self.listWidget.currentRow()
        col = self.listObj[row].color
        color = QColorDialog.getColor(col, self)
        pix.fill(color)
        item.setIcon(QIcon(pix))
        self.listObj[row].setColor(color)


if __name__=="__main__":
    import sys
    from plaque import *
    pl = Plaque("sortiesrealplex/ref-machine.txt")
    for g in pl.listGene:
        setattr(g, 'color', QColor('#000000'))
    pl.listGene[2].setEnabled(Qt.Unchecked)
    app = QApplication(sys.argv)
    f = PropDialog(listObj=pl.listGene[1:])
    f.show()
    app.exec_()
