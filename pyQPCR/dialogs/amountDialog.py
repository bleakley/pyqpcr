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

class AmountDialog(QDialog):
    
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
        self.setWindowTitle("New amount")

    def reformatFloat(self, input):
        output = QString("%.2f" % float(input))
        return output

    def populateList(self):
        self.listWidget.clear()
        for it in self.project.hashAmount.keys():
            if it != '':
                item = QListWidgetItem(it)
                item.setStatusTip(it)
                self.listWidget.addItem(item)

    def add(self):
        dialog = AddAmDialog(self)
        if dialog.exec_():
            am = dialog.am.text()
            amname = self.reformatFloat(am)
            if not self.project.hashAmount.has_key(amname):
                self.project.hashAmount[amname] = float(am)
                self.populateList()
            else:
                QMessageBox.warning(self, "Already exist",
                            "The amount %s is already defined !" % amname)

    def edit(self):
        if len(self.listWidget.selectedItems()) == 0:
            return
        am_before = self.listWidget.currentItem().statusTip()

        dialog = AddAmDialog(self, am=am_before)
        if dialog.exec_():
            am = dialog.am.text()
            amname = self.reformatFloat(am)
            if self.project.dicoAmount.has_key(am_before) and amname != am_before:
                ind = self.project.dicoAmount.index(am_before)
                self.project.dicoAmount.insert(ind, amname, 
                          self.project.dicoAmount[am_before])
                ind = self.project.hashAmount.index(am_before)
                self.project.hashAmount.insert(ind, amname, 
                          self.project.hashAmount[am_before])
                for well in self.project.dicoAmount[amname]:
                    well.setAmount(float(am))
                self.project.dicoAmount.__delitem__(am_before)
                self.project.hashAmount.__delitem__(am_before)
                self.project.unsaved = True
            self.populateList()

    def remove(self):
        ams = []
        if len(self.listWidget.selectedItems()) == 0:
            return
        for it in self.listWidget.selectedItems():
            am = it.statusTip()
            ams.append(am)

        reply = QMessageBox.question(self, "Remove",
                        "Remove %s ?" % ams,
                        QMessageBox.Yes|QMessageBox.No)
        if reply == QMessageBox.Yes:
            for am in ams:
                if self.project.dicoAmount.has_key(am):
                    for well in self.project.dicoAmount[am]:
                        well.setAmount('')
                    self.project.setDicoAm()
                    self.project.hashAmount.__delitem__(am)
            self.project.unsaved = True
            self.populateList()


class AddAmDialog(QDialog):
    
    def __init__(self, parent=None, am=None):
        self.parent = parent
        QDialog.__init__(self, parent)
        lab = QLabel("Amount:")
        if am is not None:
            self.am = QLineEdit(am)
        else:
            self.am = QLineEdit()
        self.am.setValidator(QDoubleValidator(self))
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|
                                     QDialogButtonBox.Cancel)

        layout = QGridLayout()
        layout.addWidget(lab, 0, 0)
        layout.addWidget(self.am, 0, 1)
        layout.addWidget(buttonBox, 1, 0, 1, 2)
        self.setLayout(layout)

        self.connect(buttonBox, SIGNAL("accepted()"), self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"), self, SLOT("reject()"))
        self.setWindowTitle("New sample")

if __name__=="__main__":
    app = QApplication(sys.argv)
    pl = project('toto.csv')
    f = AmountDialog(project=pl)
    f.show()
    app.exec_()
