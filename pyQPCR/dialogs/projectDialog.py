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

from PyQt4.QtGui import *
from PyQt4.QtCore import *
#from pyQPCR.utils.odict import *
from odict import *

__author__ = "$Author: tgastine $"
__date__ = "$Date: 2010-01-24 11:25:34 +0100 (dim. 24 janv. 2010) $"
__version__ = "$Rev: 128 $"

class NewProjectDialog(QDialog):
    
    def __init__(self, parent=None, pwd=None):
        self.parent = parent
        self.pwd = pwd
        self.fileNames = OrderedDict()
        QDialog.__init__(self, parent)

        lab1 = QLabel("1. &Project name")
        self.edt = QLineEdit()
        lab1.setBuddy(self.edt)

        lab2 = QLabel("2. &Machine type")
        cbox = QComboBox()
        lab2.setBuddy(cbox)
        cbox.addItem("Eppendorf")

        lab3 = QLabel("3. Plates files")
        self.listFiles = QListWidget()
        self.listFiles.setAlternatingRowColors(True)
        self.listFiles.setSelectionMode(3)

        btnAdd = QPushButton("&Add")
        btnRemove = QPushButton("&Remove")
        vLay = QVBoxLayout()
        vLay.addWidget(btnAdd)
        vLay.addWidget(btnRemove)
        vLay.addStretch()
        hLay = QHBoxLayout()
        hLay.addWidget(self.listFiles)
        hLay.addLayout(vLay)

        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|
                                     QDialogButtonBox.Cancel)

        finalLayout = QVBoxLayout()
        finalLayout.addWidget(lab1)
        finalLayout.addWidget(self.edt)
        finalLayout.addWidget(lab2)
        finalLayout.addWidget(cbox)
        finalLayout.addWidget(lab3)
        finalLayout.addLayout(hLay)
        finalLayout.addWidget(buttonBox)

        self.setLayout(finalLayout)

        self.connect(buttonBox, SIGNAL("accepted()"), self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"), self, SLOT("reject()"))
        self.connect(btnAdd, SIGNAL("clicked()"), self.addPlate)
        self.connect(btnRemove, SIGNAL("clicked()"), self.removePlate)
        self.setWindowTitle("New project")

    def populateList(self):
        self.listFiles.clear()
        for fname in self.fileNames.keys():
            item = QListWidgetItem(fname)
            self.listFiles.addItem(item)

    def addPlate(self):
        dir = self.pwd if self.pwd is not None else "."
        formats =[u"*.txt", u"*.csv"]
        fileNames = QFileDialog.getOpenFileNames(self,
                       "pyQPCR - Choose plates", dir,
                       "Input files (%s)" % " ".join(formats))
        if fileNames:
            for fname in fileNames:
                self.fileNames[QFileInfo(fname).fileName()] = fname
            self.populateList()

    def removePlate(self):
        if len(self.listFiles.selectedItems()) == 0:
            return
        for it in self.listFiles.selectedItems():
            self.fileNames.pop(it.text())
        self.populateList()

    def accept(self):
        if self.edt.text() == '':
            QMessageBox.warning(self, "No project name",
               "<b>Warning</b>: you must give a project name ! " )
        else:
            if self.edt.text().endsWith('xml') or self.edt.text().endsWith('XML'):
                self.projectName = self.edt.text()
            else:
                self.projectName = QString("%s.xml" % self.edt.text())
            QDialog.accept(self)


if __name__=="__main__":
    import sys
    app = QApplication(sys.argv)
    f = NewProjectDialog()
    f.show()
    app.exec_()
