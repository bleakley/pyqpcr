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
import pyQPCR.qrc_resources
from pyQPCR.utils.odict import *
import os

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"

class NewProjectDialog(QDialog):
    
    def __init__(self, parent=None, pwd=None):
        self.parent = parent
        self.pwd = pwd
        self.fileNames = OrderedDict()
        QDialog.__init__(self, parent)

        lab1 = QLabel("<b>1. &Project name</b>")
        self.edt = QLineEdit()
        lab1.setBuddy(self.edt)

        lab2 = QLabel("<b>2. &Machine type</b>")
        self.machBox = QComboBox()
        lab2.setBuddy(self.machBox)
        self.machBox.addItem("Eppendorf")
        self.machBox.addItem("Applied StepOne")
        self.machBox.addItem("Applied 7000")
        self.machBox.addItem("Applied 7500")
        self.machBox.addItem("Roche LightCycler 480")

        lab3 = QLabel("<b>3. Plates files</b>")
        self.listFiles = QListWidget()
        self.listFiles.setAlternatingRowColors(True)
        self.listFiles.setSelectionMode(QAbstractItemView.ExtendedSelection)

        lab4 = QLabel("<b>4. Destination directory</b>")
        self.file = QLineEdit()
        lab1.setBuddy(self.file)
        self.file.setReadOnly(True)
        btn = QToolButton()
        ic = QIcon(":/fileopen")
        btn.setIcon(ic)
        hLay2 = QHBoxLayout()
        hLay2.addWidget(self.file)
        hLay2.addWidget(btn)

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
        finalLayout.addWidget(self.machBox)
        finalLayout.addWidget(lab3)
        finalLayout.addLayout(hLay)
        finalLayout.addWidget(lab4)
        finalLayout.addLayout(hLay2)
        finalLayout.addWidget(buttonBox)

        self.setLayout(finalLayout)

        self.connect(buttonBox, SIGNAL("accepted()"), self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"), self, SLOT("reject()"))
        self.connect(btnAdd, SIGNAL("clicked()"), self.addPlate)
        self.connect(btnRemove, SIGNAL("clicked()"), self.removePlate)
        self.connect(btn, SIGNAL("clicked()"), self.setFilePath)

        self.setWindowTitle("New project")

    def populateList(self):
        self.listFiles.clear()
        for fname in self.fileNames.keys():
            item = QListWidgetItem(fname)
            self.listFiles.addItem(item)

    def addPlate(self):
        dir = self.pwd if self.pwd is not None else "."
        if self.machBox.currentText() == 'Eppendorf':
            formats =[u"*.txt", u"*.csv"]
            type = 'Eppendorf machines'
        elif self.machBox.currentText() == 'Applied StepOne':
            formats =[u"*.txt", u"*.csv"]
            type = 'Applied StepOne machines'
        elif self.machBox.currentText() == 'Applied 7000':
            formats =[u"*.csv"]
            type = 'Applied 7000 machines'
        elif self.machBox.currentText() == 'Applied 7500':
            formats =[u"*.txt", u"*.csv"]
            type = 'Applied 7500 machines'
        elif self.machBox.currentText() == 'Roche LightCycler 480':
            formats =[u"*.txt"]
            type = 'Roche LightCycler 480'
        fileNames = QFileDialog.getOpenFileNames(self,
                       "pyQPCR - Choose plates", dir,
                       "Input files [%s] (%s)" % (type, " ".join(formats)))
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
        elif self.file.text() == '':
            QMessageBox.warning(self, "No project directory",
               "<b>Warning</b>: you must choose a directory for the project ! " )
        else:
            if self.edt.text().endsWith('xml') or self.edt.text().endsWith('XML'):
                projectName = self.edt.text()
            else:
                projectName = QString("%s.xml" % self.edt.text())
# Gestion du / ou du \ selon l'OS utilise avec os.sep
            self.projectName = "%s%s%s" % (self.workDir, os.sep, projectName)
            self.machineType = self.machBox.currentText()
            if os.path.exists(self.projectName):
                QMessageBox.warning(self, "This project already exists",
                  """<b>Warning</b>: you must choose a project name that
                     does not exists. %s is already used""" % projectName)
            else:
                QDialog.accept(self)

    def setFilePath(self):
        dir = QFileDialog.getExistingDirectory(self, 'Choose the directory')
        if dir:
            self.workDir = dir
            self.file.setText(self.workDir)


if __name__=="__main__":
    import sys
    app = QApplication(sys.argv)
    f = NewProjectDialog()
    f.show()
    app.exec_()
