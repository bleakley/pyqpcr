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
    """
    This object is used to display a wizard that helps to create
    a new project. It allows to enter the name of the project, the
    type of calculation, the type of PCR device and the plates 
    you want to import.

    :attribute pwd: the path of the project
    :attribute fileNames: a list that contains the name of the plates
                          that will be imported
    :attribute typeCalc: the type of calculation (absolute or relative
                         quantification
    :attribute edt: the name of the project
    :attribute machbox: the type of PCR device used in the experiment
    :attribute listFiles: a list that contains the name of the plates
    :attribute file: the path where the project is going to be saved
    :attribute workdir: the directory where the project will be saved
    """
    
    def __init__(self, parent=None, pwd=None, machine='Eppendorf'):
        """
        Constructor of NewProjectDialog

        :param parent: the QWidget parent
        :type parent: PyQt4.QtGui.QWidget
        :param pwd: the path of the new project
        :type pwd: PyQt4.QtCore.QString
        :param machine: the type of PCR device used in the experiment
        :type machine: PyQt4.QtCore.QString
        """
        self.pwd = pwd
        self.fileNames = OrderedDict()
        QDialog.__init__(self, parent)

        lab0 = QLabel("<b>1. &Calculation type</b>")
        self.typeCalc = QComboBox()
        self.typeCalc.addItems(['Relative quantification', 'Absolute quantification'])
        lab0.setBuddy(self.typeCalc)


        lab1 = QLabel("<b>2. &Project name</b>")
        self.edt = QLineEdit()
        lab1.setBuddy(self.edt)

        lab2 = QLabel("<b>3. &Machine type</b>")
        self.machBox = QComboBox()
        lab2.setBuddy(self.machBox)
        self.machBox.addItems(['Eppendorf', 'Applied StepOne', 'Applied 7000',
                               'Applied 7500', 'Applied 7700', 'Applied 7900',
                               'Biorad C1000', 'Biorad MyIQ', 'Cepheid SmartCycler', 
                               'Qiagen Corbett', 'Roche LightCycler 480',
                               'Stratagene Mx3000'])
        if machine == 'Eppendorf':
            self.machBox.setCurrentIndex(0)
        elif machine == 'Applied StepOne':
            self.machBox.setCurrentIndex(1)
        elif machine == 'Applied 7000':
            self.machBox.setCurrentIndex(2)
        elif machine == 'Applied 7500':
            self.machBox.setCurrentIndex(3)
        elif machine == 'Applied 7700':
            self.machBox.setCurrentIndex(4)
        elif machine == 'Applied 7900':
            self.machBox.setCurrentIndex(5)
        elif machine == 'Biorad C1000':
            self.machBox.setCurrentIndex(6)
        elif machine == 'Biorad MyIQ':
            self.machBox.setCurrentIndex(7)
        elif machine == 'Cepheid SmartCycler':
            self.machBox.setCurrentIndex(8)
        elif machine == 'Qiagen Corbett':
            self.machBox.setCurrentIndex(9)
        elif machine == 'Roche LightCycler 480':
            self.machBox.setCurrentIndex(10)
        elif machine == 'Stratagene Mx3000':
            self.machBox.setCurrentIndex(11)

        lab3 = QLabel("<b>4. Plates files</b>")
        self.listFiles = QListWidget()
        self.listFiles.setAlternatingRowColors(True)
        self.listFiles.setSelectionMode(QAbstractItemView.ExtendedSelection)

        lab4 = QLabel("<b>5. Destination directory</b>")
        self.file = QLineEdit()
        lab1.setBuddy(self.file)
        self.file.setReadOnly(True)

        btn = QToolButton()
        ic = QIcon(":/fileopen.png")
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
        finalLayout.addWidget(lab0)
        finalLayout.addWidget(self.typeCalc)
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
        """
        A method to populate the list of files displayed in the QDialog.
        """
        self.listFiles.clear()
        for fname in self.fileNames.keys():
            item = QListWidgetItem(fname)
            self.listFiles.addItem(item)

    def addPlate(self):
        """
        A method to add a plate to the project.
        """
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
        elif self.machBox.currentText() == 'Applied 7700':
            formats =[u"*.csv"]
            type = 'Applied 7700 machines'
        elif self.machBox.currentText() == 'Applied 7900':
            formats =[u"*.txt"]
            type = 'Applied 7900 machines'
        elif self.machBox.currentText() == 'Biorad C1000':
            formats =[u"*.txt"]
            type = 'Biorad C1000 machines'
        elif self.machBox.currentText() == 'Biorad MyIQ':
            formats =[u"*.csv"]
            type = 'Biorad MyIQ machines'
        elif self.machBox.currentText() == 'Cepheid SmartCycler':
            formats =[u"*.csv"]
            type = 'Cepheid SmartCycler machines'
        elif self.machBox.currentText() == 'Qiagen Corbett':
            formats =[u"*.csv"]
            type = 'Qiagen Corbett machines'
        elif self.machBox.currentText() == 'Roche LightCycler 480':
            formats =[u"*.txt"]
            type = 'Roche LightCycler 480'
        elif self.machBox.currentText() == 'Stratagene Mx3000':
            formats =[u"*.txt"]
            type = 'Stratagene Mx3000'
        fileNames = QFileDialog.getOpenFileNames(self,
                       "pyQPCR - Choose plates", dir,
                       "Input files [%s] (%s)" % (type, " ".join(formats)))
        if fileNames:
            for fname in fileNames:
                self.fileNames[QFileInfo(fname).fileName()] = fname
            self.populateList()

    def removePlate(self):
        """
        A method to remove a plate from the project.
        """
        if len(self.listFiles.selectedItems()) == 0:
            return
        for it in self.listFiles.selectedItems():
            self.fileNames.pop(it.text())
        self.populateList()

    def accept(self):
        """
        Overload of the 'accept' method.
        """
        if self.edt.text() == '':
            QMessageBox.warning(self, "No project name",
               "<b>Warning</b>: you must give a project name ! " )
        elif self.file.text() == '':
            QMessageBox.warning(self, "No project directory",
               "<b>Warning</b>: you must choose a directory to save the project! " )
        else:
            if self.edt.text().endsWith('xml') or self.edt.text().endsWith('XML'):
                projectName = self.edt.text()
            else:
                projectName = QString("%s.xml" % self.edt.text())
            # Gestion du / ou du \ selon l'OS utilise avec os.sep
            self.projectName = "%s%s%s" % (self.workDir, os.sep, projectName)
            self.machineType = self.machBox.currentText()
            self.calculationType = self.typeCalc.currentText()
            if os.path.exists(self.projectName):
                QMessageBox.warning(self, "This project already exists",
                  """<b>Warning</b>: you must choose a project name that
                     does not exists. %s is already used""" % projectName)
            else:
                QDialog.accept(self)

    def setFilePath(self):
        """
        A method to set where the project is going to be saved.
        """
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
