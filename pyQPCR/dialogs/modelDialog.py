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
from pyQPCR.project import *
from pyQPCR.qrc_resources import *
import os
import copy

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"

class ModelDialog(QDialog):
    
    def __init__(self, parent=None, pr=None, pwd=None):
        """
        Constructor of ModelDialog

        :param parent: the QWidget parent
        :type parent: PyQt4.QtGui.QWidget
        """
        QDialog.__init__(self, parent)
        self.pwd = pwd
        self.project = copy.deepcopy(pr)
        self.setWindowTitle("Assistant model")

        lab0 = QLabel("<b>1. Choose a model:</b>")

        self.inPlate = QRadioButton("Plate from the &current project:")
        self.inPlate.setChecked(Qt.Checked)
        self.cboxInPlate = QComboBox()

        self.outPlate = QRadioButton("Plate from &another project:")
        self.outWidget = QWidget()
        self.file = QLineEdit()
        self.file.setReadOnly(True)
        btn = QToolButton()
        ic = QIcon(":/fileopen.png")
        btn.setIcon(ic)
        hLayout = QHBoxLayout()
        hLayout.addWidget(self.file)
        hLayout.addWidget(btn)

        self.cboxOutPlate = QComboBox()

        vLay = QVBoxLayout()
        vLay.addLayout(hLayout)
        vLay.addWidget(self.cboxOutPlate)
        self.outWidget.setLayout(vLay)
        self.outWidget.setVisible(False)

        lab1 = QLabel("<b>2. Choose a target:</b>")
        self.plateCalc = QComboBox()
        lab1.setBuddy(self.plateCalc)


        lab2 = QLabel("<b>3. Summary:</b>")

        self.info = QLabel("The plate ... will be modfied")

        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|
                                     QDialogButtonBox.Cancel)

        vLayout = QVBoxLayout()
        vLayout.addWidget(self.inPlate)
        vLayout.addWidget(self.cboxInPlate)
        vLayout.addWidget(self.outPlate)
        vLayout.addWidget(self.outWidget)
        vLayout.addStretch(1)

        finalLayout = QVBoxLayout()
        finalLayout.addWidget(lab0)
        finalLayout.addLayout(vLayout)
        finalLayout.addWidget(lab1)
        finalLayout.addWidget(self.plateCalc)
        finalLayout.addWidget(lab2)
        finalLayout.addWidget(self.info)
        finalLayout.addWidget(buttonBox)
        finalLayout.addStretch(1)
        finalLayout.setSizeConstraint(QLayout.SetFixedSize)

        self.setLayout(finalLayout)

        self.populateCbox()

        self.connect(buttonBox, SIGNAL("accepted()"), self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"), self, SLOT("reject()"))
        self.connect(self.inPlate, SIGNAL("clicked()"), self.inPlateClicked)
        self.connect(self.outPlate, SIGNAL("clicked()"), self.outPlateClicked)
        self.connect(btn, SIGNAL("clicked()"), self.setFilePath)

    def outPlateClicked(self):
        self.outWidget.setVisible(True)
        self.cboxInPlate.setVisible(False)

    def inPlateClicked(self):
        self.cboxInPlate.setVisible(True)
        self.outWidget.setVisible(False)

    def populateCbox(self):
        for plate in self.project.dicoPlates.keys():
            self.cboxInPlate.addItem(plate)
            self.plateCalc.addItem(plate)

        if len(self.project.dicoPlates.keys()) >= 2:
            self.cboxInPlate.setCurrentIndex(0)
            self.plateCalc.setCurrentIndex(1)

    def populateCboxOut(self):
        """
        This method populates the combo box if an external
        project has been selected.
        """
        self.cboxOutPlate.clear()
        for plate in self.projectOut.dicoPlates.keys():
            self.cboxOutPlate.addItem(plate)

        if len(self.projectOut.dicoPlates.keys()) >= 2:
            self.cboxOutPlate.setCurrentIndex(0)

    def setFilePath(self):
        """
        A method to select an external project.
        """
        dir = self.pwd if self.pwd is not None else "."
        formats = [u"*.xml"]
        fileName = QFileDialog.getOpenFileName(self, 'Choose a pyQPCR project', dir,
                                               "Input files (%s)" % (" ".join(formats)))
        if fileName:
            self.file.setText(fileName)
            try:
                self.projectOut = Project(fileName)
                self.populateCboxOut()
            except ProjectError,e:
                QMessageBox.warning(self, "Problem in import !", "%s" % str(e))

    def accept(self):
        """
        Overload of the 'accept' method.
        """
        if self.inPlate.isChecked() is True:
            model = self.cboxInPlate.currentText()
            target = self.plateCalc.currentText()

            if self.project.dicoPlates[target].type == self.project.dicoPlates[model].type:
                brokenWell = []
                for well in self.project.dicoPlates[target].listePuits:
                    try:
                        wellmodel = getattr(self.project.dicoPlates[model], str(well.name))
                        well.setType(wellmodel.type)
                        well.setGene(wellmodel.gene)
                        well.setEch(wellmodel.ech)
                        well.setAmount(wellmodel.amount)
                    except AttributeError:
                        brokenWell.append(well.name)

                self.project.dicoPlates[target].geneRef = self.project.dicoPlates[model].geneRef
                self.project.dicoPlates[target].echRef = self.project.dicoPlates[model].echRef
                self.project.dicoPlates[target].setDicoGene()
                self.project.dicoPlates[target].setDicoEch()
                self.project.setDicoAm()

                if len(brokenWell) != 0:
                    QMessageBox.information(self, "Some wells have not been changed",
                            "The following wells: <b>%s</b> do not exist in the model." % brokenWell)
            else:
                QMessageBox.warning(self, "Incompatible model and target",
                   "<b>Warning</b>: model and target must have the same number of wells ! " )
                return

        else: # an external plate has been chosen
            model = self.cboxOutPlate.currentText()
            target = self.plateCalc.currentText()
            if model == '':
                QMessageBox.warning(self, "No file selected",
               "<b>Warning</b>: you have to select a valid project ! " )
                return

            if self.project.dicoPlates[target].type == self.projectOut.dicoPlates[model].type:
                brokenWell = []
                for well in self.project.dicoPlates[target].listePuits:
                    try:
                        wellmodel = getattr(self.projectOut.dicoPlates[model], str(well.name))
                        well.setType(wellmodel.type)
                        well.setGene(wellmodel.gene)
                        well.setEch(wellmodel.ech)
                        well.setAmount(wellmodel.amount)
                    except AttributeError:
                        brokenWell.append(well.name)

                self.project.dicoPlates[target].geneRef = self.projectOut.dicoPlates[model].geneRef
                self.project.dicoPlates[target].echRef = self.projectOut.dicoPlates[model].echRef
                self.project.dicoPlates[target].setDicoGene()
                self.project.dicoPlates[target].setDicoEch()
                self.project.setDicoAm()

                self.project.initLocGene(plate=self.project.dicoPlates[target])
                self.project.initLocEch(plate=self.project.dicoPlates[target])
                self.project.initLocAm(plate=self.project.dicoPlates[target])

                if len(brokenWell) != 0:
                    QMessageBox.information(self, "Some wells have not been changed",
                            "The following wells: <b>%s</b> do not exist in the model." % brokenWell)
            else:
                QMessageBox.warning(self, "Incompatible model and target",
                   "<b>Warning</b>: model and target must have the same number of wells ! " )
                return
        QDialog.accept(self)

if __name__=="__main__":
    import sys
    from project import *
    pr = Project("2plates_color.xml")
    app = QApplication(sys.argv)
    f = ModelDialog(pr=pr)
    f.show()
    app.exec_()
