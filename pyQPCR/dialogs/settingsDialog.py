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

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from pyQPCR.utils.odict import *

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Revision$"

class SuffixComboBox(QComboBox):
    """
    This class is an overloading of QComboBox. It define a comboBox
    with a suffix (that is '%' in our case)
    """

    def __init__(self, parent=None):
        """
        Constructor of SuffixComboBox

        :param parent: the QWidget parent
        :type parent: PtQt4.QtCore.QWidget
        """
        QComboBox.__init__(self, parent)

    def addItem(self, item, *args):
        """
        Overload of addItem: add a suffix to the item
        """
        suffix = QString("%")
        item += suffix
        QComboBox.addItem(self, item, *args)

    def addItems(self, items):
        """
        Overload of addItems: add a suffix to the item
        """
        for item in items:
            self.addItem(item)

    def currentText(self):
        """
        Overload of currentText: return the string without
        the suffix.
        """
        cText = QComboBox.currentText(self)
        return cText[:-1]

class SettingsDialog(QDialog):
    """
    This dialog is the configuration dialog of pyQPCR. It is open when 
    the Settings is open. It allows to set the different settings of
    pyQPCR (type of machine and calculations, ...)
    """

    def __init__(self, parent=None, ect=0.3, ctmin=35, confidence=0.9,
                 errtype="normal", machine='Eppendorf', 
                 typeCalc='Relative quantification'):
        QDialog.__init__(self, parent)

        labCalc = QLabel("<b>Calculation :</b>")
        lab0 = QLabel("&Type of calculation :")
        self.typeCalc = QComboBox()
        self.typeCalc.addItems(['Relative quantification', 
                                'Absolute quantification'])
        lab0.setBuddy(self.typeCalc)
        if typeCalc == 'Relative quantification':
            self.typeCalc.setCurrentIndex(0)
        elif typeCalc == 'Absolute quantification':
            self.typeCalc.setCurrentIndex(1)

        labTit = QLabel("<b>Quality Control:</b>")
        lab1 = QLabel("&E(ct) maximum :")
        self.ectLineEdit = QLineEdit("%.2f" % ect)
        self.ectLineEdit.setValidator(QDoubleValidator(self))
        lab1.setBuddy(self.ectLineEdit)
        lab2 = QLabel("&Negative ct maximum :")
        self.ctMinLineEdit = QLineEdit("%.2f" % ctmin)
        self.ctMinLineEdit.setValidator(QDoubleValidator(self))
        lab2.setBuddy(self.ctMinLineEdit)

        labConf = QLabel("<b>Confidence interval :</b>")
        lab3 = QLabel("&Distribution type :")
        self.typeCbox = QComboBox()
        self.types = {}
        self.types[QString('Gaussian')] = 'normal'
        self.types[QString('Student t-test')] = 'student'
        self.typeCbox.addItems(self.types.keys())
        if errtype == "student":
            self.typeCbox.setCurrentIndex(0)
        else:
            self.typeCbox.setCurrentIndex(1)
        lab3.setBuddy(self.typeCbox)

        lab4 = QLabel("&Confidence level :")
        self.confCbox = SuffixComboBox()

        conf = '%.2f' % (100*confidence)
        liste = OrderedDict()
        liste['75.00'] = 0
        liste['80.00'] = 1
        liste['85.00'] = 2
        liste['90.00'] = 3
        liste['95.00'] = 4
        liste['97.50'] = 5
        liste['99.00'] = 6
        liste['99.50'] = 7
        liste['99.75'] = 8
        liste['99.90'] = 9
        liste['99.95'] = 10
        self.confCbox.addItems(liste.keys())
        try:
            self.confCbox.setCurrentIndex(liste[conf])
        except KeyError:
            self.confCbox.setCurrentIndex(4)
        lab4.setBuddy(self.confCbox)

        labMachine = QLabel("<b>PCR device :</b>")
        lab5 = QLabel("&Machine type : ")
        self.machBox = QComboBox()
        self.machBox.addItems(['Eppendorf', 'Applied StepOne', 'Applied 7000', 
                               'Applied 7500', 'Applied 7700', 'Applied 7900',
                               'Biorad MyIQ', 'Cepheid SmartCycler', 
                               'Qiagen Corbett', 'Roche LightCycler 480'])
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
        elif machine == 'Biorad MyIQ':
            self.machBox.setCurrentIndex(6)
        elif machine == 'Cepheid SmartCycler':
            self.machBox.setCurrentIndex(7)
        elif machine == 'Qiagen Corbett':
            self.machBox.setCurrentIndex(8)
        elif machine == 'Roche LightCycler 480':
            self.machBox.setCurrentIndex(9)
        lab5.setBuddy(self.machBox)

        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|
                                     QDialogButtonBox.Cancel)

        gLayout = QGridLayout()
        gLayout.addWidget(labCalc, 0, 0, 1, 2)
        gLayout.addWidget(lab0, 1, 0)
        gLayout.addWidget(self.typeCalc, 1, 1)
        gLayout.addWidget(labTit, 2, 0, 1, 2)
        gLayout.addWidget(lab1, 3, 0)
        gLayout.addWidget(self.ectLineEdit, 3, 1)
        gLayout.addWidget(lab2, 4, 0)
        gLayout.addWidget(self.ctMinLineEdit, 4, 1)
        gLayout.addWidget(labConf, 5, 0, 1, 2)
        gLayout.addWidget(lab3, 6, 0)
        gLayout.addWidget(self.typeCbox, 6, 1)
        gLayout.addWidget(lab4, 7, 0)
        gLayout.addWidget(self.confCbox, 7, 1)
        gLayout.addWidget(labMachine, 8, 0, 1, 2)
        gLayout.addWidget(lab5, 9, 0)
        gLayout.addWidget(self.machBox, 9, 1)
        gLayout.addWidget(buttonBox, 10, 0, 1, 2)

        self.setLayout(gLayout)

        self.setWindowTitle("%s Settings" % QApplication.applicationName())
        # Connections
        self.connect(buttonBox, SIGNAL("accepted()"), self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"), self, SLOT("reject()"))


if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    form = SettingsDialog()
    form.show()
    app.exec_()
