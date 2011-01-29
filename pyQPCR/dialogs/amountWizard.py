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

import sys
from PyQt4.QtCore import *
from PyQt4.QtGui import *

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"

class AmountWizard(QDialog):
    """
    :attribute editStart: a QLineEdit that contains the starting quantity.
    :attribute editDilution: a QLineEdit that contains the dilution factor.
    :attribute editNumber: a QLineEdit that contains the number of dilutions.
    :attribute amounts: a list that contains the amounts computed
                        with this wizard.
    """

    def __init__(self, parent=None):
        """
        Constructor of AmountWizard

        :param parent: the QWidget parent
        :type parent: PyQt4.QtGui.QWidget
        """
        QDialog.__init__(self, parent)

        lab1 = QLabel("&Starting amount:")
        lab2 = QLabel("&Dilution factor:")
        lab3 = QLabel("&Number of dilutions:")

        self.editStart = QLineEdit("100")
        self.editDilution = QLineEdit("4")
        self.editNumber = QLineEdit("4")

        self.editStart.setValidator(QDoubleValidator(self))
        self.editDilution.setValidator(QDoubleValidator(self))
        self.editNumber.setValidator(QIntValidator(self))

        lab1.setBuddy(self.editStart)
        lab2.setBuddy(self.editDilution)
        lab3.setBuddy(self.editNumber)

        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|
                                     QDialogButtonBox.Cancel)

        layout = QGridLayout()
        layout.addWidget(lab1, 0, 0)
        layout.addWidget(self.editStart, 0, 1)
        layout.addWidget(lab2, 1, 0)
        layout.addWidget(self.editDilution, 1, 1)
        layout.addWidget(lab3, 2, 0)
        layout.addWidget(self.editNumber, 2, 1)
        layout.addWidget(buttonBox, 3, 0, 1, 2)
        self.setLayout(layout)

        self.connect(buttonBox, SIGNAL("accepted()"), self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"), self, SLOT("reject()"))
        self.setWindowTitle("Amount helper")

    def accept(self):
        """
        Overload the accept function. Compute the dilution.
        """
        start, success = self.editStart.text().toDouble()
        dilution, success = self.editDilution.text().toDouble()
        number, success = self.editNumber.text().toInt()

        self.amounts = [start]
        for i in range(number):
            self.amounts.append(start/(dilution*(i+1.)))
        QDialog.accept(self)

if __name__=="__main__":
    app = QApplication(sys.argv)
    f = AmountWizard()
    f.show()
    app.exec_()

