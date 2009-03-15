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

__author__ = "$Author:$"
__date__ = "$Date:$"
__version__ = "$Revision:$"


class SettingsDialog(QDialog):

    def __init__(self, parent=None):
        self.parent = parent
        QDialog.__init__(self, parent)

        value = 0.3
        ctmin = 35
        labTit = QLabel("<b>Quality Control:</b>")
        lab1 = QLabel("E(ct) >")
        edit1 = QLineEdit("%.2f" % value)
        lab1.setBuddy(edit1)
        lab2 = QLabel("Negative ct >")
        edit2 = QLineEdit("%.2f" % ctmin)
        lab2.setBuddy(edit2)
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|
                                     QDialogButtonBox.Cancel)

        gLayout = QGridLayout()
        gLayout.addWidget(lab1, 0, 0)
        gLayout.addWidget(edit1, 0, 1)
        gLayout.addWidget(lab2, 1, 0)
        gLayout.addWidget(edit2, 1, 1)

        layout = QVBoxLayout()
        layout.addWidget(labTit)
        layout.addLayout(gLayout)
        layout.addWidget(buttonBox)
        self.setLayout(layout)

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
