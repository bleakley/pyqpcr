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
from amountWizard import *
import pyQPCR.qrc_resources

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"

class AmountDialog(QDialog):
    """
    This object is a dialog that allows to Add/Edit/Remove the different amounts
    of your plates.

    :attribute listWidget: the main displayed list of items that contains the 
                           different amounts.
    :attribute project: a copy (copy.deepcopy) of the project
    """
    
    def __init__(self, parent=None, project=None):
        """
        Constructor of AmountDialog

        :param parent: the QWidget parent
        :type parent: PyQt4.QtGui.QWidget
        :param project: the current project
        :type project: pyQPCR.project.Project
        """
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
        buttonWizard = QPushButton("&Automatic dilutions...")
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|
                                     QDialogButtonBox.Cancel)

        vlayout = QVBoxLayout()
        vlayout.addWidget(buttonAdd)
        vlayout.addWidget(buttonEdit)
        vlayout.addWidget(buttonRemove)
        vlayout.addWidget(buttonWizard)
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
        self.connect(buttonWizard, SIGNAL("clicked()"), self.automatic)
        self.setWindowTitle("New amount")

    def reformatFloat(self, input):
        """
        A method to reformat a float into a string with the correct format.

        :param input: the input string (representing a float or an integer)
        :type input: PyQt4.QtCore.QString
        """
        if float(input) >= 1e-2 and float(input) <= 1e3:
            output = QString("%.2f" % float(input))
        else:
            output = QString("%.2e" % float(input))
        return output

    def populateList(self):
        """
        A method to fill up the QListWidget that contains the different amounts.
        """
        self.listWidget.clear()
        for it in self.project.hashAmount.keys():
            if it != '':
                item = QListWidgetItem(it)
                item.setStatusTip(it)
                self.listWidget.addItem(item)

    def add(self):
        """
        This method calls the AddAmDialog wizard that allows to add a new amount.
        A warning is displayed if the amount already exists.
        """
        dialog = AddAmDialog(self)
        if dialog.exec_():
            am = dialog.am.text()
            amname = self.reformatFloat(am)
            if not self.project.hashAmount.has_key(amname):
                self.project.hashAmount[amname] = float(am)
                self.project.dicoAmount[amname] = []
                self.populateList()
            else:
                QMessageBox.warning(self, "Already exist",
                            "The amount <b>%s</b> has been already defined !" % amname)


    def automatic(self):
        """
        This method calls the AmountWizard dialog that helps to generate quantities
        automatically.
        """
        dialog = AmountWizard(self)
        if dialog.exec_():
            alreadyDefined = []
            for am in dialog.amounts:
                amname = self.reformatFloat(am)
                if not self.project.hashAmount.has_key(amname):
                    self.project.hashAmount[amname] = am
                    self.project.dicoAmount[amname] = []
                    self.populateList()
                else:
                    alreadyDefined.append(str(amname))

            if len(alreadyDefined) > 0:
                QMessageBox.warning(self, "Already exist",
                    "The amount <b>%s</b> have been already defined !" % alreadyDefined)

    def edit(self):
        """
        This method call the AddAmDialog wizard that allows to edit the existing amounts.
        """
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
                self.project.hashAmount[amname] = float(am)
                self.project.unsaved = True
            self.populateList()

    def remove(self):
        """
        This is used to remove an existing amount.
        """
        ams = [] ; sts = []
        if len(self.listWidget.selectedItems()) == 0:
            return
        for it in self.listWidget.selectedItems():
            am = it.statusTip()
            if self.project.hashAmount[am] >= 1e-2 and self.project.hashAmount[am] <= 1e3:
                st = "%.2f" % self.project.hashAmount[am]
            else:
                st = "%.2e" % self.project.hashAmount[am]
            sts.append(st)
            ams.append(am)

        reply = QMessageBox.question(self, "Remove",
                "Remove the following amounts: <b>%s</b> ?" % sts,
                        QMessageBox.Yes|QMessageBox.No)
        if reply == QMessageBox.Yes:
            for am in ams:
                delete = False
                if self.project.hashAmount.has_key(am):
                    self.project.hashAmount.__delitem__(am)
                    delete = True
            if delete:
                for am in ams:
                    if self.project.dicoAmount.has_key(am):
                        for well in self.project.dicoAmount[am]:
                            well.setAmount('')
                self.project.setDicoAm()

            self.project.unsaved = True
            self.populateList()


class AddAmDialog(QDialog):
    """
    This object is a small dialog used to Add/Edit a new amount in 
    your experiment.

    :attribute am: the amount edited
    """
    
    def __init__(self, parent=None, am=None):
        """
        Constructor of AddAmDialog

        :param parent: the QWidget parent
        :type parent: PyQt4.QtGui.QWidget
        :param am: the current amount we want to Edit (or Add)
        :type project: pyQPCR.project.Project
        """
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
