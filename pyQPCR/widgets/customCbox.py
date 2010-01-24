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
#from pyQPCR.wellGeneSample import *
import pyQPCR.qrc_resources

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"


class HeaderModel(QStandardItemModel):

    def __init__(self, parent=None):
        QStandardItemModel.__init__(self, parent)

    def flags(self, index):
        flags = QStandardItemModel.flags(self, index)
        if not index.isValid():
            return flags
        value = index.data(Qt.UserRole).toString()
        if value == "header":
            flags &= Qt.ItemIsSelectable
            flags &= Qt.ItemIsEnabled
        return flags

class ComboDelegate(QItemDelegate):

    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)

    def paint(self, painter, option, index):
        text = index.model().data(index, Qt.DisplayRole).toString()
        if not (index.model().flags(index) and Qt.ItemIsEnabled):
            fontBold = QFont()
            fontBold.setBold(True)
            fontBold.setItalic(True)
            painter.setFont(fontBold)
            palette = QApplication.palette()
            painter.drawText(option.rect, Qt.AlignLeft|Qt.TextSingleLine,
                             text)
        else:
            QItemDelegate.paint(self, painter, option, index)

class GeneEchComboBox(QComboBox):

    def __init__(self, parent=None):
        QComboBox.__init__(self, parent)
        self.setModel(HeaderModel(self))
        self.setItemDelegate(ComboDelegate(self))

    def addItem(self, obj, *args):
        if hasattr(obj, "name"):
            item = obj.name
        else:
            try:
                item =  QString("%.2f" % obj)
            except TypeError:
                item = QString(obj)
        QComboBox.addItem(self, item, *args)

    def addItems(self, hashObj, editDialog = False):
        self.hashObj = hashObj
        if editDialog is False:
            for key in self.hashObj.keys():
                if key != '':
                    obj = self.hashObj[key]
                    self.addItem(obj)
                    if hasattr(obj, "isRef") and obj.isRef == Qt.Checked:
                        self.setItemIcon(self.hashObj.index(key), 
                                         QIcon(":/reference.png"))
        else:
            for key in self.hashObj.keys():
                obj = self.hashObj[key]
                self.addItem(obj)
                if hasattr(obj, "isRef") and obj.isRef == Qt.Checked:
                    self.setItemIcon(self.hashObj.index(key), 
                                     QIcon(":/reference.png"))

    def currentObj(self):
        objName = QComboBox.currentText(self)
        if objName != "header":
            return self.hashObj[objName]
