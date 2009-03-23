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
import pyQPCR.qrc_resources

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Revision$"


class HelpDialog(QDialog):

    def __init__(self, page, parent=None):
        self.parent = parent
        QDialog.__init__(self, parent)

        backAction = QAction(QIcon(":/back.png"), "&Back", self)
        backAction.setShortcut(QKeySequence.Back)
        nextAction = QAction(QIcon(":/next.png"), "&Next", self)
        homeAction = QAction(QIcon(":/home.png"), "&Home", self)
        homeAction.setShortcut("Home")
        self.pageLabel = QLabel()

        toolBar = QToolBar()
        toolBar.addAction(backAction)
        toolBar.addAction(nextAction)
        toolBar.addAction(homeAction)
        toolBar.addWidget(self.pageLabel)
        toolBar.setIconSize(QSize(22, 22))
        self.textBrowser = QTextBrowser()

        layout = QVBoxLayout()
        layout.addWidget(toolBar)
        layout.addWidget(self.textBrowser)
        self.setLayout(layout)

        self.connect(backAction, SIGNAL("triggered()"),
                     self.textBrowser, SLOT("backward()"))
        self.connect(nextAction, SIGNAL("triggered()"),
                     self.textBrowser, SLOT("forward()"))
        self.connect(homeAction, SIGNAL("triggered()"),
                     self.textBrowser, SLOT("home()"))
        self.connect(self.textBrowser, SIGNAL("sourceChanged(QUrl)"),
                     self.updatePageTitle)

        self.textBrowser.setSearchPaths([":/help"])
        self.textBrowser.setSource(QUrl(page))
        self.resize(600, 600)
        self.setWindowTitle("%s Help" % QApplication.applicationName())

    def updatePageTitle(self):
        self.pageLabel.setText(self.textBrowser.documentTitle())

if __name__ == "__main__":
    import sys

    app = QApplication(sys.argv)
    form = HelpDialog("index.html")
    form.show()
    app.exec_()
