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
from pyQPCR.widgets.matplotlibWidget import MatplotlibWidget, NavToolBar
import pyQPCR.qrc_resources

__author__ = "$Author: tgastine $"
__date__ = "$Date: 2010-07-11 20:00:00 +0200 (dim. 11 juil. 2010) $"
__version__ = "$Rev: 282 $"



class MplStdWidget(QWidget):
    """
    This class is used to define the Standard plot widget. It contains both the plot
    and informations concerning the standard curve calculation (linear regression,
    Pearsson's coefficient and so on). It also contains a QComboBox to change the 
    displayed target.
    """

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        layout = QVBoxLayout()
        layout.addStretch()
        self.geneStdBox = QComboBox()
        lab1 = QLabel("<b>&Gene :</b>")
        lab1.setBuddy(self.geneStdBox)
        lab2 = QLabel("<b>&Linear Regression :</b>")
        self.labEquation = QLabel()
        lab2.setBuddy(self.labEquation)
        lab3 = QLabel("<b>R&#178; :</b>")
        self.labR2 =  QLabel()
        lab3.setBuddy(self.labR2)
        self.labEff =  QLabel()
        lab4 = QLabel("<b>Efficiency :</b>")
        lab4.setBuddy(self.labEff)

        layout.addWidget(lab1)
        layout.addWidget(self.geneStdBox)
        layout.addWidget(lab2)
        layout.addWidget(self.labEquation)
        layout.addWidget(lab3)
        layout.addWidget(self.labR2)
        layout.addWidget(lab4)
        layout.addWidget(self.labEff)
        layout.addStretch()

        vLayout = QVBoxLayout()
        self.mplCanStd = MatplotlibWidget(self, width=5,
                                          height=4, dpi=100)
        toolBar = NavToolBar(self.mplCanStd, self)
        vLayout.addWidget(toolBar)
        vLayout.addWidget(self.mplCanStd)

        hLayout = QHBoxLayout()
        hLayout.addLayout(layout)
        hLayout.addLayout(vLayout)

        self.setLayout(hLayout)

        self.connect(self.geneStdBox, SIGNAL("activated(int)"),
                     self.updatePlot)
                                    
    def updatePlot(self):
        """
        Update the plot
        """
        self.plotStd()

    def plotStd(self, data=None, nameGene=None):
        """
        This method allows to plot the standard curves.

        :param data: a dictionnary which contains the different plot
                     parameters
        :type data: pyQPCR.utils.odict.OrderedDict
        :param nameGene: the name of a gene we want to plot
        :type nameGene: PyQt4.QtCore.QString
        """
        if data is not None:
            self.data = data
        if nameGene is None:
            geneName = self.geneStdBox.currentText()
        else:
            geneName = nameGene
        obj = self.data[geneName]
        self.mplCanStd.axes.cla()
        self.mplCanStd.axes.scatter(obj.x, obj.y, marker='o')
        self.mplCanStd.axes.plot(obj.x, obj.yest)

        self.labEquation.setText('ct = %.2f log q0 + %.2f' \
                                % (obj.slope, obj.orig))
        self.labR2.setText('%.3f' % obj.R2)
        self.labEff.setText('%.2f%% %s %.2f' % (obj.eff, unichr(177), obj.stdeff))
        self.mplCanStd.draw()
