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
from pyQPCR.widgets.matplotlibWidget import MatplotlibWidget, NavToolBar, DraggableLegend
from pyQPCR.utils.odict import OrderedDict
from pyQPCR.dialogs.objDialog import PropDialog
import pyQPCR.qrc_resources

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"



class MplUnknownWidget(QWidget):

    def __init__(self, parent=None, barWth=0.1, barSpac=0.1, labelFt=10, labelRot=0):
        self.barWth = barWth
        self.barSpacing = barSpac
        self.labelFontSize = labelFt
        self.labelRotation = labelRot
        QWidget.__init__(self, parent)

        vLay = QVBoxLayout()
        self.cboxPlate = QComboBox()
        self.cboxSens = QComboBox()
        lab1 = QLabel("&Plot axis:")
        lab1.setBuddy(self.cboxSens)
        self.cboxSens.addItems(["Target vs Sample", "Sample vs Target"])
        self.spinWidth = QDoubleSpinBox()
        self.spinWidth.setLocale(QLocale(QLocale.English, 
                                         QLocale.UnitedStates))
        self.spinWidth.setValue(self.barWth)
        self.spinWidth.setRange(0.01, 0.5)
        self.spinWidth.setSingleStep(0.02)
        lab2 = QLabel("Bar &width:")
        lab2.setBuddy(self.spinWidth)
        self.spinSpacing = QDoubleSpinBox()
        self.spinSpacing.setLocale(QLocale(QLocale.English, 
                                           QLocale.UnitedStates))
        self.spinSpacing.setValue(self.barSpacing)
        self.spinSpacing.setRange(0.1, 2)
        self.spinSpacing.setSingleStep(0.1)
        lab3 = QLabel("Bar &spacing:")
        lab3.setBuddy(self.spinSpacing)
        self.btnPlot = QPushButton("&Colors and order...")
        lab4 = QLabel("&Font size:")
        self.cboxFontsize = QSpinBox()
        lab4.setBuddy(self.cboxFontsize)
        self.cboxFontsize.setValue(self.labelFontSize)
        self.cboxFontsize.setRange(4, 16)
        lab5 = QLabel("Labels &rotation:")
        self.cboxRot = QSpinBox()
        lab5.setBuddy(self.cboxRot)
        self.cboxRot.setValue(self.labelRotation)
        self.cboxRot.setRange(0, 45)
        self.cboxRot.setSingleStep(5)

        vLay.addStretch()
        vLay.addWidget(self.cboxPlate)
        vLay.addWidget(lab1)
        vLay.addWidget(self.cboxSens)
        vLay.addWidget(lab2)
        vLay.addWidget(self.spinWidth)
        vLay.addWidget(lab3)
        vLay.addWidget(self.spinSpacing)
        vLay.addWidget(lab4)
        vLay.addWidget(self.cboxFontsize)
        vLay.addWidget(lab5)
        vLay.addWidget(self.cboxRot)
        vLay.addWidget(self.btnPlot)
        vLay.addStretch()
        vLayout = QVBoxLayout()
        self.mplCanUnknown = MatplotlibWidget(self, width=5,
                                              height=4, dpi=100)
        self.toolBar = NavToolBar(self.mplCanUnknown, self)
        vLayout.addWidget(self.toolBar)
        vLayout.addWidget(self.mplCanUnknown)
        hLayout = QHBoxLayout()
        hLayout.addLayout(vLay)
        hLayout.addLayout(vLayout)
        self.setLayout(hLayout)

        self.connect(self.toolBar.combo, SIGNAL("activated(int)"), self.changeAxesScale)
        self.connect(self.cboxFontsize, SIGNAL("valueChanged(int)"),
                     self.changeFontsize)
        self.connect(self.cboxPlate, SIGNAL("activated(int)"),
                     self.updatePlot)
        self.connect(self.cboxSens, SIGNAL("activated(int)"),
                     self.updatePlot)
        self.connect(self.spinWidth, SIGNAL("valueChanged(double)"),
                     self.updatePlot)
        self.connect(self.spinSpacing, SIGNAL("valueChanged(double)"),
                     self.updatePlot)
        self.connect(self.cboxFontsize, SIGNAL("valueChanged(int)"),
                     self.changeFontsize)
        self.connect(self.cboxRot, SIGNAL("valueChanged(int)"),
                     self.changeLabelsRotation)
        self.connect(self.btnPlot, SIGNAL("clicked()"),
                     self.setPlotColor)

    def updatePlot(self):
        self.plotUnknown()

    def setPlotColor(self):
        """
        This method is called when the user want to change the colors
        used in the histograms. It open a widget that lets the user change
        the colors.
        """
        dialog = PropDialog(self, hashGene=self.project.hashGene,
                            hashEch=self.project.hashEch)
        if dialog.exec_():
            self.project.hashGene = dialog.hashGene
            self.project.hashEch = dialog.hashEch
            self.plotUnknown()

    def changeFontsize(self, idraw=True):
        """
        A method to change the matplotlib axes font sizes.
        """
        size = int(self.cboxFontsize.value())
        for t in self.leg.get_texts():
            t.set_fontsize(size)
        for ytick in self.mplCanUnknown.axes.get_yticklabels():
            ytick.set_fontsize(size)
        for xtick in self.mplCanUnknown.axes.get_xticklabels():
            xtick.set_fontsize(size)
        if idraw:
            self.mplCanUnknown.draw()

    def changeLabelsRotation(self):
        """
        A method to change the matplotlib xlabels orientation.
        """
        size = int(self.cboxRot.value())
        for xtick in self.mplCanUnknown.axes.get_xticklabels():
            xtick.set_rotation(size)
        self.mplCanUnknown.draw()

    def changeAxesScale(self):
        if self.toolBar.combo.currentText() == 'Linear':
            self.mplCanUnknown.axes.set_yscale('linear')
        if self.toolBar.combo.currentText() == 'Logarithmic':
            self.mplCanUnknown.axes.set_yscale('log')
        self.mplCanUnknown.draw()

    def plotUnknown(self, project=None):
        """
        A method to plot the unknown histograms
        """
        if project is not None:
            self.project = project
        if self.cboxPlate.currentText()  == 'All plates':
            platesToPlot = self.project.dicoPlates.keys()
        else:
            platesToPlot = [self.cboxPlate.currentText()]

        size = int(self.cboxFontsize.value())
        self.mplCanUnknown.axes.cla()
        width = self.spinWidth.value()
        spacing = self.spinSpacing.value()
        colors = [QColor(Qt.blue), QColor(Qt.red), QColor(Qt.green), 
                  QColor(Qt.yellow), QColor(Qt.magenta),
                  QColor(Qt.cyan), QColor(Qt.gray),
                  QColor(Qt.darkBlue), QColor(Qt.darkRed), 
                  QColor(Qt.darkGreen), QColor(Qt.darkYellow),
                  QColor(Qt.darkMagenta), QColor(Qt.darkCyan),
                  QColor(Qt.darkGray), QColor(Qt.lightGray), 
                  QColor(Qt.black)]

        # color attributions
        ind = 0
        for gene in self.project.hashGene.values()[1:]:
            if not hasattr(gene, 'color'):
                gene.setColor(colors[ind])
            if ind < len(colors)-1:
                ind += 1
            else:
                ind = 0

        ind = 0
        for ech in self.project.hashEch.values()[1:]:
            if not hasattr(ech, 'color'):
                ech.setColor(colors[ind])
            if ind < len(colors)-1:
                ind += 1
            else:
                ind = 0

        legPos = [] ; legName = [] ; xlabel = []
        dicoAbs = OrderedDict()

        xmin, ymin = 0, 0
        xmax, ymax = 0, 0
        # Gene vs Ech
        if self.cboxSens.currentIndex() == 0:
            self.project.findBars(width, spacing, 'geneEch', platesToPlot)
            for g in self.project.hashGene.keys()[1:]:
                NRQ = [] ; NRQerror = [] ; valx = []
                for ech in self.project.hashEch.keys()[1:]:
                    for pl in platesToPlot:
                        if self.project.dicoTriplicat[pl].has_key(g) and \
                          self.project.hashGene[g].enabled == Qt.Checked and \
                          self.project.dicoTriplicat[pl][g].has_key(ech) and \
                          self.project.hashEch[ech].enabled == Qt.Checked and \
                          hasattr(self.project.dicoTriplicat[pl][g][ech], 'NRQ'):
                            NRQ.append(\
                                  self.project.dicoTriplicat[pl][g][ech].NRQ)
                            NRQerror.append(\
                                  self.project.dicoTriplicat[pl][g][ech].NRQerror)
                            valx.append(self.project.barWidth[ech])
                            self.project.barWidth[ech] += width
                color = self.project.hashGene[g].color.name()
                if len(valx) != 0:
                    p = self.mplCanUnknown.axes.bar(valx, 
                            NRQ, width, color=str(color), bottom=1e-10,
                            yerr=NRQerror, ecolor='k',
                            label=str(g), align='center')
                    xmin = min(xmin, min(valx))
                    ymin = min(ymin, min(NRQ))
                    xmax = max(xmax, max(valx))
                    ymax = max(ymax, max(NRQ))

        # Ech vs Gene
        elif self.cboxSens.currentIndex() == 1:
            self.project.findBars(width, spacing, 'echGene', platesToPlot)
            for ech in self.project.hashEch.keys()[1:]:
                NRQ = [] ; NRQerror = [] ; valx = []
                for g in self.project.hashGene.keys()[1:]:
                    for pl in platesToPlot:
                        if self.project.dicoTriplicat[pl].has_key(g) and \
                          self.project.hashGene[g].enabled == Qt.Checked and \
                          self.project.dicoTriplicat[pl][g].has_key(ech) and \
                          self.project.hashEch[ech].enabled == Qt.Checked and \
                          hasattr(self.project.dicoTriplicat[pl][g][ech], 'NRQ'):
                            NRQ.append(\
                                  self.project.dicoTriplicat[pl][g][ech].NRQ)
                            NRQerror.append(\
                                  self.project.dicoTriplicat[pl][g][ech].NRQerror)
                            valx.append(self.project.barWidth[g])
                            self.project.barWidth[g] += width
                color = self.project.hashEch[ech].color.name()
                if len(valx) != 0:
                    p = self.mplCanUnknown.axes.bar(valx, 
                            NRQ, width, color=str(color), bottom=1e-10, 
                            yerr=NRQerror, ecolor='k',
                            label=str(ech), align='center')
                    xmin = min(xmin, min(valx))
                    ymin = min(ymin, min(NRQ))
                    xmax = max(xmax, max(valx))
                    ymax = max(ymax, max(NRQ))

        # plot
        self.mplCanUnknown.axes.set_xticks(self.project.barXticks.values())
        self.mplCanUnknown.axes.set_xticklabels(self.project.barXticks.keys(), 
                                                fontsize=size, 
                                                rotation=int(self.cboxRot.value()))
        # Legend + xlim
        self.leg = self.mplCanUnknown.axes.legend(loc='upper right', 
                              shadow=True, labelspacing=0.005,
                              fancybox=True)
        legend = DraggableLegend(self.leg)
        # matplotlib 0.99.1 workaround :
        self.leg.set_axes(self.mplCanUnknown.axes)
        #
        self.leg.get_frame().set_alpha(0.2)
        legendWidth = 0.3 * xmax
        self.changeFontsize(idraw=False)
        self.mplCanUnknown.axes.set_xlim((0., xmax+legendWidth))
        self.mplCanUnknown.axes.set_ylim(ymin=1.e-5)
        #self.mplCanUnknown.axes.set_ylim(ymin=0.)
        self.mplCanUnknown.draw()
