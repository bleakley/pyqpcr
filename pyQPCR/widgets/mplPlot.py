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

__author__ = "$Author: tgastine $"
__date__ = "$Date: 2010-10-04 09:07:07 +0200 (lun. 04 oct. 2010) $"
__version__ = "$Rev: 314 $"



class MplUnknownWidget(QWidget):
    """
    This class is used to define the relative quantification plot widget. It contains both 
    the plot and buttons allowing to customize the appearance of the plot (bars, fontsizes,
    orientation of labels, legend, ...).
    """

    def __init__(self, parent=None, barWth=0.1, barSpac=0.1, labelFt=10, labelRot=0):
        """
        Constructor of MplUnknownWidget.

        :param barWth: the width of the bars
        :type barWth: float
        :param barSpac: the spacing between bars
        :type barSpec: float
        :param labelFt: the fontsize of the labels
        :type labelFt: integer
        :param labelRot: the rotation of the labels
        :type labelRot: integer
        """
        self.barWth = barWth
        self.barSpacing = barSpac
        self.labelFontSize = labelFt
        self.labelRotation = labelRot
        QWidget.__init__(self, parent)

        vLay = QVBoxLayout()

        labPl = QLabel("&Displayed plate(s):")
        self.cboxPlate = QComboBox()
        labPl.setBuddy(self.cboxPlate)
        self.cboxPlate.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLength)
        self.cboxSens = QComboBox()
        lab1 = QLabel("&Plot axis:")
        lab1.setBuddy(self.cboxSens)
        self.cboxSens.addItems(["Target vs Sample", "Sample vs Target"])

        self.cboxScale = QComboBox()
        labScale = QLabel("Plot &scale:")
        labScale.setBuddy(self.cboxScale)
        self.cboxScale.addItems(["Linear", "Logarithmic"])

        self.spinWidth = QDoubleSpinBox()
        self.spinWidth.setLocale(QLocale(QLocale.English, 
                                         QLocale.UnitedStates))
        self.spinWidth.setValue(self.barWth)
        self.spinWidth.setRange(0.01, 0.5)
        self.spinWidth.setSingleStep(0.02)
        self.spinWidth.setVisible(False)
        self.lab2 = QLabel("Bar &width:")
        self.lab2.setBuddy(self.spinWidth)
        self.lab2.setVisible(False)

        self.spinSpacing = QDoubleSpinBox()
        self.spinSpacing.setLocale(QLocale(QLocale.English, 
                                           QLocale.UnitedStates))
        self.spinSpacing.setValue(self.barSpacing)
        self.spinSpacing.setRange(0.1, 2)
        self.spinSpacing.setSingleStep(0.1)
        self.spinSpacing.setVisible(False)

        self.labNcol = QLabel("&Legend columns")
        self.ncolLegend = QSpinBox()
        self.labNcol.setBuddy(self.ncolLegend)
        self.labNcol.setVisible(False)
        self.ncolLegend.setMaximum(5)
        self.ncolLegend.setMinimum(1)
        self.ncolLegend.setVisible(False)

        self.lab3 = QLabel("Bar &spacing:")
        self.lab3.setBuddy(self.spinSpacing)
        self.lab3.setVisible(False)
        self.btnPlot = QPushButton("&Colors and order...")
        self.lab4 = QLabel("&Font size:")
        self.cboxFontsize = QSpinBox()
        self.lab4.setBuddy(self.cboxFontsize)
        self.lab4.setVisible(False)
        self.cboxFontsize.setValue(self.labelFontSize)
        self.cboxFontsize.setRange(4, 16)
        self.cboxFontsize.setVisible(False)
        self.lab5 = QLabel("Labels &rotation:")
        self.cboxRot = QSpinBox()
        self.lab5.setBuddy(self.cboxRot)
        self.lab5.setVisible(False)
        self.cboxRot.setValue(self.labelRotation)
        self.cboxRot.setRange(0, 45)
        self.cboxRot.setSingleStep(5)
        self.cboxRot.setVisible(False)

        labLegend = QLabel("Display &legend")
        self.hideLeg = QCheckBox()
        self.hideLeg.setCheckState(Qt.Checked)
        labLegend.setBuddy(self.hideLeg)
        layLeg = QHBoxLayout()
        layLeg.addWidget(labLegend)
        layLeg.addStretch()
        layLeg.addWidget(self.hideLeg)

        line = QFrame(self)
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)

        lab = QLabel("Advanced settings...")
        self.ref = QCheckBox()
        layAdv = QHBoxLayout()
        layAdv.addWidget(lab)
        layAdv.addStretch()
        layAdv.addWidget(self.ref)

        wid = QWidget()
        vLay.addStretch()
        vLay.addWidget(labPl)
        vLay.addWidget(self.cboxPlate)
        vLay.addWidget(lab1)
        vLay.addWidget(self.cboxSens)
        vLay.addWidget(labScale)
        vLay.addWidget(self.cboxScale)
        vLay.addLayout(layLeg)
        vLay.addWidget(self.btnPlot)
        vLay.addWidget(line)
        vLay.addLayout(layAdv)
        vLay.addWidget(self.lab2)
        vLay.addWidget(self.spinWidth)
        vLay.addWidget(self.lab3)
        vLay.addWidget(self.spinSpacing)
        vLay.addWidget(self.lab4)
        vLay.addWidget(self.cboxFontsize)
        vLay.addWidget(self.lab5)
        vLay.addWidget(self.cboxRot)
        vLay.addWidget(self.labNcol)
        vLay.addWidget(self.ncolLegend)

        wid.setLayout(vLay)

        scrollArea = QScrollArea()
        scrollArea.setWidget(wid) 
        scrollArea.setSizePolicy(QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)) 
        scrollArea.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        scrollArea.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        vLay.setSizeConstraint(QLayout.SetFixedSize)

        vLayout = QVBoxLayout()
        self.mplCanUnknown = MatplotlibWidget(self, width=5,
                                              height=4, dpi=100)
        toolBar = NavToolBar(self.mplCanUnknown, self)
        vLayout.addWidget(toolBar)
        vLayout.addWidget(self.mplCanUnknown)
        hLayout = QHBoxLayout()
        hLayout.addWidget(scrollArea)
        hLayout.addLayout(vLayout)

        self.setLayout(hLayout)

        self.connect(self.ref, SIGNAL("stateChanged(int)"), self.unHide)
        self.connect(self.hideLeg, SIGNAL("stateChanged(int)"), self.hideLegend)
        self.connect(self.cboxScale, SIGNAL("activated(int)"), self.changeAxesScale)
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
        self.connect(self.ncolLegend, SIGNAL("valueChanged(int)"),
                     self.updatePlot)
        self.connect(self.cboxRot, SIGNAL("valueChanged(int)"),
                     self.changeLabelsRotation)
        self.connect(self.btnPlot, SIGNAL("clicked()"),
                     self.setPlotColor)

    def unHide(self):
        """
        This method is used to unhide the advanced configurations tools.
        """
        self.lab2.setVisible(self.ref.isChecked())
        self.lab3.setVisible(self.ref.isChecked())
        self.lab4.setVisible(self.ref.isChecked())
        self.lab5.setVisible(self.ref.isChecked())
        self.spinWidth.setVisible(self.ref.isChecked())
        self.spinSpacing.setVisible(self.ref.isChecked())
        self.cboxFontsize.setVisible(self.ref.isChecked())
        self.cboxRot.setVisible(self.ref.isChecked())
        self.labNcol.setVisible(self.ref.isChecked())
        self.ncolLegend.setVisible(self.ref.isChecked())

    def hideLegend(self):
        """
        This method is used to hide the legend of the current plot. The x limits of
        the plot are then modified.
        """
        self.leg.set_visible(self.hideLeg.isChecked())
        bool = int(self.hideLeg.checkState())/2
        xmax = 0.3 * self.xmax*bool+self.xmax + 0.05 * self.xmax*(1-bool)
        self.mplCanUnknown.axes.set_xlim((self.xmin, xmax))
        self.mplCanUnknown.draw()

    def updatePlot(self):
        """
        Simple method to update the plot.
        """
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
            # If the color or activation has changed then update
            if self.project.hashGene != dialog.hashGene:
                self.project.hashGene = dialog.hashGene
                for pl in self.project.dicoPlates.values():
                    for g in self.project.hashGene.values()[1:]:
                        if pl.dicoGene.has_key(g.name):
                            for well in pl.dicoGene[g.name]:
                                well.setGene(g)

            # If the color or activation has changed then update
            if self.project.hashEch != dialog.hashEch:
                self.project.hashEch = dialog.hashEch
                for pl in self.project.dicoPlates.values():
                    for e in self.project.hashEch.values()[1:]:
                        if pl.dicoEch.has_key(e.name):
                            for well in pl.dicoEch[e.name]:
                                well.setEch(e)
            self.plotUnknown()

    def changeFontsize(self, idraw=True):
        """
        A method to change the matplotlib axes font sizes.

        :param idraw: a boolean to indicate if the plot should be updated or not
        :type idraw: logical
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
        """
        A method to change the type of vertical scaling. It can be either
        'Linear' or 'Logarithmic'.
        """
        if self.cboxScale.currentText() == 'Linear':
            self.mplCanUnknown.axes.set_yscale('linear')
        if self.cboxScale.currentText() == 'Logarithmic':
            self.mplCanUnknown.axes.set_yscale('log')
        bool = int(self.hideLeg.checkState())/2
        xmax = 0.3*self.xmax*bool+self.xmax+0.05*self.xmax*(1-bool)
        self.mplCanUnknown.axes.set_xlim((self.xmin, xmax))
        if self.ymin not in (1e10, 0):
            self.mplCanUnknown.axes.set_ylim(ymin=self.ymin/100)
        else:
            self.mplCanUnknown.axes.set_ylim(ymin=1e-5)
        if self.cboxScale.currentText() == 'Logarithmic'and self.ymax > 0:
            self.mplCanUnknown.axes.set_ylim(ymax=self.ymax*50)
        # matplotlib 1.0.0 workaround :
        elif self.cboxScale.currentText() == 'Linear'and self.ymax > 0:
            self.mplCanUnknown.axes.set_ylim(ymax=self.ymax*1.4)
        self.mplCanUnknown.draw()

    def plotUnknown(self, project=None):
        """
        A method to plot the histograms that correspond to the relative 
        quantifications.  The errorbars displayed correspond to the 
        standard-error of NRQ. You can plot either targets vs samples, 
        or samples vs targets.

        :param project: the project you are working on
        :type project: pyQPCR.project.Project
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
        legPos = [] ; legName = [] ; xlabel = []
        dicoAbs = OrderedDict()

        self.ymin = 1e10
        self.xmax, self.ymax= 0, 0
        # Gene vs Ech
        if self.cboxSens.currentIndex() == 0:
            self.project.findBars(width, spacing, 'geneEch', platesToPlot)
            for g in self.project.hashGene.keys()[1:]:
                NRQ = [] ; NRQerror = [] ; valx = []
                for ech in self.project.hashEch.keys()[1:]:
                    for pl in platesToPlot:
                        if self.project.dicoTriplicat[pl].has_key(g) and \
                          self.project.hashGene[g].enabled and \
                          self.project.dicoTriplicat[pl][g].has_key(ech) and \
                          self.project.hashEch[ech].enabled and \
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
                    if min(NRQ) != 0:
                        self.ymin = min(self.ymin, min(NRQ))
                    self.xmax = max(self.xmax, max(valx))
                    self.ymax = max(self.ymax, max(NRQ))

        # Ech vs Gene
        elif self.cboxSens.currentIndex() == 1:
            self.project.findBars(width, spacing, 'echGene', platesToPlot)
            for ech in self.project.hashEch.keys()[1:]:
                NRQ = [] ; NRQerror = [] ; valx = []
                for g in self.project.hashGene.keys()[1:]:
                    for pl in platesToPlot:
                        if self.project.dicoTriplicat[pl].has_key(g) and \
                          self.project.hashGene[g].enabled and \
                          self.project.dicoTriplicat[pl][g].has_key(ech) and \
                          self.project.hashEch[ech].enabled and \
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
                    if min(NRQ) != 0:
                        self.ymin = min(self.ymin, min(NRQ))
                    self.xmax = max(self.xmax, max(valx))
                    self.ymax = max(self.ymax, max(NRQ))

        # plot
        self.mplCanUnknown.axes.set_xticks(self.project.barXticks.values())
        self.mplCanUnknown.axes.set_xticklabels(self.project.barXticks.keys(), 
                                                fontsize=size, 
                                                rotation=int(self.cboxRot.value()))
        # Legend + xlim
        ncol = self.ncolLegend.value()
        self.leg = self.mplCanUnknown.axes.legend(loc='upper right', 
                              shadow=True, labelspacing=0.005,
                              fancybox=True, ncol=ncol)
        legend = DraggableLegend(self.leg)
        # matplotlib 0.99.1 workaround :
        self.leg.set_axes(self.mplCanUnknown.axes)
        #
        self.leg.get_frame().set_alpha(0.2)
        self.xmin = 0
        bool = int(self.hideLeg.checkState())/2
        xmax = 0.3*self.xmax*bool+self.xmax+0.05*self.xmax*(1-bool)
        self.changeFontsize(idraw=False)
        self.changeAxesScale()
        self.mplCanUnknown.axes.set_xlim((self.xmin, xmax))

        if self.ymin not in (1e10, 0):
            self.mplCanUnknown.axes.set_ylim(ymin=self.ymin/100)
        else:
            self.mplCanUnknown.axes.set_ylim(ymin=1e-5)
        if self.cboxScale.currentText() == 'Logarithmic'and self.ymax > 0:
            self.mplCanUnknown.axes.set_ylim(ymax=self.ymax*50)
        # matplotlib 1.0.0 workaround :
        elif self.cboxScale.currentText() == 'Linear'and self.ymax > 0:
            self.mplCanUnknown.axes.set_ylim(ymax=self.ymax*1.4)

        self.mplCanUnknown.draw()
