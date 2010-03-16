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
import pyQPCR.qrc_resources
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4 import NavigationToolbar2QT
from matplotlib.figure import Figure

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"


class MatplotlibWidget(FigureCanvas):
    """
    This class allows to define a MatplotlibWidget as a classical
    PyQt widget and use it in QDialog or QMainWindow (as any other
    widget).

    @ivar fig: the matplotlib figure
    @type fig: matplotlib.figure.Figure
    @ivar axes: the matplotlib axes
    @type axes: matplotlib.axes.AxesSubplot
    """

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        """
        Constructor of the matplotlib widget

        @param parent: the parent QWidget
        @type parent: PyQt4.QtGui.QWidget
        @param width: the width of the widget
        @type width: float
        @param height: the height of the widget
        @type height: float
        @param dpi: the DPI of the matplotlib widget
        @type dpi: int
        """
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        self.fig.subplots_adjust(right=0.98, top=0.95, left=0.09)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def sizeHint(self):
        """
        size hint of the new MatplotlibWidget
        """
        w, h = self.get_width_height()
        return QSize(w, h)

    def minimumSizeHint(self):
        """
        minimum size hint of the new MatplotlibWidget
        """
        return QSize(10, 10)

class NavToolBar(NavigationToolbar2QT):
    """
    This class overloads the NavigationToolbar of matplotlib.
    The original icons are replaced by oxygen ones.
    """

    def __init__(self, canvas, parent, coordinates=True):
        """
        Constructor

        @param canvas: the matplotlib canvas
        @type canvas: matplotlib.backends.backend_qt4agg.FigureCanvasQTAgg
        @param parent: the parent QWidget
        @type parent: PyQt4.QtGui.QWidget
        @param coordinates: boolean coordinates
        @type coordinates: bool
        """
        NavigationToolbar2QT.__init__(self, canvas, parent, coordinates)
        self.setIconSize(QSize(16, 16))

    def _init_toolbar(self):
        """
        This function is a simple modification of the toolbar 
        definition in order to get the oxygen icons in the matplotib toolbar
        """
        a = self.addAction(QIcon(':/home.png'), 'Home', self.home)
        a.setToolTip('Reset original view')
        a = self.addAction(QIcon(':/undo.png'), 'Back', self.back)
        a.setToolTip('Back to previous view')
        a = self.addAction(QIcon(':/redo.png'), 'Forward', self.forward)
        a.setToolTip('Forward to next view')
        self.addSeparator()
        a = self.addAction(QIcon(':/move.png'), 'Pan', self.pan)
        a.setToolTip('Pan axes with left mouse, zoom with right')
        a = self.addAction(QIcon(':/zoom.png'), 'Zoom', self.zoom)
        a.setToolTip('Zoom to rectangle')
        self.addSeparator()
        a = self.addAction(QIcon(':/settings.png'), 'Subplots',
                self.configure_subplots)
        a.setToolTip('Configure subplots')
        a = self.addAction(QIcon(':/filesave.png'), 'Save',
                self.save_figure)
        a.setToolTip('Save the figure')

        self.buttons = {}

        if self.coordinates:
            self.locLabel = QLabel( "", self )
            self.locLabel.setAlignment(Qt.AlignRight|Qt.AlignTop )
            self.locLabel.setSizePolicy( QSizePolicy(QSizePolicy.Expanding,
                                  QSizePolicy.Ignored))
            labelAction = self.addWidget(self.locLabel)
            labelAction.setVisible(True)

        # reference holder for subplots_adjust window
        self.adj_window = None

