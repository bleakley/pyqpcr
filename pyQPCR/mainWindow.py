#!/usr/bin/env python
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
from pyQPCR.dialogs import *
from pyQPCR.widgets.matplotlibWidget import MatplotlibWidget, NavToolBar
from pyQPCR.widgets.tableMainWindow import PlateWidget, ResultWidget
from pyQPCR.widgets.customCbox import GeneEchComboBox
from pyQPCR.plate import Plaque, ReplicateError, PlateError
from pyQPCR.wellGeneSample import WellError
from pyQPCR.project import Project, NRQError
import matplotlib
from numpy import linspace, log10, log, sqrt, sum, mean, polyfit, polyval, \
        asarray, append, array, delete
from scipy.stats import t, norm
import os
import copy

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"
__progversion__ = "0.4"

class Qpcr_qt(QMainWindow):

    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
 
        self.setMinimumSize(640, 480)
        self.filename = None
        self.printer = None

        self.tree = QTreeWidget()
        self.tree.setHeaderLabel("Parameters")

        self.project = Project()
        self.onglet = QTabWidget()
        self.pileTables = OrderedDict()
        self.createProjWidget()
        self.onglet.addTab(self.projWidget, "Plates")
#
        self.createMplUnknownWiget()
        self.createMplStdWiget()
#
        self.createResultWidget()
        self.pileResults = OrderedDict()
#
        self.vSplitter = QSplitter(Qt.Horizontal)
        self.vSplitter.addWidget(self.tree)
        self.vSplitter.addWidget(self.onglet)
        self.mainSplitter = QSplitter(Qt.Vertical)
        self.mainSplitter.addWidget(self.vSplitter)
        self.mainSplitter.addWidget(self.resulWidget)

        self.setCentralWidget(self.mainSplitter)

        self.vSplitter.setStretchFactor(0, 1)
        self.vSplitter.setStretchFactor(1, 9)
        self.mainSplitter.setStretchFactor(0, 2)
        self.mainSplitter.setStretchFactor(1, 1)

        status = self.statusBar()

# Undo / Redo
        self.projectStack = []
        self.undoInd = -1
# Toolbar et Menus
        self.createMenusAndToolbars()

# Slots
        self.connect(self.geneComboBox, SIGNAL("activated(int)"),
                     self.modifyGene)
        self.connect(self.echComboBox, SIGNAL("activated(int)"),
                     self.modifyEch)
        self.connect(self.amComboBox, SIGNAL("activated(int)"),
                     self.modifyAm)
        self.connect(self.typeComboBox, SIGNAL("activated(int)"),
                     self.setType)
        self.connect(self.cboxSens, SIGNAL("activated(int)"),
                     self.plotUnknown)
        self.connect(self.spinWidth, SIGNAL("valueChanged(double)"),
                     self.plotUnknown)
        self.connect(self.spinSpacing, SIGNAL("valueChanged(double)"),
                     self.plotUnknown)
        self.connect(self.cboxFontsize, SIGNAL("valueChanged(int)"),
                     self.changeFontsize)
        self.connect(self.cboxRot, SIGNAL("valueChanged(int)"),
                     self.changeLabelsRotation)
        self.connect(self.btnPlot, SIGNAL("clicked()"),
                     self.setPlotColor)
        self.connect(self.geneStdBox, SIGNAL("activated(int)"),
                     self.plotStd)
        self.connect(self.tabulPlates, SIGNAL("currentChanged(int)"),
                     self.changeCurrentIndex)

# Settings pour sauvegarde de l'application
        settings = QSettings()
        self.recentFiles = settings.value("RecentFiles").toStringList()
        self.ectMax, status = settings.value("EctMax").toDouble()
        if not status:
            self.ectMax = 0.3
        self.ctMin, status = settings.value("ctMin").toDouble()
        if not status:
            self.ctMin = 35.
        self.confidence, status = settings.value("confidence").toDouble()
        self.errtype = settings.value("errtype").toString()
        self.machine = settings.value("machine").toString()
        if not status:
            self.confidence = 0.9
            self.errtype = "normal"
            self.machine = 'Eppendorf'
        geom = settings.value("Geometry").toByteArray()
        self.restoreGeometry(geom)
        self.restoreState(settings.value("MainWindow/State").toByteArray())
        self.vSplitter.restoreState(
                settings.value("VerticalSplitter").toByteArray())
        self.mainSplitter.restoreState(
                settings.value("MainSplitter").toByteArray())
        self.setWindowTitle("pyQPCR")
        self.updateFileMenu()

    def createProjWidget(self):
        self.projWidget = QWidget(self)
        vLay = QVBoxLayout()
        self.tabulPlates = QTabWidget()
        self.tabulPlates.setTabPosition(QTabWidget.West)
        vLay.addWidget(self.tabulPlates)
        self.projWidget.setLayout(vLay)

    def createResultWidget(self):
        self.resulWidget = QWidget(self)
        vLay = QVBoxLayout()
        self.tabulResults = QTabWidget()
        self.tabulResults.setTabPosition(QTabWidget.West)
        vLay.addWidget(self.tabulResults)
        self.resulWidget.setLayout(vLay)

    def createMplUnknownWiget(self):
        self.plotUnknownWidget = QWidget()
        vLay = QVBoxLayout()
        self.cboxSens = QComboBox()
        lab1 = QLabel("&Plot axis:")
        lab1.setBuddy(self.cboxSens)
        self.cboxSens.addItems(["Target vs Sample", "Sample vs Target"])
        self.spinWidth = QDoubleSpinBox()
        self.spinWidth.setLocale(QLocale(QLocale.English, 
                                         QLocale.UnitedStates))
        self.spinWidth.setValue(0.1)
        self.spinWidth.setRange(0.01, 0.5)
        self.spinWidth.setSingleStep(0.02)
        lab2 = QLabel("Bar &width:")
        lab2.setBuddy(self.spinWidth)
        self.spinSpacing = QDoubleSpinBox()
        self.spinSpacing.setLocale(QLocale(QLocale.English, 
                                           QLocale.UnitedStates))
        self.spinSpacing.setValue(0.3)
        self.spinSpacing.setRange(0.1, 2)
        self.spinSpacing.setSingleStep(0.1)
        lab3 = QLabel("Bar &spacing:")
        lab3.setBuddy(self.spinSpacing)
        self.btnPlot = QPushButton("&Colors and order...")
        lab4 = QLabel("&Font size:")
        self.cboxFontsize = QSpinBox()
        lab4.setBuddy(self.cboxFontsize)
        self.cboxFontsize.setValue(10)
        self.cboxFontsize.setRange(4, 16)
        lab5 = QLabel("Labels &rotation:")
        self.cboxRot = QSpinBox()
        lab5.setBuddy(self.cboxRot)
        self.cboxRot.setValue(0)
        self.cboxRot.setRange(0, 45)
        self.cboxRot.setSingleStep(5)

        vLay.addStretch()
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
        self.mplCanUnknown = MatplotlibWidget(self.plotUnknownWidget, width=5,
                                              height=4, dpi=100)
        toolBar = NavToolBar(self.mplCanUnknown, self)
        vLayout.addWidget(toolBar)
        vLayout.addWidget(self.mplCanUnknown)
        hLayout = QHBoxLayout()
        hLayout.addLayout(vLay)
        hLayout.addLayout(vLayout)
        self.plotUnknownWidget.setLayout(hLayout)

    def createMplStdWiget(self):
        layout = QVBoxLayout()
        layout.addStretch()
        self.geneStdBox = QComboBox()
        lab1 = QLabel("<b>&Gene:</b>")
        lab1.setBuddy(self.geneStdBox)
        lab2 = QLabel("<b>&Linear Regression:</b>")
        self.labEquation = QLabel()
        lab2.setBuddy(self.labEquation)
        lab3 = QLabel("<b>R^2:</b>")
        self.labR2 =  QLabel()
        lab3.setBuddy(self.labR2)
        self.labEff =  QLabel()
        lab4 = QLabel("<b>Efficiency:</b>")
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

        self.plotStdWidget = QWidget()
        vLayout = QVBoxLayout()
        self.mplCanStd = MatplotlibWidget(self.plotUnknownWidget, width=5, 
                                          height=4, dpi=100)
        toolBar = NavToolBar(self.mplCanStd, self)
        vLayout.addWidget(toolBar)
        vLayout.addWidget(self.mplCanStd)

        hLayout = QHBoxLayout()
        hLayout.addLayout(layout)
        hLayout.addLayout(vLayout)
        self.plotStdWidget.setLayout(hLayout)

    def createMenusAndToolbars(self):
        fileOpenAction = self.createAction("&Open...", self.fileOpen, 
                QKeySequence.Open, "fileopen", "Open an existing project")
        fileNewAction = self.createAction("&New project...", self.fileNew, 
                QKeySequence.New, "filenew", "Create a new project")
        self.fileImportAction = self.createAction("&Import...", self.fileImport, 
                "Ctrl+I", "fileimport", "Import/add an existing plate")
        self.closeTabAction = self.createAction("&Close a plate", self.closePlate, 
                QKeySequence.Close, "closeplate", "Close an existing plate")
        self.filePrintAction = self.createAction("&Print", self.filePrint,
                QKeySequence.Print, "fileprint", "Print results")
        self.exportAction = self.createAction("&Export as PDF", self.fileExport,
                "Ctrl+D", "pdf", "Export results in a PDF file")
        self.fileSaveAction = self.createAction("&Save", self.fileSave,
                QKeySequence.Save, "filesave", "Save the file")
        self.fileSaveAction.setEnabled(False)
        self.fileSaveAsAction = self.createAction("Save &As...",
                self.fileSaveAs, icon="filesaveas",
                tip="Save the file using a new name")
        fileQuitAction = self.createAction("&Quit", self.close, 
                "Ctrl+Q", "filequit", "Close the application")
        self.editAction = self.createAction("Edit wells", self.editWell, 
                "Ctrl+E", "edit", "Edit selected wells")
        self.undoAction = self.createAction("Undo", self.undo, 
                QKeySequence.Undo, "undo", "Undo")
        self.redoAction = self.createAction("Redo", self.redo,
                QKeySequence.Redo, "redo", "Redo")
        self.addGeneAction = self.createAction("Add &Target...", self.addGene,
                "Ctrl+T", "addgene", "Add a new target")
        self.addEchAction = self.createAction("Add &Sample...", self.addEch,
                "Ctrl+G", "addech", "Add a new sample")
        self.addAmAction = self.createAction("Add A&mount...", self.addAmount,
                "Ctrl+M", "addamount", "Add a new amount")
        self.plotAction = self.createAction("Quantifications", 
                             self.computeUnknown, "Ctrl+Shift+U", 
                             "plotUnknown", "Plot results")
        self.plotStdAction = self.createAction("Standard curves", 
                              self.computeStd, "Ctrl+Shift+S", 
                              "plotStandard", "Plot standard curves")
        self.enableAction = self.createAction("Enable wells", self.enable, 
                     None, "enable", "Enable selected wells")
        self.disableAction = self.createAction("Disable wells", self.disable,
                     None, "disable", "Disable selected wells")
        settingsAction = self.createAction("&Configure pyQPCR...", 
                                           self.configure, icon="settings")
        helpAboutAction = self.createAction("&About pyQPCR", self.helpAbout,
                icon="about")
        helpHelpAction = self.createAction("&Help", self.helpHelp,
                QKeySequence.HelpContents, icon="help")
# Menus
        fileMenu = self.menuBar().addMenu("&File")
        self.addActions(fileMenu, (fileOpenAction, fileNewAction, 
                                   self.fileImportAction, self.closeTabAction))
        self.recentFileMenu = fileMenu.addMenu(QIcon(":/filerecent.png"),
                "Open recent files")
        fileMenu.addSeparator()
        self.addActions(fileMenu, (self.filePrintAction, self.exportAction, None, 
                        self.fileSaveAction, self.fileSaveAsAction, 
                        None, fileQuitAction))
        editMenu = self.menuBar().addMenu("&Edit")
        editMenu.addAction(self.editAction)
        editMenu.addSeparator()
        self.addActions(editMenu, (self.undoAction, self.redoAction))
        editMenu.addSeparator()
        self.addActions(editMenu, (self.addEchAction, self.addGeneAction,
                                   self.addAmAction))
        calculMenu = self.menuBar().addMenu("&Computations")
        self.addActions(calculMenu, (self.enableAction, self.disableAction,
                                  None, self.plotStdAction, self.plotAction))
        settingsMenu = self.menuBar().addMenu("&Settings")
        settingsMenu.addAction(settingsAction)
        helpMenu = self.menuBar().addMenu("&Help")
        self.addActions(helpMenu, (helpAboutAction, helpHelpAction))

# Le menu doit afficher les fichiers recemment ouverts
        self.connect(self.recentFileMenu, SIGNAL("aboutToShow()"),
                self.updateFileMenu)
# Toolbars
        fileToolbar = self.addToolBar("File")
        fileToolbar.setObjectName("FileToolBar")
        self.addActions(fileToolbar, (fileOpenAction, fileNewAction,
                        self.fileImportAction, self.closeTabAction, 
                        self.filePrintAction, self.exportAction, 
                        self.fileSaveAction, self.fileSaveAsAction))
        fileToolbar.setIconSize(QSize(22, 22))

        editToolbar = self.addToolBar("Edit")
        editToolbar.setObjectName("Edit ToolBar")
        self.addActions(editToolbar, (self.undoAction, self.redoAction))
        editToolbar.addSeparator()
        self.typeComboBox = QComboBox()
        self.typeComboBox.addItems(["unknown", "standard", "negative"])
        pix = QPixmap(32,32)
        bleu = QColor(116, 167, 227) ; pix.fill(bleu)
        self.typeComboBox.setItemIcon(0, QIcon(pix))
        rouge = QColor(233, 0, 0) ; pix.fill(rouge)
        self.typeComboBox.setItemIcon(1, QIcon(pix))
        jaune = QColor(255, 250, 80) ; pix.fill(jaune)
        self.typeComboBox.setItemIcon(2, QIcon(pix))
        self.typeComboBox.setToolTip("List of types")
        self.typeComboBox.setStatusTip(self.typeComboBox.toolTip())
        self.typeComboBox.setFocusPolicy(Qt.NoFocus)
        editToolbar.addWidget(self.typeComboBox)
        editToolbar.addSeparator()
        self.geneComboBox = GeneEchComboBox()
        self.geneComboBox.setToolTip("List of targets")
        self.geneComboBox.setStatusTip(self.geneComboBox.toolTip())
        self.geneComboBox.setFocusPolicy(Qt.NoFocus)
        self.echComboBox = GeneEchComboBox()
        self.echComboBox.setToolTip("List of samples")
        self.echComboBox.setStatusTip(self.echComboBox.toolTip())
        self.echComboBox.setFocusPolicy(Qt.NoFocus)
        self.amComboBox = GeneEchComboBox()
        self.amComboBox.setToolTip("List of amounts")
        self.amComboBox.setStatusTip(self.amComboBox.toolTip())
        self.amComboBox.setFocusPolicy(Qt.NoFocus)
        editToolbar.addWidget(self.geneComboBox)
        editToolbar.addAction(self.addGeneAction)
        editToolbar.addSeparator()
        editToolbar.addWidget(self.echComboBox)
        editToolbar.addAction(self.addEchAction)
        editToolbar.addSeparator()
        editToolbar.addWidget(self.amComboBox)
        editToolbar.addAction(self.addAmAction)
        editToolbar.addSeparator()
        editToolbar.addAction(self.editAction)
        editToolbar.setIconSize(QSize(22, 22))

        plotToolbar = self.addToolBar("Plot")
        plotToolbar.setObjectName("PlotToolBar")
        self.addActions(plotToolbar, (self.plotStdAction, self.plotAction))
        plotToolbar.setIconSize(QSize(22, 22))
# ContextMenu
        self.projWidget.setContextMenuPolicy(Qt.ActionsContextMenu)
        self.resulWidget.setContextMenuPolicy(Qt.ActionsContextMenu)
        self.addActions(self.projWidget, (fileOpenAction, fileNewAction, 
                       self.fileImportAction, self.closeTabAction, 
                       self.editAction, self.undoAction, self.redoAction, 
                       self.addGeneAction, self.addEchAction, 
                       self.enableAction, self.disableAction))
        self.addActions(self.resulWidget, (self.filePrintAction, self.exportAction,
                                      self.fileSaveAction, self.fileSaveAsAction,
                                      self.plotStdAction, self.plotAction))
# Desactivation par defaut
        self.activateDesactivate(False)

    def createAction(self, text, slot=None, shortcut=None, icon=None,
                     tip=None, checkable=None, signal="triggered()"):
        """
        This method highly simplifies the creation of QAction
        """
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut  is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action

    def addActions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def updateStatus(self, message, time=5000):
        self.statusBar().showMessage(message, time)

    def fileOpen(self):
        if not self.okToContinue():
            return
        dir = os.path.dirname(self.filename) if self.filename is not None \
                else "."
        formats =[u"*.xml"]
        fname = unicode(QFileDialog.getOpenFileName(self,
                                                    "pyQPCR - Choose a file", 
                dir, "Input files (%s)" % " ".join(formats)))
        if fname:
            self.loadFile(fname)

    def loadFile(self, fname=None):
        if fname is None:
            action = self.sender()
            if isinstance(action, QAction):
                fname = unicode(action.data().toString())
                if not self.okToContinue():
                    return
        if fname:
            self.setWindowTitle("pyQPCR - %s[*]" % QFileInfo(fname).fileName())
            self.addRecentFile(fname)
            self.filename = fname
            message = "Loaded %s" % QFileInfo(fname).fileName()
            self.updateStatus(message)

            self.cleanBeforeOpen()

            self.project = Project(fname)

            for key in self.project.dicoPlates.keys():
                pl = self.project.dicoPlates[key]
                self.appendPlate(pl, key)
                self.appendResult(pl, key)

# Pile de plaques pour le Undo/Redo
            self.projectStack.append(copy.deepcopy(self.project))
            self.updateUi()

    def fileNew(self):
        if not self.okToContinue():
            return
        dir = os.path.dirname(self.filename) if self.filename is not None \
                else "."
        dialog = NewProjectDialog(self, pwd=dir)
        if dialog.exec_():
            self.cleanBeforeOpen()
            self.project = Project(dialog.projectName)
            self.projectStack.append(copy.deepcopy(self.project))
            self.filename = dialog.projectName
            self.setWindowTitle("pyQPCR - %s[*]" % QFileInfo(self.filename).fileName())
            for fname in dialog.fileNames.values():
                self.machine = dialog.machineType
                self.addPlate(fname)

    def fileImport(self):
        dir = os.path.dirname(self.filename) if self.filename is not None \
                else "."
        if self.machine == 'Eppendorf':
            formats =[u"*.txt", u"*.csv"]
            type = 'Eppendorf machines'
        elif self.machine == 'Applied StepOne':
            formats =[u"*.txt"]
            type = 'Applied StepOne machines'
        elif self.machine == 'Applied 7000':
            formats =[u"*.csv"]
            type = 'Applied 7000 machines'
        else:
            formats =[u"*.txt", u"*.csv"]
            type = 'Eppendorf machines'
        fileNames = QFileDialog.getOpenFileNames(self,
                       "pyQPCR - Choose a file", dir, 
                       "Input files [%s] (%s)" % (type, " ".join(formats)))
        if fileNames:
            for file in fileNames:
                if not self.project.dicoPlates.has_key(QFileInfo(file).fileName()):
                    self.addPlate(file)
                else:
                    name = QFileInfo(file).fileName()
                    QMessageBox.warning(self, "Import cancelled",
                      "<b>Warning</b>: the plate %s is already in the project %s ! \
                       As a consequence, it has not been added to the project."       
                                    % (name, QFileInfo(self.filename).fileName()))

    def addPlate(self, fname=None):
        if fname is None:
            action = self.sender()
            if isinstance(action, QAction):
                fname = unicode(action.data().toString())
                if not self.okToContinue():
                    return
        if fname:
            message = "Loaded %s" % QFileInfo(fname).fileName()
            self.updateStatus(message)
# Nettoyage des tableaux avant l'eventuel remplissage
            try:
                plaque = Plaque(fname, self.machine)
                self.project.addPlate(plaque)
                key = QFileInfo(fname).fileName()

                self.appendPlate(plaque, key)
                self.appendResult(plaque, key)

                self.updateUi()
                self.project.unsaved = True
                self.fileSaveAction.setEnabled(True)
                self.projectStack.append(copy.deepcopy(self.project))
            except KeyError:
                st = "<b>Warning:</b>: an error occurs during import ! "
                st += "It probably comes from your file which has not been correctly parsed. "
                st += "Is it a raw <b>%s</b> file ?" % self.machine
                QMessageBox.warning(self, "Warning file import", st)
            except PlateError, e:
                QMessageBox.warning(self, "Warning file import", str(e))

    def closePlate(self):
        reply = QMessageBox.question(self, "Remove a plate",
                            "Are you sure to remove %s ?" % self.currentPlate,
                            QMessageBox.Yes|QMessageBox.No)
        plToDestroy = self.currentPlate
        if reply == QMessageBox.Yes:
            message = "Closed %s" % plToDestroy
            index = self.project.dicoPlates.index(plToDestroy)
            self.project.removePlate(plToDestroy)
# Nettoyage des onglets et des tableaux
            self.tabulPlates.removeTab(index)
            self.pileTables.__delitem__(plToDestroy)
            self.tabulResults.removeTab(index)
            self.pileResults.__delitem__(plToDestroy)
# Maj des cbox et Remplissage du tree
            self.updateUi()
            self.project.unsaved = True
            self.fileSaveAction.setEnabled(True)
            self.projectStack.append(copy.deepcopy(self.project))

    def activateDesactivate(self, bool):
        """
        This method allows to enable/disable several QAction

        @param bool: the boolean value
        @type bool: bool
        """
        self.addGeneAction.setEnabled(bool)
        self.addEchAction.setEnabled(bool)
        self.addAmAction.setEnabled(bool)
        self.editAction.setEnabled(bool)
        self.plotAction.setEnabled(bool)
        self.plotStdAction.setEnabled(bool)
        self.typeComboBox.setEnabled(bool)
        self.geneComboBox.setEnabled(bool)
        self.echComboBox.setEnabled(bool)
        self.amComboBox.setEnabled(bool)
        self.fileSaveAsAction.setEnabled(bool)
        self.filePrintAction.setEnabled(bool)
        self.exportAction.setEnabled(bool)
        self.undoAction.setEnabled(bool)
        self.redoAction.setEnabled(bool)
        self.enableAction.setEnabled(bool)
        self.disableAction.setEnabled(bool)
        self.fileImportAction.setEnabled(bool)
        self.closeTabAction.setEnabled(bool)

    def populateTree(self):
        self.tree.clear()
        if self.filename is not None:
            ancestor = QTreeWidgetItem(self.tree, 
                                       [QFileInfo(self.filename).fileName()])
        else:
            ancestor = QTreeWidgetItem(self.tree, ["untitled project"])
        for key in self.project.dicoPlates.keys():
            item = QTreeWidgetItem(ancestor, [key])
            pl = self.project.dicoPlates[key]
            itemQuant = QTreeWidgetItem(item, ["Quantification"])
            itemRefGene = QTreeWidgetItem(itemQuant , ["Reference Target"])
            itemRefEch = QTreeWidgetItem(itemQuant , ["Reference Sample"])
            item = QTreeWidgetItem(itemRefGene, [pl.geneRef])
            item = QTreeWidgetItem(itemRefEch, [pl.echRef])
        itemStd = QTreeWidgetItem(ancestor, ["Standard"])
        for gene in self.project.hashGene.values()[1:]:
            eff = "%s:%.2f%s%.2f" % (gene.name, gene.eff, unichr(177), gene.pm)
            item = QTreeWidgetItem(itemStd, [eff])
        self.tree.expandAll()

    def populateCbox(self, cbox, items, name="Target"):
        cbox.clear()
        cbox.addItem(name, QVariant("header"))
        cbox.addItems(items)

    def fileSave(self):
        self.project.exportXml(self.filename)
        self.updateStatus("Saved %s" % self.filename)
        self.project.unsaved = False
        self.fileSaveAction.setEnabled(False)

    def fileSaveAs(self):
        formats =[u"*.xml"]
        fname = self.filename if self.filename is not None else "."
        fname = unicode(QFileDialog.getSaveFileName(self, 
                "pyQPCR - Save a file", fname,
                "Result files (%s)" % " ".join(formats)))
        if fname:
            self.setWindowTitle("pyQPCR - %s[*]" % QFileInfo(fname).fileName())
            self.filename = fname
            self.fileSave()

    def generateHTML(self):
        dialog = PrintingDialog(self)
        if dialog.exec_():
            isTable = dialog.btnRes.isChecked()
            isStd = dialog.btnStd.isChecked()
            isQuant = dialog.btnQuant.isChecked()
        else:
            return
        if not hasattr(self, "project"):
            return
        html = u""
        css = ("<html><head>\n"
               '<style type="text/css">\n'
               "table {border-color:black; border-style:solid;}\n"
               "th, td {font-size:6pt;}"
               "</style>\n"
               "</head>\n")
        html += css
        html += "<h1 align=center> qPCR results </h1><br><br>"
        if isTable:
            for key in self.project.dicoPlates.keys():
                html += "<br><h2>Results table (%s)</h2><br>" % key
                html += self.project.dicoPlates[key].writeHtml()
        if isStd and self.nplotStd !=0:
            html += "<p style='page-break-before:always;'>"
            html += "<br><h2>Standard curves</h2><br>"
            self.geneStdBox.addItems(self.project.dicoStd.keys())
            for index in range(len(self.project.dicoStd.keys())):
                if (not index % 3 and index !=0):
                    html += "<p style='page-break-before:always;'>"
                else:
                    html += "<p>"
                html += "<table border=0 width=100%>\n"
                self.geneStdBox.setCurrentIndex(index)
                self.plotStd()
                fig = self.mplCanStd.figure.savefig("output%i.png" % index, 
                                                    dpi=100)
                html += "<tr valign=middle>\n"
                html += ("<th align=center>"
                         "<table width=100% border=0>")
                html += "<tr><th><font size 10pt><b>Gene:</b> %s</th></tr>" %\
                        self.geneStdBox.currentText()
                html += "<tr><th><b>Linear Regression:</b> %s</th></tr>" %\
                        self.labEquation.text() 
                html += "<tr><th><b>Efficiency:</b> %s</th></tr>" % \
                        self.labEff.text()
                html += "<tr><th><b>R<SUP>2</SUP>:</b> %s</th></tr>" % \
                        self.labR2.text()
                html += ("</table>"
                         "</th>")

                html += "<th align=center>"
                html += "<p><img src='output%i.png' width=300></p>" % \
                        index
                html += ("</th>"
                         "</tr>\n")
                html += "</table>"
                html += "</p>"
            html += "</p>"
        if isQuant:
            html += "<br><h2>Quantification curves</h2>"
            fig = self.mplCanUnknown.figure.savefig("output.png", dpi=100)
            html += "<p><img src='output.png' width=500></p>"
        html += "</html>"
        return html

    def filePrint(self):
        """
        A method to print results thanks to a QPrinter object
        """
        html = self.generateHTML()
        if html is None:
            return
        if self.printer is None:
            self.printer = QPrinter(QPrinter.HighResolution)
            self.printer.setPageSize(QPrinter.A4)
        form = QPrintDialog(self.printer, self)
        if form.exec_():
            document = QTextDocument()
            document.setHtml(html)
            document.print_(self.printer)
            self.cleanPngs()

    def fileExport(self):
        """
        A method to export results in a PDF file
        """
        html = self.generateHTML()
        if html is None:
            return
        if self.printer is None:
            self.printer = QPrinter(QPrinter.HighResolution)
            self.printer.setPageSize(QPrinter.A4)
            self.printer.setOutputFormat(QPrinter.PdfFormat)
        formats =[u"*.pdf"]
        fname = unicode(QFileDialog.getSaveFileName(self, 
                "pyQPCR - Export results in a pdf file", "results.pdf",
                "PDF - Portable Document Format (%s)" % " ".join(formats)))
        if fname:
            self.printer.setOutputFileName(fname)
            document = QTextDocument()
            document.setHtml(html)
            document.print_(self.printer)
            self.cleanPngs()

    def cleanPngs(self):
        """
        Remove the png files needed after print or export
        """
        for file in os.listdir('.'):
            if file.startswith('output'):
                os.remove(file)

    def configure(self):
        dialog = SettingsDialog(self, ect=self.ectMax,
                                ctmin=self.ctMin,
                                confidence=self.confidence,
                                errtype=self.errtype,
                                machine=self.machine)
        if dialog.exec_():
            self.ectMax, st = dialog.ectLineEdit.text().toFloat()
            self.ctMin, st = dialog.ctMinLineEdit.text().toFloat()
            self.confidence, st = dialog.confCbox.currentText().toFloat()
            errtype = dialog.typeCbox.currentText()
            self.machine = dialog.machBox.currentText()
            self.errtype = dialog.types[errtype]
            self.confidence /= 100

    def helpAbout(self):
        import platform
        QMessageBox.about(self, "About pyQPCR",
        """<b>pyQPCR</b> v %s
        <p>Copyright &copy; 2008- Thomas Gastine - Magali Hennion
        All rights reserved.
        <p>This application can be used to perform
        simple qPCR analysis.
        <p> It may be used, copied and modified with no restriction
        <p>Python %s - PyQt %s - Matplotlib %s
        on %s""" % (__progversion__, platform.python_version(),
        PYQT_VERSION_STR, matplotlib.__version__, platform.system()))

    def helpHelp(self):
        f = HelpDialog("index.html", self)
        f.show()

    def addRecentFile(self, fname):
        if fname is None:
            return
        if not self.recentFiles.contains(fname):
            self.recentFiles.prepend(QString(fname))
            while self.recentFiles.count() > 9:
                self.recentFiles.takeLast()

    def updateFileMenu(self):
        self.recentFileMenu.clear()
        current = QString(self.filename) if self.filename is not None else None
        recentFiles = []
        for fname in self.recentFiles:
            if fname != current and QFile.exists(fname):
                recentFiles.append(fname)
        if recentFiles:
            for i, fname in enumerate(recentFiles):
                action = QAction(QIcon(":/filerecent.png"), "&%d %s" % \
                        (i+1, QFileInfo(fname).fileName()), self)
                action.setData(QVariant(fname))
                self.connect(action, SIGNAL("triggered()"), self.loadFile)
                self.recentFileMenu.addAction(action)

    def closeEvent(self, event):
# Widget de confirmation
        if self.okToContinue():
            settings = QSettings()
            filename = QVariant(QString(self.filename)) \
                  if self.filename is not None else QVariant()
            settings.setValue("Last File", filename)
            recentFiles = QVariant(self.recentFiles) if self.recentFiles \
                  else QVariant()
            settings.setValue("RecentFiles", recentFiles)
            ectMax = QVariant(self.ectMax) if self.ectMax \
                  else QVariant()
            settings.setValue("EctMax", ectMax)
            ctMin = QVariant(self.ctMin) if self.ctMin \
                  else QVariant()
            settings.setValue("ctMin", ctMin)
            confidence = QVariant(self.confidence) if self.confidence \
                  else QVariant()
            errtype = QVariant(self.errtype) if self.errtype \
                  else QVariant()
            machine = QVariant(self.machine) if self.machine \
                  else QVariant()
            settings.setValue("confidence", confidence)
            settings.setValue("errtype", errtype)
            settings.setValue("machine", machine)
            settings.setValue("Geometry", QVariant(self.saveGeometry()))
            settings.setValue("MainWindow/State", QVariant(self.saveState()))
            settings.setValue("VerticalSplitter", 
                    QVariant(self.vSplitter.saveState()))
            settings.setValue("MainSplitter", 
                    QVariant(self.mainSplitter.saveState()))
        else:
            event.ignore()

    def okToContinue(self):
        if hasattr(self, "project"):
            if self.project.unsaved:
                reponse = QMessageBox.question(self,
                        "pyQPCR - Unsaved Changes",
                        "Save unsaved changes?",
                        QMessageBox.Yes|QMessageBox.No|QMessageBox.Cancel)
                if reponse == QMessageBox.Cancel:
                    return False
                elif reponse == QMessageBox.Yes:
                    self.fileSave()
        return True

    def redo(self):
        if self.undoInd < -1:
            self.undoInd += 1
        self.project = copy.deepcopy(self.projectStack[self.undoInd])
        self.updateUi()
#
        for key in self.project.dicoPlates.keys():
            pl = self.project.dicoPlates[key]
            if not self.pileTables.has_key(key):
                self.appendPlate(pl, key)
                self.appendResult(pl, key)
            self.pileTables[key].populateTable(pl)
            self.pileResults[key].populateResult(pl)

        for key in self.pileTables.keys():
            if not self.project.dicoPlates.has_key(key):
                index = self.pileTables.index(key)
                self.tabulPlates.removeTab(index)
                self.tabulResults.removeTab(index)
                self.pileTables.__delitem__(key)
                self.pileResults.__delitem__(key)

        if self.undoInd == 0 or self.undoInd == -len(self.projectStack):
            self.project.unsaved = False
            self.fileSaveAction.setEnabled(False)
        else:
            self.project.unsaved = True
            self.fileSaveAction.setEnabled(True)

    def undo(self):
        if abs(self.undoInd) < abs(len(self.projectStack)):
            self.undoInd -= 1
        self.project = copy.deepcopy(self.projectStack[self.undoInd])
        self.updateUi()
#
        for key in self.project.dicoPlates.keys():
            pl = self.project.dicoPlates[key]
            if not self.pileTables.has_key(key):
                self.appendPlate(pl, key)
                self.appendResult(pl, key)
            self.pileTables[key].populateTable(pl)
            self.pileResults[key].populateResult(pl)

        for key in self.pileTables.keys():
            if not self.project.dicoPlates.has_key(key):
                index = self.pileTables.index(key)
                self.tabulPlates.removeTab(index)
                self.tabulResults.removeTab(index)
                self.pileTables.__delitem__(key)
                self.pileResults.__delitem__(key)

# Si on remonte tous les undo pas besoin de sauvegarder
        if self.undoInd == 0 or self.undoInd == -len(self.projectStack):
            self.project.unsaved = False
            self.fileSaveAction.setEnabled(False)
        else:
            self.project.unsaved = True
            self.fileSaveAction.setEnabled(True)

    def editWell(self):
        setType = set()
        setGene = set()
        setEch = set()
        setAm = set()
        selected = [setType, setGene, setEch, setAm]
        for it in self.pileTables[self.currentPlate].selectedItems():
            nom = it.statusTip()
            well = getattr(self.project.dicoPlates[self.currentPlate], str(nom))
            setType.add(well.type)
            setEch.add(well.ech.name)
            setGene.add(well.gene.name)
            setAm.add(well.amount)
        dialog = EditDialog(self, project=self.project, selected=selected)
        #
        if dialog.exec_():
            ge = dialog.cboxGene.currentObj()
            if dialog.cboxType.currentText() == QString('unknown'):
                ech = dialog.cboxSample.currentObj()
                for it in self.pileTables[self.currentPlate].selectedItems():
                    nom = it.statusTip()
                    well = getattr(self.project.dicoPlates[self.currentPlate],
                                   str(nom))
                    well.setType(QString('unknown'))
                    if ge.name != "": well.setGene(ge)
                    if ech.name != "": well.setEch(ech)
                for pl in self.project.dicoPlates.values():
                    pl.setDicoEch()
            if dialog.cboxType.currentText() == QString('standard'):
                am = dialog.cboxAm.currentText()
                for it in self.pileTables[self.currentPlate].selectedItems():
                    nom = it.statusTip()
                    well = getattr(self.project.dicoPlates[self.currentPlate],
                                   str(nom))
                    well.setType(QString('standard'))
                    if ge.name != '':
                        well.setGene(ge)
                    if am != '':
                        well.setAmount(float(am))
                self.project.setDicoAm()
            self.project.unsaved = True
            self.fileSaveAction.setEnabled(True)
            for pl in self.project.dicoPlates.values():
                pl.setDicoGene()
            self.projectStack.append(copy.deepcopy(self.project))
            self.pileTables[self.currentPlate].populateTable( \
                   self.project.dicoPlates[self.currentPlate])
            self.pileResults[self.currentPlate].populateResult( \
                   self.project.dicoPlates[self.currentPlate])

    def addGene(self):
        dialog = GeneDialog(self, project=self.project)
        if dialog.exec_():
            project = dialog.project
            self.populateCbox(self.geneComboBox, project.hashGene, "Target")
            self.project = project
            self.fileSaveAction.setEnabled(self.project.unsaved)
            for pl in self.project.dicoPlates.values():
                pl.setDicoGene()
            self.projectStack.append(copy.deepcopy(self.project))
            for key in self.project.dicoPlates.keys():
                pl = self.project.dicoPlates[key]
                self.pileTables[key].populateTable(pl)
                self.pileResults[key].populateResult(pl)
            self.populateTree()

    def addEch(self):
        dialog = EchDialog(self, project=self.project)
        if dialog.exec_():
            project = dialog.project
            self.populateCbox(self.echComboBox, project.hashEch, "Sample")
            self.project = project
            self.fileSaveAction.setEnabled(self.project.unsaved)
            for pl in self.project.dicoPlates.values():
                pl.setDicoEch()
            self.projectStack.append(copy.deepcopy(self.project))
            for key in self.project.dicoPlates.keys():
                pl = self.project.dicoPlates[key]
                self.pileTables[key].populateTable(pl)
                self.pileResults[key].populateResult(pl)
            self.populateTree()

    def addAmount(self):
        dialog = AmountDialog(self, project=self.project)
        if dialog.exec_():
            project = dialog.project
            self.populateCbox(self.amComboBox, project.hashAmount, "Amount")
            self.project = project
            self.fileSaveAction.setEnabled(self.project.unsaved)
            self.project.setDicoAm()
            self.projectStack.append(copy.deepcopy(self.project))
            for key in self.project.dicoPlates.keys():
                pl = self.project.dicoPlates[key]
                self.pileTables[key].populateTable(pl)
                self.pileResults[key].populateResult(pl)
            self.populateTree()

    def modifyGene(self):
        for it in self.pileTables[self.currentPlate].selectedItems():
            gene = self.geneComboBox.currentObj()
            nom = it.statusTip()
            well = getattr(self.project.dicoPlates[self.currentPlate], str(nom))
            well.setGene(gene)
        self.project.unsaved = True
        self.fileSaveAction.setEnabled(True)
        self.project.dicoPlates[self.currentPlate].setDicoGene()
        self.projectStack.append(copy.deepcopy(self.project))
        self.pileTables[self.currentPlate].populateTable( \
                      self.project.dicoPlates[self.currentPlate])
        self.pileResults[self.currentPlate].populateResult( \
                      self.project.dicoPlates[self.currentPlate])

    def modifyEch(self):
        for it in self.pileTables[self.currentPlate].selectedItems():
            ech = self.echComboBox.currentObj()
            nom = it.statusTip()
            well = getattr(self.project.dicoPlates[self.currentPlate], str(nom))
            well.setEch(ech)
        self.project.unsaved = True
        self.fileSaveAction.setEnabled(True)
        self.project.dicoPlates[self.currentPlate].setDicoEch()
        self.projectStack.append(copy.deepcopy(self.project))
        self.pileTables[self.currentPlate].populateTable( \
                     self.project.dicoPlates[self.currentPlate])
        self.pileResults[self.currentPlate].populateResult( \
                     self.project.dicoPlates[self.currentPlate])

    def modifyAm(self):
        for it in self.pileTables[self.currentPlate].selectedItems():
            am = self.amComboBox.currentText()
            nom = it.statusTip()
            well = getattr(self.project.dicoPlates[self.currentPlate], str(nom))
            well.setAmount(float(am))
        self.project.unsaved = True
        self.fileSaveAction.setEnabled(True)
        self.project.setDicoAm()
        self.projectStack.append(copy.deepcopy(self.project))
        self.pileTables[self.currentPlate].populateTable( \
                    self.project.dicoPlates[self.currentPlate])
        self.pileResults[self.currentPlate].populateResult( \
                    self.project.dicoPlates[self.currentPlate])

    def setType(self):
        for it in self.pileTables[self.currentPlate].selectedItems():
            type = self.typeComboBox.currentText()
            nom = it.statusTip()
            well = getattr(self.project.dicoPlates[self.currentPlate], str(nom))
            well.setType(type)
        self.project.unsaved = True
        self.fileSaveAction.setEnabled(True)
        self.projectStack.append(copy.deepcopy(self.project))
        self.pileTables[self.currentPlate].populateTable( \
                    self.project.dicoPlates[self.currentPlate])
        self.pileResults[self.currentPlate].populateResult( \
                    self.project.dicoPlates[self.currentPlate])

    def enable(self):
        for it in self.pileTables[self.currentPlate].selectedItems():
            ech = self.echComboBox.currentText()
            nom = it.statusTip()
            well = getattr(self.project.dicoPlates[self.currentPlate], str(nom))
            well.setEnabled(True)
        self.project.unsaved = True
        self.fileSaveAction.setEnabled(True)
        self.projectStack.append(copy.deepcopy(self.project))
        self.pileTables[self.currentPlate].populateTable( \
                    self.project.dicoPlates[self.currentPlate])
        self.pileResults[self.currentPlate].populateResult( \
                    self.project.dicoPlates[self.currentPlate])

    def disable(self):
        for it in self.pileTables[self.currentPlate].selectedItems():
            ech = self.echComboBox.currentText()
            nom = it.statusTip()
            well = getattr(self.project.dicoPlates[self.currentPlate], str(nom))
            well.setEnabled(False)
            well.setCtmean('')
            well.setCtdev('')
            well.setNRQ('')
            well.setNRQerror('')
        self.project.unsaved = True
        self.fileSaveAction.setEnabled(True)
        self.projectStack.append(copy.deepcopy(self.project))
        self.pileTables[self.currentPlate].populateTable( \
                        self.project.dicoPlates[self.currentPlate])
        self.pileResults[self.currentPlate].populateResult( \
                        self.project.dicoPlates[self.currentPlate])

    def displayWarnings(self):
        for key in self.project.dicoPlates.keys():
            pl = self.project.dicoPlates[key]
            self.pileTables[key].populateTable(pl)
            self.pileResults[key].populateResult(pl)

    def checkNegative(self, ctMin):
        """
        A method to check negative samples quality
        """
        st = "<b>Warning</b>: ct of the following wells are lower than %.2f :" % \
                ctMin
        st += '<ul>'
        bad = False
        for pl in self.project.dicoPlates:
            for well in self.project.dicoPlates[pl].listePuits:
                if well.type == QString('negative'):
                    if well.ct <= ctMin:
                        st += '<li><b>%s</b> : %.2f </li>' % (well.name, well.ct)
                        bad = True
        st += '</ul>'
        if bad:
            QMessageBox.warning(self, "Warning Negative", st)

    def setRefs(self):
        """
        Determine the reference target and sample
        """
        bad = True
        for plname in self.project.dicoPlates.keys():
            pl = self.project.dicoPlates[plname]
            if pl.contUkn:
                bad = False
                if pl.geneRef == '':
                    QMessageBox.warning(self, "Warning",
                        "Reference target undefined for plate <b>%s</b>!" % plname)
                    raise ValueError
                elif not pl.dicoGene.has_key(pl.geneRef):
                    QMessageBox.warning(self, "Warning",
                        """Wrong reference target for plate <b>%s</b>!
                           This plate does not contain <b>%s</b>.""" % (plname, pl.geneRef))
                    raise ValueError
                if pl.echRef == '':
                    QMessageBox.warning(self, "Warning",
                        "Reference sample undefined for plate <b>%s</b>!" % plname)
                    raise ValueError
                elif not pl.dicoEch.has_key(pl.echRef):
                    QMessageBox.warning(self, "Warning",
                        """Wrong reference sample for plate <b>%s</b>!
                           This plate does not contain <b>%s</b>.""" % (plname, pl.echRef))
                    raise ValueError
        if bad:
            QMessageBox.warning(self, "Warning",
                "The plates does not contain 'unknown'-type wells."
                "The quantification can't be proceeded.")
            raise ValueError

    def computeUnknown(self):
# On verifie que chaque plaque contient des unknown
        self.project.findUnknown()
# On fixe le gene de reference et le triplicat de reference
        try:
            self.setRefs()
        except ValueError:
            return
# On verifie la qualite des negative control
        self.checkNegative(self.ctMin)
# On construit tous les triplicats
        try:
            self.project.findTrip(self.ectMax, self.confidence,
                                      self.errtype)

        except ReplicateError, e:
            st = "<b>Warning</b>: E(ct) of the following replicates is greater"
            st += " than %.2f :" % self.ectMax
            st += str(e)
            QMessageBox.warning(self, "Warning Replicates", st)

        except WellError, e:
            QMessageBox.warning(self, "Problem occurs in calculation !",
                "<b>Warning</b>: A problem occured in the calculations. " \
                "It seems to come from the well(s) <b>%s</b>. " \
                "Check whether ct are correctly defined." \
                % e.brokenWells)  
            self.displayWarnings()
            return

        if self.nplotGene == 0:
            self.onglet.addTab(self.plotUnknownWidget, "Quantification")

# On calcule NRQ
        try:
            self.project.calcNRQ()
        except NRQError, e:
            QMessageBox.warning(self, "Problem occurs in calculation !",
                "<b>Warning</b>: A problem occured in the following replicates : " \
                "%s It probably comes from the reference target or sample which" \
                " is undefined for this gene or sample." \
                "<p> As a consequence these replicates have not been plotted." \
                % str(e))

# On reremplit la table de resultats
        for key in self.project.dicoPlates.keys():
            pl = self.project.dicoPlates[key]
            self.pileResults[key].populateResult(pl)
# On trace le resultat
        self.plotUnknown()
        self.project.unsaved = True
        self.fileSaveAction.setEnabled(True)

    def computeStd(self):
# On cherche les std
        try:
            self.project.findStd(self.ectMax, self.confidence,
                                self.errtype)

        except ReplicateError, e:
            st = "<b>Warning</b>: E(ct) of the following replicates is greater"
            st += " than %.2f" % self.ectMax
            st += "<ul>"
            for trip in e.listRep:
                st += "<li><b>%s, %.2f</b> : E(ct)=%.2f </li>" % (trip.gene, \
                                                  trip.amList[0], trip.ctdev)
            st += "</ul>"
            QMessageBox.warning(self, "Warning Replicates", st)

        except WellError, e:
            QMessageBox.warning(self, "Problem occurs in standard calculation !",
                "<b>Warning</b>: A problem occured in the calculations. " \
                "It seems to come from the well(s) <b>%s</b>. " \
                "Check whether ct are correctly defined." \
                % e.brokenWells)  
            self.displayWarnings()
            self.displayWarnings()
            return

# On trace le resultat on rajoute un onglet si c'est la premiere fois
        if len(self.project.dicoStd.keys()) != 0:
            if self.nplotStd == 0:
                self.onglet.addTab(self.plotStdWidget, "Standard curves")
            self.geneStdBox.clear()
            self.geneStdBox.addItems(self.project.dicoStd.keys())
            # Calcul des courbes standards
            self.project.calcStd(self.confidence, self.errtype)
            self.project.unsaved = True
            self.fileSaveAction.setEnabled(True)
            self.projectStack.append(copy.deepcopy(self.project))
            for key in self.project.dicoPlates.keys():
                pl = self.project.dicoPlates[key]
                self.pileResults[key].populateResult(pl)
            self.populateTree()
            self.plotStd()
        else:
            QMessageBox.warning(self, "Warning",
                "The plates does not contain 'standard'-type wells."
                "The standard curves can't be proceeded.")


    def plotUnknown(self):
        """
        A method to plot the unknown histograms
        """
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

# Gene vs Ech
        valmax = 0
        if self.cboxSens.currentIndex() == 0:
            self.project.findBars(width, spacing, 'geneEch')
            for g in self.project.hashGene.keys()[1:]:
                NRQ = [] ; NRQerror = [] ; valx = []
                for ech in self.project.hashEch.keys()[1:]:
                    for pl in self.project.dicoPlates.keys():
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
                            NRQ, width, color=str(color), 
                            yerr=NRQerror, ecolor='k',
                            label=str(g), align='center')
                    valmax = max(valmax, max(valx))
            self.nplotGene += 1

# Ech vs Gene
        elif self.cboxSens.currentIndex() == 1:
            self.project.findBars(width, spacing, 'echGene')
            for ech in self.project.hashEch.keys()[1:]:
                NRQ = [] ; NRQerror = [] ; valx = []
                for g in self.project.hashGene.keys()[1:]:
                    for pl in self.project.dicoPlates.keys():
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
                            NRQ, width, color=str(color), 
                            yerr=NRQerror, ecolor='k',
                            label=str(ech), align='center')
                    valmax = max(valmax, max(valx))
            self.nplotEch += 1

# plot
        self.mplCanUnknown.axes.set_xticks(self.project.barXticks.values())
        self.mplCanUnknown.axes.set_xticklabels(self.project.barXticks.keys(), 
                                                fontsize=size, 
                                                rotation=int(self.cboxRot.value()))
# Legend + xlim
        self.leg = self.mplCanUnknown.axes.legend(loc='upper right', 
                              shadow=True, labelspacing=0.005)
        legendWidth = 0.3 * valmax
        self.changeFontsize(idraw=False)
        self.mplCanUnknown.axes.set_xlim((0., valmax+legendWidth))
        self.mplCanUnknown.axes.set_ylim(ymin=0.)
        self.mplCanUnknown.draw()

    def plotStd(self):
        """
        A method to plot the standard curves
        """
        self.mplCanStd.axes.cla()
        geneName = self.geneStdBox.currentText()
        x = array([])
        y = array([])
        for trip in self.project.dicoStd[geneName].values():
            x = append(x, asarray(trip.amList))
            y = append(y, asarray(trip.ctList))
        x = log10(x)
        self.mplCanStd.axes.scatter(x, y, marker='o')
        slope, orig = polyfit(x, y, 1)
        yest = polyval([slope, orig], x)
        seps = sqrt(sum((yest-y)**2)/(len(y)-2)) # Formule 2
        sx = sqrt(sum((x-x.mean())**2)/(len(x))) # Formule 3
        stderr = seps / (sx*sqrt(len(x))) # Formule 4 corrigee
        if self.errtype == "student":
            talpha = t.ppf(1.-(1.-self.confidence)/2., len(x)-2) # Student
        elif self.errtype == "normal":
            talpha = norm.ppf(1.-(1.-self.confidence)/2.) # Gaussian
        slopeerr = talpha * stderr
        eff = (10**(-1./slope)-1)*100 # Formule 5 adaptee
        # Erreur(Eff) = (Eff+100) * slopeerr / slope**2 
        stdeff = (eff+100)*log(10)*slopeerr/slope**2 # Formule 6 adaptee
        # Coefficient de Pearsson de correlation
        R2 = 1 - sum((y-yest)**2)/sum((y-mean(y))**2)

        self.mplCanStd.axes.plot(x, yest)

        self.labEquation.setText('ct = %.2f log q0 + %.2f' \
                                % (slope, orig))
        self.labR2.setText('%.3f' % R2)
        self.labEff.setText('%.2f%% %s %.2f' % (eff, unichr(177), stdeff))
        self.mplCanStd.draw()
        self.nplotStd += 1

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
        A method to change the matplotlib xlabels orientation
        """
        size = int(self.cboxRot.value())
        for xtick in self.mplCanUnknown.axes.get_xticklabels():
            xtick.set_rotation(size)
        self.mplCanUnknown.draw()

    def setPlotColor(self):
        dialog = PropDialog(self, hashGene=self.project.hashGene,
                            hashEch=self.project.hashEch)
        if dialog.exec_():
            self.project.hashGene = dialog.hashGene
            self.project.hashEch = dialog.hashEch
            self.plotUnknown()

    def changeCurrentIndex(self):
        """
        A simple method which determines the current plate.
        """
        self.currentIndex = self.tabulPlates.currentIndex()
        self.currentPlate = self.tabulPlates.tabText(self.currentIndex)

    def appendPlate(self, plaque, key):
        mytab = PlateWidget()
        mytab.populateTable(plaque)
        self.tabulPlates.addTab(mytab, key)
        self.pileTables[key] = mytab
        self.connect(mytab, SIGNAL("cellDoubleClicked(int,int)"),
                     self.editWell)

    def appendResult(self, plaque, key):
        mytab = ResultWidget()
        mytab.populateResult(plaque)
        self.tabulResults.addTab(mytab, key)
        self.pileResults[key] = mytab

    def updateUi(self):
        """
        This methods allows to clean up and populate the UI when
        changing some elements.
        """
        self.populateCbox(self.geneComboBox, self.project.hashGene, "Target")
        self.populateCbox(self.echComboBox, self.project.hashEch, "Sample")
        self.populateCbox(self.amComboBox, self.project.hashAmount, "Amount")
        self.populateTree()

    def cleanBeforeOpen(self):
        """
        This methods allows to clean up the UI before a new project.
        """
        # clean-up the QTabWidget
        for ind in range(self.tabulPlates.count()):
            self.tabulPlates.removeTab(0)
            self.tabulResults.removeTab(0)
        self.onglet.removeTab(1)
        self.onglet.removeTab(1)
        # counters=0
        self.nplotGene = 0
        self.nplotStd = 0
        self.nplotEch = 0

        # descativate actions
        self.activateDesactivate(True)
        # undo/redo buffer
        self.projectStack = []


def run():
    import sys
    app = QApplication(sys.argv)
    app.setApplicationName("pyQPCR")
    app.setOrganizationName("pyqpcr")
    app.setOrganizationDomain("pyqpcr.sourceforge.net")
    app.setWindowIcon(QIcon(":/logo.png"))
    f = Qpcr_qt()
    f.show()
    app.exec_()
