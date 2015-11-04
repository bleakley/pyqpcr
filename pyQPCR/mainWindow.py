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
from pyQPCR.widgets.mplPlot import MplUnknownWidget
from pyQPCR.widgets.stdPlot import MplStdWidget
from pyQPCR.widgets.tableMainWindow import PlateWidget, ResultWidget
from pyQPCR.widgets.customCbox import GeneEchComboBox
from pyQPCR.plate import Plaque, ReplicateError, PlateError
from pyQPCR.wellGeneSample import WellError
from pyQPCR.project import Project, NRQError, QabsError, ProjectError
import matplotlib
from numpy import linspace, log10, log, sqrt, sum, mean, polyfit, polyval, \
        asarray, append, array, delete
from scipy.stats import t, norm
import os
import copy
import time

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"
__progversion__ = "0.10dev"

class Qpcr_qt(QMainWindow):
    """
    This class is the build of the main window of pyQPCR. It inheritates
    from QMainWindow.
    """

    def __init__(self, parent=None):
        """
        Construction of the main window. QSettings setup and connexions 
        SIGNAL/SLOT.

        :param parent: the parent QWidget
        :type parent: PyQt4.QtGui.QWidget
        """
        QMainWindow.__init__(self, parent)
 
        self.setMinimumSize(640, 480)
        self.filename = None
        self.printer = None
        self.nplotGene = 0
        self.subIndex = 1

        self.tree = QTreeWidget()
        self.tree.setHeaderLabel("Parameters")

        self.project = Project()
        self.onglet = QTabWidget()
        self.pileTables = OrderedDict()
        self.createProjWidget()
        self.onglet.addTab(self.projWidget, "Plates")
 
        settings = QSettings()
        self.labelRotation, status  = settings.value("Mpl/labelRotation", 
                                                     QVariant(0)).toInt()
        self.labelFontSize, status  = settings.value("Mpl/labelFontSize", 
                                                     QVariant(10)).toInt()
        self.barWth, status  = settings.value("Mpl/barWidth", 
                                              QVariant(0.1)).toDouble()
        self.barSpacing, status = settings.value("Mpl/barSpacing", 
                                                 QVariant(0.3)).toDouble()

        self.topMplCan, status = settings.value("MplCan/top", 
                                                QVariant(0.95)).toDouble()
        self.botMplCan, status = settings.value("MplCan/bottom", 
                                                QVariant(0.10)).toDouble()
        self.rightMplCan, status = settings.value("MplCan/right", 
                                                 QVariant(0.98)).toDouble()
        self.leftMplCan, status = settings.value("MplCan/left", 
                                                 QVariant(0.09)).toDouble()
 
        self.createResultWidget()
        self.pileResults = OrderedDict()
 
        self.vSplitter = QSplitter(Qt.Horizontal)
        self.vSplitter.addWidget(self.tree)
        self.vSplitter.addWidget(self.onglet)
        self.mainSplitter = QSplitter(Qt.Vertical)
        self.mainSplitter.addWidget(self.vSplitter)
        self.mainSplitter.addWidget(self.resulWidget)

        self.setCentralWidget(self.mainSplitter)

        self.vSplitter.setStretchFactor(0, 1)
        self.vSplitter.setStretchFactor(1, 9)
        self.mainSplitter.setStretchFactor(0, 1)
        self.mainSplitter.setStretchFactor(1, 1)

        # Status bar
        status = self.statusBar()
        # Undo / Redo
        self.projectStack = []
        self.undoInd = -1
        # Toolbar and Menus
        self.createMenusAndToolbars()

        # Connexions
        self.connect(self.geneComboBox, SIGNAL("activated(int)"),
                     self.modifyGene)
        self.connect(self.echComboBox, SIGNAL("activated(int)"),
                     self.modifyEch)
        self.connect(self.amComboBox, SIGNAL("activated(int)"),
                     self.modifyAm)
        self.connect(self.typeComboBox, SIGNAL("activated(int)"),
                     self.setType)
        self.connect(self.tabulPlates, SIGNAL("currentChanged(int)"),
                     self.changeCurrentIndex)

        # Settings pour sauvegarde de l'application
        self.recentFiles = settings.value("RecentFiles").toStringList()
        self.ectMax, status = settings.value("EctMax", QVariant(0.3)).toDouble()
        self.ctMin, status = settings.value("ctMin", QVariant(35.)).toDouble()
        self.confidence, status = settings.value("Error/confidence", 
                                                 QVariant(0.9)).toDouble()
        self.errtype = settings.value("Error/errtype", 
                                      QVariant('normal')).toString()
        self.typeCalc = settings.value("Calc/typeCalc", 
                             QVariant('Relative quantification')).toString()
        self.machine = settings.value("machine",
                                      QVariant('Eppendorf Mastercycler')).toString()

        size = settings.value("MainWindow/Size",
                              QVariant(QSize(1024, 768))).toSize()
        self.resize(size)
        position = settings.value("MainWindow/Position",
                                  QVariant(QPoint(0, 0))).toPoint()
        self.move(position)
        self.restoreState(settings.value("MainWindow/State").toByteArray())
        self.vSplitter.restoreState(
                settings.value("MainWindow/VerticalSplitter").toByteArray())
        self.mainSplitter.restoreState(
                settings.value("MainWindow/MainSplitter").toByteArray())
        self.setWindowTitle("pyQPCR")
        self.updateFileMenu()

    def createProjWidget(self):
        """
        This methods create the plate widget.
        """
        self.projWidget = QWidget(self)
        vLay = QVBoxLayout()
        self.tabulPlates = QTabWidget()
        self.tabulPlates.setTabPosition(QTabWidget.West)
        vLay.addWidget(self.tabulPlates)
        self.projWidget.setLayout(vLay)

    def createResultWidget(self):
        """
        This methods create the result widget.
        """
        self.resulWidget = QWidget(self)
        vLay = QVBoxLayout()
        self.tabulResults = QTabWidget()
        self.tabulResults.setTabPosition(QTabWidget.West)
        vLay.addWidget(self.tabulResults)
        self.resulWidget.setLayout(vLay)

    def createMenusAndToolbars(self):
        """
        This methods create the menus and toolbars of pyQPCR.
        """
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
        self.applyModelAction = self.createAction("Apply &model...",
                self.applyModel, icon="applymodel",
                tip="Apply a model")
        fileQuitAction = self.createAction("&Quit", self.close, 
                "Ctrl+Q", "filequit", "Close the application")
        self.editAction = self.createAction("Edit wells", self.editWell, 
                "Ctrl+E", "edit", "Edit selected wells")
        self.undoAction = self.createAction("Undo", self.undo, 
                QKeySequence.Undo, "undo", "Undo")
        self.undoAction.setEnabled(False)
        self.redoAction = self.createAction("Redo", self.redo,
                QKeySequence.Redo, "redo", "Redo")
        self.redoAction.setEnabled(False)
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
        self.extractAction = self.createAction("E&xtract sub-plates",
                             self.extractSubplate, "Ctrl+X",
                           "extract", "Extract a sub-plate from an existing plate")
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
                        self.fileSaveAction, self.fileSaveAsAction, self.applyModelAction,
                        None, fileQuitAction))

        editMenu = self.menuBar().addMenu("&Edit")
        editMenu.addAction(self.editAction)
        editMenu.addSeparator()
        self.addActions(editMenu, (self.undoAction, self.redoAction))
        editMenu.addSeparator()
        self.addActions(editMenu, (self.addEchAction, self.addGeneAction,
                                   self.addAmAction))
        editMenu.addSeparator()
        editMenu.addAction(self.extractAction)

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
                        self.fileSaveAction, self.fileSaveAsAction,
                        self.applyModelAction))
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
        editToolbar.addAction(self.extractAction)
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
                                      self.applyModelAction,
                                      self.plotStdAction, self.plotAction))
        
        # Desactivation par defaut
        self.activateDesactivate(False)

    def createAction(self, text, slot=None, shortcut=None, icon=None,
                     tip=None, checkable=None, signal="triggered()"):
        """
        This method highly simplifies the creation of QAction.
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
        """
        A simple method to display a message in the status bar

        :param message: the message to be displayed
        :type message: PyQt4.QtCore.QString
        :param time: the time the message is displayed (in seconds)
        :type time: int
        """
        self.statusBar().showMessage(message, time)

    def fileOpen(self):
        if not self.okToContinue():
            return
        dir = os.path.dirname(self.filename) if self.filename is not None \
                else "."
        formats =[u"*.xml"]
        fname = unicode(QFileDialog.getOpenFileName(self,
                                                    "pyQPCR - Choose a file", 
                dir, "pyQPCR files (%s)" % " ".join(formats)))
        if fname:
            self.loadFile(fname)

    def loadFile(self, fname=None):
        """
        A method that loads an XML file

        :param fname: the name of the file
        :type fname: PyQt4.QtCore.QString
        """
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

            try:
                self.project = Project(fname)
            except ProjectError, e:
                QMessageBox.warning(self, "Problem in import !", "%s" % str(e))
                return

            for key in self.project.dicoPlates.keys():
                pl = self.project.dicoPlates[key]
                self.appendPlate(pl, key)
                self.appendResult(pl, key)

            # Pile de plaques pour le Undo/Redo
            self.projectStack.insert(len(self.projectStack) + self.undoInd + 1, 
                                     copy.deepcopy(self.project))
            if self.undoInd != -1:
                del self.projectStack[self.undoInd+1:]
                self.undoInd = -1
                self.redoAction.setEnabled(False)
            self.updateUi()

    def fileNew(self):
        """
        This method is used when one wants to create a new project.
        It opens the NewProjectDialog, a wizard that helps to define your
        new project.
        """
        if not self.okToContinue():
            return
        dir = os.path.dirname(self.filename) if self.filename is not None \
                else "."
        dialog = NewProjectDialog(self, pwd=dir, machine=self.machine)
        if dialog.exec_():
            self.cleanBeforeOpen()
            self.project = Project(dialog.projectName, open=False)
            self.filename = dialog.projectName
            self.setWindowTitle("pyQPCR - %s[*]" % QFileInfo(self.filename).fileName())
            for fname in dialog.fileNames.values():
                self.machine = dialog.machineType
                self.typeCalc = dialog.calculationType
                self.addPlate(fname, fileNew=True)

            self.projectStack.insert(len(self.projectStack) + self.undoInd + 1,
                                     copy.deepcopy(self.project))
            if self.undoInd != -1:
                del self.projectStack[self.undoInd+1:]
                self.undoInd = -1
                self.redoAction.setEnabled(False)

    def fileImport(self):
        """
        This method allows to import raw data file and XML pyQPCR file format.
        It lets the user choose the files corresponding to the device chosen in
        the preference dialog.
        """
        dir = os.path.dirname(self.filename) if self.filename is not None \
                else "."
        if self.machine == 'Applied StepOne':
            formats =[u"*.txt", u"*.csv"]
            type = 'Applied StepOne machines'
        elif self.machine == 'Applied 7000':
            formats =[u"*.txt", u"*.csv"]
            type = 'Applied 7000 machines'
        elif self.machine == 'Applied 7500':
            formats =[u"*.txt", u"*.csv"]
            type = 'Applied 7500 machines'
        elif self.machine == 'Applied 7700':
            formats =[u"*.csv"]
            type = 'Applied 7700 machines'
        elif self.machine == 'Applied 7900':
            formats =[u"*.txt"]
            type = 'Applied 7900 machines'
        if self.machine == 'Applied Viia7':
            formats =[u"*.txt", u"*.csv"]
            type = 'Applied Viia7 machines'
        elif self.machine == 'Biorad C1000':
            formats =[u"*.txt"]
            type = 'Biorad C1000 machines'
        elif self.machine == 'Biorad MyIQ':
            formats =[u"*.csv"]
            type = 'Biorad MyIQ machines'
        elif self.machine == 'Biorad Opticon':
            formats =[u"*.txt"]
            type = 'Biorad Opticon machines'
        elif self.machine == 'Cepheid SmartCycler':
            formats =[u"*.csv"]
            type = 'Cepheid SmartCycler machines'
        if self.machine == 'Eppendorf Mastercycler':
            formats =[u"*.txt", u"*.csv"]
            type = 'Eppendorf Mastercycler machines'
        elif self.machine == 'Esco Spectrum 48':
            formats =[u"*.csv"]
            type = 'Esco Spectrum 48 machines'
        elif self.machine == 'Illumina Eco':
            formats =[u"*.csv"]
            type = 'Illumina Eco machines'
        elif self.machine == 'Qiagen Corbett':
            formats =[u"*.csv"]
            type = 'Qiagen Corbett machines'
        elif self.machine == 'Roche LightCycler 480':
            formats =[u"*.txt"]
            type = 'Roche LightCycler 480'
        elif self.machine == 'Stratagene Mx3000':
            formats =[u"*.txt"]
            type = 'Stratagene Mx3000'
        else:
            formats =[u"*.txt", u"*.csv"]
            type = 'Eppendorf Mastercycler machines'
        fileNames = QFileDialog.getOpenFileNames(self,
                       "pyQPCR - Choose a file", dir, 
                       "Input files [%s] (%s);;pyQPCR file (*.xml)" % (type, 
                                " ".join(formats)))
        if fileNames:
            for file in fileNames:
                if file.endsWith('xml') or file.endsWith('XML'):
                    st = '<ul>'
                    warn = False
                    try:
                        prtmp = Project(file, initIndexColor=self.project.indexColor)
                    except ProjectError, e:
                        QMessageBox.warning(self, "Problem in import !", "%s" % str(e))
                        return
                    for plate in prtmp.dicoPlates.keys():
                        if not self.project.dicoPlates.has_key(plate):
                            self.project.addPlate(prtmp.dicoPlates[plate], key=plate)
                            self.appendPlate(prtmp.dicoPlates[plate], plate)
                            self.appendResult(prtmp.dicoPlates[plate], plate)
                        else:
                            st += '<li><b>%s</b></li>\n' % plate
                            warn = True
                    st += '</ul>'
                    if warn:
                        QMessageBox.warning(self, "Import not complete",
                         "<b>Warning</b>: the plate(s): %s is (are) already in the project ! \
                         As a consequence, it (they) has (have) not been added."       
                                        % (st))

                    self.updateUi()
                    self.project.unsaved = True
                    self.fileSaveAction.setEnabled(True)
                    self.undoAction.setEnabled(True)
                    self.projectStack.insert(len(self.projectStack)+self.undoInd+1,
                                             copy.deepcopy(self.project))
                    if self.undoInd != -1:
                        del self.projectStack[self.undoInd+1:]
                        self.undoInd = -1
                        self.redoAction.setEnabled(False)
                else:
                    if not self.project.dicoPlates.has_key(QFileInfo(file).fileName()):
                        self.addPlate(file)
                    else:
                        name = QFileInfo(file).fileName()
                        QMessageBox.warning(self, "Import cancelled",
                          "<b>Warning</b>: the plate %s is already in the project %s ! \
                           As a consequence, it has not been added to the project."       
                                        % (name, QFileInfo(self.filename).fileName()))

    def addPlate(self, fname=None, fileNew=False):
        """
        This method is used to add a plate to the project.

        :param fname: the name of the plate
        :type fname: PyQt4.QtCore.QString
        :param fileNew: a boolean to indicate wheter it is a new project
                        or an additional plate to an existing project
        :type fileNew: logical
        """
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
                if not fileNew:
                    self.undoAction.setEnabled(True)
                    self.projectStack.insert(len(self.projectStack) + self.undoInd + 1,
                                             copy.deepcopy(self.project))
                    if self.undoInd != -1:
                        del self.projectStack[self.undoInd+1:]
                        self.undoInd = -1
                        self.redoAction.setEnabled(False)
            except KeyError:
                st = "<b>Warning:</b>: an error occurs during import ! "
                st += "It probably comes from your file which has not been correctly parsed. "
                st += "Is it a raw <b>%s</b> file ?" % self.machine
                QMessageBox.warning(self, "Warning file import", st)
            except PlateError, e:
                QMessageBox.warning(self, "Warning file import", str(e))

    def closePlate(self):
        """
        This method is used to close a plate from your project. A dialog
        will warn you that you are going to close a plate.
        """
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
            self.undoAction.setEnabled(True)
            self.projectStack.insert(len(self.projectStack) + self.undoInd + 1,
                                     copy.deepcopy(self.project))
            if self.undoInd != -1:
                del self.projectStack[self.undoInd+1:]
                self.undoInd = -1
                self.redoAction.setEnabled(False)

    def activateDesactivate(self, bool):
        """
        This method allows to enable/disable several QAction

        :param bool: the boolean value
        :type bool: bool
        """
        self.addGeneAction.setEnabled(bool)
        self.addEchAction.setEnabled(bool)
        self.addAmAction.setEnabled(bool)
        self.editAction.setEnabled(bool)
        self.plotAction.setEnabled(bool)
        self.plotStdAction.setEnabled(bool)
        self.extractAction.setEnabled(bool)
        self.typeComboBox.setEnabled(bool)
        self.geneComboBox.setEnabled(bool)
        self.echComboBox.setEnabled(bool)
        self.amComboBox.setEnabled(bool)
        self.fileSaveAsAction.setEnabled(bool)
        self.applyModelAction.setEnabled(bool)
        self.filePrintAction.setEnabled(bool)
        self.exportAction.setEnabled(bool)
        self.enableAction.setEnabled(bool)
        self.disableAction.setEnabled(bool)
        self.fileImportAction.setEnabled(bool)
        self.closeTabAction.setEnabled(bool)

    def populateTree(self):
        """
        This method allows to populate the QTreeWidget which displays informations
        about the current project (reference target and samples, targets efficiencies, ...)
        """
        self.tree.clear()
        if self.filename is not None:
            ancestor = QTreeWidgetItem(self.tree, 
                                       [QFileInfo(self.filename).fileName()])
        else:
            ancestor = QTreeWidgetItem(self.tree, ["untitled project"])
        itemErr = QTreeWidgetItem(ancestor, ["Error calculation"])
        st = 'Error type : %s' % self.errtype
        item = QTreeWidgetItem(itemErr, [st])
        st = 'Confidence level : %s' % self.confidence
        item = QTreeWidgetItem(itemErr, [st])
        for key in self.project.dicoPlates.keys():
            item = QTreeWidgetItem(ancestor, [key])
            pl = self.project.dicoPlates[key]
            itemQuant = QTreeWidgetItem(item, ["Quantification"])
            itemRefGene = QTreeWidgetItem(itemQuant , ["Reference Target"])
            itemRefEch = QTreeWidgetItem(itemQuant , ["Reference Sample"])
            for geneName in pl.geneRef:
                item = QTreeWidgetItem(itemRefGene, [geneName])
            item = QTreeWidgetItem(itemRefEch, [pl.echRef])
        itemStd = QTreeWidgetItem(ancestor, ["Standard"])
        for gene in self.project.hashGene.values()[1:]:
            try:
                eff = "%s:%.2f%s%.2f" % (gene.name, gene.eff, unichr(177), gene.pm)
            except TypeError:
                eff = '%s:' % gene.name
            item = QTreeWidgetItem(itemStd, [eff])
        self.tree.expandAll()

    def populateCbox(self, cbox, items, name="Target"):
        """
        A simple method to populate a QComboBox object.

        :param cbox: the comboBox we want to populate
        :type cbox: PyQt4.QtGui.QComboBox
        :param items: the items we want to add
        :type items: list
        :param name: the header of the QComboBox
        :type name: string
        """
        cbox.clear()
        cbox.addItem(name, QVariant("header"))
        cbox.addItems(items)

    def fileSave(self):
        """
        Save the project in a XML format.
        """
        self.project.exportXml(self.filename)
        self.updateStatus("Saved %s" % self.filename)
        self.project.unsaved = False
        self.fileSaveAction.setEnabled(False)

    def fileSaveAs(self):
        """
        Save as wizard.
        """
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
        """
        This method is used to generate an HTML output that corresponds to the
        result table. It is used to generate the PDF output of the result table.
        """
        dialog = PrintingDialog(self)
        if dialog.exec_():
            isPlate = dialog.btnPlate.isChecked()
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
               "h1 {background-color:#9bcd9b; color:black;}"
               "h2 {background-color:SlateGrey; color:white;}"
               "</style>\n"
               "</head>\n")
        html += css
        html += "<h1 align=center> qPCR results </h1><br><br>"
        html += "<br><h2>Setup</h2><br>\n"
        html += "<ul>\n"
        html +=  "<li> <b>Calculation type :</b> %s</li>\n" % self.typeCalc
        if self.errtype == "student":
            html +=  "<li> <b>Error type :</b> Student t-test </li>\n" 
        elif self.errtype == "normal":
            html +=  "<li> <b>Error type :</b> Gaussian </li>\n" 
        html +=  "<li> <b>Confidence interval :</b> %.2f &#37;</li>\n" % (100*self.confidence)
        html +=  "<li> <b>Machine :</b> %s</li>\n" % self.machine
        html +=  "<li> <b>Maximum E(Ct) :</b> %.2f</li>\n" % self.ectMax
        html +=  "<li> <b>Minimum Ct :</b> %.2f</li>\n" % self.ctMin
        html +=  "<li> <b>Date :</b> %s</li>\n" % time.strftime('%d/%m/%Y',time.localtime())
        html +=  "<li> <b>pyQPCR version :</b> %s</li>\n" % __progversion__
        html +=  "</ul>\n"

        if isPlate:
            for index, key in enumerate(self.project.dicoPlates.keys()):
                if index == 1:
                    html += "<p style='page-break-before:always;'>"
                    html += "<br><h2>Plate setup (%s)</h2><br>" % key
                    html += self.project.dicoPlates[key].writeHtmlPlateMap()
                    html += "</p>"
                else:
                    html += "<br><h2>Plate setup (%s)</h2><br>" % key
                    html += "<p>"
                    html += self.project.dicoPlates[key].writeHtmlPlateMap()
                    html += "</p>"

        if isTable:
            for key in self.project.dicoPlates.keys():
                html += "<br><h2>Results table (%s)</h2><br>" % key
                html += self.project.dicoPlates[key].writeHtml(self.ctMin, self.ectMax,
                                                               self.typeCalc)

        if isStd and self.nplotStd !=0:
            html += "<p style='page-break-before:always;'>"
            html += "<br><h2>Standard curves</h2><br>"
            #self.plotStdWidget.geneStdBox.addItems(self.project.dicoStd.keys())
            #for index in range(len(self.project.dicoStd.keys())):
            for index, data in enumerate(self.project.dicoPlotStd.keys()):
                if (not index % 3 and index != 0):
                    html += "<p style='page-break-before:always;'>"
                else:
                    html += "<p>"
                html += "<table border=0 width=100%>\n"
                #self.plotStdWidget.geneStdBox.setCurrentIndex(index)
                self.plotStdWidget.plotStd(nameGene=data)
                fig = self.plotStdWidget.mplCanStd.figure.savefig("output%i.png" % index, 
                                                    dpi=250)
                html += "<tr valign=middle>\n"
                html += ("<th align=center>"
                         "<table width=100% border=0>")
                html += "<tr><th><font size 10pt><b>Gene : </b> %s</th></tr>" %\
                        data
                html += "<tr><th><b>Linear Regression : </b> %s</th></tr>" %\
                        self.plotStdWidget.labEquation.text() 
                html += "<tr><th><b>Efficiency : </b> %s</th></tr>" % \
                        self.plotStdWidget.labEff.text()
                html += "<tr><th><b>R&#178; : </b> %s</th></tr>" % \
                        self.plotStdWidget.labR2.text()
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

        if isQuant and self.nplotGene != 0:
            html += "<br><h2>Quantification curves</h2>"
            fig = self.mplUknWidget.mplCanUnknown.figure.savefig("output.png", dpi=500)
            html += "<p align=center><img src='output.png' width=500></p>"
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
            self.printer.setColorMode(QPrinter.Color)
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
        """
        This method is called when the user want to configure pyQPCR settings.
        It opens a dialog in which the user can change the preferences.
        """
        dialog = SettingsDialog(self, ect=self.ectMax,
                                ctmin=self.ctMin,
                                confidence=self.confidence,
                                errtype=self.errtype,
                                machine=self.machine,
                                typeCalc=self.typeCalc)
        if dialog.exec_():
            self.ectMax, st = dialog.ectLineEdit.text().toDouble()
            self.ctMin, st = dialog.ctMinLineEdit.text().toDouble()
            self.confidence, st = dialog.confCbox.currentText().toDouble()
            errtype = dialog.typeCbox.currentText()
            self.machine = dialog.machBox.currentText()
            self.typeCalc = dialog.typeCalc.currentText()
            self.errtype = dialog.types[errtype]
            self.confidence /= 100
            self.populateTree()

    def helpAbout(self):
        """
        This method is called when the user clicks on "About pyQPCR" in the 
        Help menu. 
        """
        import platform
        QMessageBox.about(self, "About pyQPCR",
        """<b>pyQPCR</b> v %s
        <p>Copyright &copy; 2008-2012 Thomas Gastine - Magali Hennion - 
        All rights reserved.
        <p>This application can be used to perform
        simple qPCR analysis.
        <p> It may be used, copied and modified with no restriction
        <p>Python %s - PyQt %s - Matplotlib %s
        on %s""" % (__progversion__, platform.python_version(),
                    PYQT_VERSION_STR, matplotlib.__version__, platform.system()))

    def helpHelp(self):
        """
        A method which display the help dialog.
        """
        f = HelpDialog("index.html", self)
        f.show()

    def addRecentFile(self, fname):
        """
        A simple method to add the file name to the 'recent files' menu entry.

        :param fname: the file name
        :type fname: string
        """
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
        """
        This action is called when one wants to close pyQPCR.
        It ensures that all the preferences of the software are saved
        (QSettings and QVariant stuff).
        """
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
            typeCalc = QVariant(self.typeCalc) if self.typeCalc \
                  else QVariant()
            machine = QVariant(self.machine) if self.machine \
                  else QVariant()

            if self.nplotGene != 0:
                top = QVariant(self.mplUknWidget.mplCanUnknown.fig.subplotpars.top)
                bottom = QVariant(self.mplUknWidget.mplCanUnknown.fig.subplotpars.bottom)
                right = QVariant(self.mplUknWidget.mplCanUnknown.fig.subplotpars.right)
                left = QVariant(self.mplUknWidget.mplCanUnknown.fig.subplotpars.left)
                rot = QVariant(int(self.mplUknWidget.cboxRot.value()))
                barW = QVariant(self.mplUknWidget.spinWidth.value())
                barS = QVariant(self.mplUknWidget.spinSpacing.value())
                fontSize = QVariant(int(self.mplUknWidget.cboxFontsize.value()))
                settings.setValue("MplCan/top", top)
                settings.setValue("MplCan/bottom", bottom)
                settings.setValue("MplCan/right", right)
                settings.setValue("MplCan/left", left)
                settings.setValue("Mpl/labelRotation", rot)
                settings.setValue("Mpl/labelFontSize", fontSize)
                settings.setValue("Mpl/barWidth", barW)
                settings.setValue("Mpl/barSpacing", barS)

            settings.setValue("Error/confidence", confidence)
            settings.setValue("Error/errtype", errtype)
            settings.setValue("Calc/typeCalc", typeCalc)
            settings.setValue("machine", machine)
            settings.setValue("MainWindow/Size", QVariant(self.size()))
            settings.setValue("MainWindow/Position",
                    QVariant(self.pos()))
            settings.setValue("MainWindow/State", QVariant(self.saveState()))
            settings.setValue("MainWindow/VerticalSplitter", 
                    QVariant(self.vSplitter.saveState()))
            settings.setValue("MainWindow/MainSplitter", 
                    QVariant(self.mainSplitter.saveState()))
        else:
            event.ignore()

    def okToContinue(self):
        """
        This small dialog is called when the user wants to exit but has
        unsaved stuff. So it asks to confirm that the user really wants
        to exit without saving.
        """
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
        """
        Called when 'redo' is clicked.
        """
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
            self.pileResults[key].populateResult(pl, self.typeCalc)

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
            self.undoAction.setEnabled(False)
        else:
            self.project.unsaved = True
            self.fileSaveAction.setEnabled(True)
            self.undoAction.setEnabled(True)
        if self.undoInd == -1:
            self.redoAction.setEnabled(False)

    def undo(self):
        """
        Called when 'undo' is clicked.
        """
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
            self.pileResults[key].populateResult(pl, self.typeCalc)

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
            self.undoAction.setEnabled(False)
        else:
            self.project.unsaved = True
            self.fileSaveAction.setEnabled(True)
            self.undoAction.setEnabled(True)
        self.redoAction.setEnabled(True)

    def applyModel(self):
        dir = os.path.dirname(self.filename) if self.filename is not None \
                else "."
        dialog = ModelDialog(self, pr=self.project, pwd=dir)
        hasChanged = False
        if dialog.exec_():
            if dialog.project != self.project:
                self.project = dialog.project
                hasChanged = True
        if hasChanged:
            self.project.unsaved = True
            self.fileSaveAction.setEnabled(True)
            self.undoAction.setEnabled(True)
            self.projectStack.insert(len(self.projectStack) + self.undoInd + 1,
                                     copy.deepcopy(self.project))
            if self.undoInd != -1:
                del self.projectStack[self.undoInd+1:]
                self.undoInd = -1
                self.redoAction.setEnabled(False)
            for key in self.project.dicoPlates.keys():
                pl = self.project.dicoPlates[key]
                self.pileTables[key].populateTable(pl)
                self.pileResults[key].populateResult(pl, self.typeCalc)
            self.updateUi()

    def editWell(self):
        """
        This method is called when the user want to edit the properties of the
        selected wells.
        """
        setType = set()
        setGene = set()
        setEch = set()
        setAm = set()
        selected = [setType, setGene, setEch, setAm]
        selectedWells = []
        for it in self.pileTables[self.currentPlate].selectedItems():
            nom = it.statusTip()
            well = getattr(self.project.dicoPlates[self.currentPlate], str(nom))
            setType.add(well.type)
            setEch.add(well.ech.name)
            setGene.add(well.gene.name)
            setAm.add(well.amount)
            selectedWells.append(well)
        if len(selectedWells) > 0:
            dialog = EditDialog(self, project=self.project, selected=selected)
             
            hasChanged = False
            if dialog.exec_():
                ge = dialog.cboxGene.currentObj()
                if dialog.project != self.project:
                    self.project = dialog.project
                    hasChanged = True
                if dialog.cboxType.currentText() == QString('unknown'):
                    ech = dialog.cboxSample.currentObj()
                    for it in self.pileTables[self.currentPlate].selectedItems():
                        nom = it.statusTip()
                        well = getattr(self.project.dicoPlates[self.currentPlate],
                                       str(nom))
                        if well.type != 'unknown':
                            well.setType(QString('unknown'))
                            hasChanged = True
                        if ge.name != '' and well.gene != ge:
                            well.setGene(ge)
                            hasChanged = True
                        if ech.name != '' and well.ech != ech:
                            well.setEch(ech)
                            hasChanged = True
                    if hasChanged:
                        for pl in self.project.dicoPlates.values():
                            pl.setDicoEch()
                if dialog.cboxType.currentText() == QString('standard'):
                    am = dialog.cboxAm.currentText()
                    for it in self.pileTables[self.currentPlate].selectedItems():
                        nom = it.statusTip()
                        well = getattr(self.project.dicoPlates[self.currentPlate],
                                       str(nom))
                        if well.type != 'standard':
                            well.setType(QString('standard'))
                            hasChanged = True
                        if ge.name != '' and well.gene != ge:
                            well.setGene(ge)
                            hasChanged = True
                        if am != '' and well.amount != float(am):
                            well.setAmount(float(am))
                            hasChanged = True
                    if hasChanged:
                        self.project.setDicoAm()
                if hasChanged:
                    self.project.unsaved = True
                    self.fileSaveAction.setEnabled(True)
                    self.undoAction.setEnabled(True)
                    for pl in self.project.dicoPlates.values():
                        pl.setDicoGene()
                    self.projectStack.insert(len(self.projectStack) + self.undoInd + 1,
                                             copy.deepcopy(self.project))
                    if self.undoInd != -1:
                        del self.projectStack[self.undoInd+1:]
                        self.undoInd = -1
                        self.redoAction.setEnabled(False)
                    for key in self.project.dicoPlates.keys():
                        pl = self.project.dicoPlates[key]
                        self.pileTables[key].populateTable(pl)
                        self.pileResults[key].populateResult(pl, self.typeCalc)
                    self.updateUi()

        else: #Please select something !!
            QMessageBox.information(self, "Please select some wells",
                "You must select at least one well to edit something.")

    def addGene(self):
        """
        This method is called when the user want to add/edit the targets of
        the project. It opens a dialog which let the user add/remove/edit the
        properties of the targets.

        >>> dialog = GeneDialog(self, project=self.project)
        """
        dialog = GeneDialog(self, project=self.project)
        if dialog.exec_():
            project = dialog.project
            if project != self.project:
                self.populateCbox(self.geneComboBox, project.hashGene, "Target")
                self.project = project
                self.fileSaveAction.setEnabled(True)
                self.undoAction.setEnabled(True)
                for pl in self.project.dicoPlates.values():
                    pl.setDicoGene()
                self.projectStack.insert(len(self.projectStack) + self.undoInd + 1,
                                         copy.deepcopy(self.project))
                if self.undoInd != -1:
                    del self.projectStack[self.undoInd+1:]
                    self.undoInd = -1
                    self.redoAction.setEnabled(False)
                for key in self.project.dicoPlates.keys():
                    pl = self.project.dicoPlates[key]
                    self.pileTables[key].populateTable(pl)
                    self.pileResults[key].populateResult(pl, self.typeCalc)
                self.populateTree()

    def addEch(self):
        """
        This method is called when the user want to add/edit the samples of
        the project. It opens a dialog which let the user add/remove/edit the
        properties of the samples.

        >>> dialog = EchDialog(self, project=self.project)
        """
        dialog = EchDialog(self, project=self.project)
        if dialog.exec_():
            project = dialog.project
            if project != self.project:
                self.populateCbox(self.echComboBox, project.hashEch, "Sample")
                self.project = project
                self.fileSaveAction.setEnabled(True)
                self.undoAction.setEnabled(True)
                for pl in self.project.dicoPlates.values():
                    pl.setDicoEch()
                self.projectStack.insert(len(self.projectStack) + self.undoInd + 1,
                                         copy.deepcopy(self.project))
                if self.undoInd != -1:
                    del self.projectStack[self.undoInd+1:]
                    self.undoInd = -1
                    self.redoAction.setEnabled(False)
                for key in self.project.dicoPlates.keys():
                    pl = self.project.dicoPlates[key]
                    self.pileTables[key].populateTable(pl)
                    self.pileResults[key].populateResult(pl, self.typeCalc)
                self.populateTree()

    def addAmount(self):
        """
        This method is called when the user want to add/edit the amounts of
        the project. It opens a dialog which let the user add/remove/edit the
        values of the amounts.

        >>> dialog = AmountDialog(self, project=self.project)
        """
        dialog = AmountDialog(self, project=self.project)
        if dialog.exec_():
            project = dialog.project
            if project != self.project:
                self.populateCbox(self.amComboBox, project.hashAmount, "Amount")
                self.project = project
                self.fileSaveAction.setEnabled(True)
                self.undoAction.setEnabled(True)
                self.project.setDicoAm()
                self.projectStack.insert(len(self.projectStack) + self.undoInd + 1,
                                         copy.deepcopy(self.project))
                if self.undoInd != -1:
                    del self.projectStack[self.undoInd+1:]
                    self.undoInd = -1
                    self.redoAction.setEnabled(False)
                for key in self.project.dicoPlates.keys():
                    pl = self.project.dicoPlates[key]
                    self.pileTables[key].populateTable(pl)
                    self.pileResults[key].populateResult(pl, self.typeCalc)
                self.populateTree()

    def extractSubplate(self):
        """
        This method is used to extract a sub-plate from an opened plate.
        """
        listWell = []
        for it in self.pileTables[self.currentPlate].selectedItems():
            nom = it.statusTip()
            well = getattr(self.project.dicoPlates[self.currentPlate], str(nom))
            listWell.append(well)

        if len(listWell) > 0:
            fname = self.currentPlate + '_sub%i' % self.subIndex
            plaque = Plaque(fname, machine=None)
            plaque.setPlateType(self.project.dicoPlates[self.currentPlate].type)
            plaque.subPlate(listWell)

            if not self.project.dicoPlates.has_key(fname):
                self.project.addPlate(plaque)
                key = QFileInfo(fname).fileName()

                self.appendPlate(plaque, key)
                self.appendResult(plaque, key)

                self.updateUi()
                self.project.unsaved = True
                self.fileSaveAction.setEnabled(True)
                self.undoAction.setEnabled(True)
                self.projectStack.insert(len(self.projectStack) + self.undoInd + 1,
                                         copy.deepcopy(self.project))
                if self.undoInd != -1:
                    del self.projectStack[self.undoInd+1:]
                    self.undoInd = -1
                    self.redoAction.setEnabled(False)

            self.subIndex += 1
        else: #Please select something !!
            QMessageBox.information(self, "Please select some wells",
                "You must select at least one well to create a new subplate.")

    def modifyGene(self):
        """
        This method is called when the user changes the gene of the selected wells.
        """
        hasChanged = False
        for it in self.pileTables[self.currentPlate].selectedItems():
            gene = self.geneComboBox.currentObj()
            nom = it.statusTip()
            well = getattr(self.project.dicoPlates[self.currentPlate], str(nom))
            if well.gene != gene:
                well.setGene(gene)
                hasChanged = True
        if hasChanged:
            self.project.unsaved = True
            self.fileSaveAction.setEnabled(True)
            self.undoAction.setEnabled(True)
            self.project.dicoPlates[self.currentPlate].setDicoGene()
            self.projectStack.insert(len(self.projectStack) + self.undoInd + 1,
                                     copy.deepcopy(self.project))
            if self.undoInd != -1:
                del self.projectStack[self.undoInd+1:]
                self.undoInd = -1
                self.redoAction.setEnabled(False)
            self.pileTables[self.currentPlate].populateTable( \
                          self.project.dicoPlates[self.currentPlate])
            self.pileResults[self.currentPlate].populateResult( \
                          self.project.dicoPlates[self.currentPlate], 
                          self.typeCalc)

    def modifyEch(self):
        """
        This method is called when the user changes the sample of the selected 
        wells.
        """
        hasChanged = False
        for it in self.pileTables[self.currentPlate].selectedItems():
            ech = self.echComboBox.currentObj()
            nom = it.statusTip()
            well = getattr(self.project.dicoPlates[self.currentPlate], str(nom))
            if well.ech != ech:
                well.setEch(ech)
                hasChanged = True
        if hasChanged:
            self.project.unsaved = True
            self.fileSaveAction.setEnabled(True)
            self.undoAction.setEnabled(True)
            self.project.dicoPlates[self.currentPlate].setDicoEch()
            self.projectStack.insert(len(self.projectStack) + self.undoInd + 1,
                                     copy.deepcopy(self.project))
            if self.undoInd != -1:
                del self.projectStack[self.undoInd+1:]
                self.undoInd = -1
                self.redoAction.setEnabled(False)
            self.pileTables[self.currentPlate].populateTable( \
                         self.project.dicoPlates[self.currentPlate])
            self.pileResults[self.currentPlate].populateResult( \
                         self.project.dicoPlates[self.currentPlate], self.typeCalc)

    def modifyAm(self):
        """
        This method is called when the user changes the amount of the selected 
        wells.
        """
        hasChanged = False
        for it in self.pileTables[self.currentPlate].selectedItems():
            am = self.amComboBox.currentText()
            nom = it.statusTip()
            well = getattr(self.project.dicoPlates[self.currentPlate], str(nom))
            if well.amount != float(am):
                well.setAmount(float(am))
                hasChanged = True
        if hasChanged:
            self.project.unsaved = True
            self.fileSaveAction.setEnabled(True)
            self.undoAction.setEnabled(True)
            self.project.setDicoAm()
            self.projectStack.insert(len(self.projectStack) + self.undoInd + 1,
                                     copy.deepcopy(self.project))
            if self.undoInd != -1:
                del self.projectStack[self.undoInd+1:]
                self.undoInd = -1
                self.redoAction.setEnabled(False)
            self.pileTables[self.currentPlate].populateTable( \
                        self.project.dicoPlates[self.currentPlate])
            self.pileResults[self.currentPlate].populateResult( \
                        self.project.dicoPlates[self.currentPlate], self.typeCalc)

    def setType(self):
        """
        This method is called when the user changes the type of the selected 
        wells.
        """
        hasChanged = False
        for it in self.pileTables[self.currentPlate].selectedItems():
            type = self.typeComboBox.currentText()
            nom = it.statusTip()
            well = getattr(self.project.dicoPlates[self.currentPlate], str(nom))
            if well.type != type:
                well.setType(type)
                hasChanged = True
        if hasChanged:
            self.project.unsaved = True
            self.fileSaveAction.setEnabled(True)
            self.undoAction.setEnabled(True)
            self.projectStack.insert(len(self.projectStack)+self.undoInd+1,
                                     copy.deepcopy(self.project))
            if self.undoInd != -1:
                del self.projectStack[self.undoInd+1:]
                self.undoInd = -1
                self.redoAction.setEnabled(False)
            self.pileTables[self.currentPlate].populateTable( \
                        self.project.dicoPlates[self.currentPlate])
            self.pileResults[self.currentPlate].populateResult( \
                        self.project.dicoPlates[self.currentPlate], self.typeCalc)

    def enable(self):
        """
        This method is called when the user want to enable the selected wells
        from the computation.
        """
        hasChanged = False
        for it in self.pileTables[self.currentPlate].selectedItems():
            ech = self.echComboBox.currentText()
            nom = it.statusTip()
            well = getattr(self.project.dicoPlates[self.currentPlate], str(nom))
            if well.enabled is False:
                well.setEnabled(True)
                hasChanged = True
        if hasChanged:
            self.project.unsaved = True
            self.fileSaveAction.setEnabled(True)
            self.undoAction.setEnabled(True)
            self.projectStack.insert(len(self.projectStack) + self.undoInd + 1,
                                     copy.deepcopy(self.project))
            if self.undoInd != -1:
                del self.projectStack[self.undoInd+1:]
                self.undoInd = -1
                self.redoAction.setEnabled(False)
            self.pileTables[self.currentPlate].populateTable( \
                        self.project.dicoPlates[self.currentPlate])
            self.pileResults[self.currentPlate].populateResult( \
                        self.project.dicoPlates[self.currentPlate], self.typeCalc)

    def disable(self):
        """
        This method is called when the user want to disable the selected wells
        from the computation.
        """
        hasChanged = False
        for it in self.pileTables[self.currentPlate].selectedItems():
            ech = self.echComboBox.currentText()
            nom = it.statusTip()
            well = getattr(self.project.dicoPlates[self.currentPlate], str(nom))
            if well.enabled is True:
                well.setEnabled(False)
                well.setCtmean('')
                well.setCtdev('')
                well.setNRQ('')
                well.setNRQerror('')
                hasChanged = True
        if hasChanged:
            self.project.unsaved = True
            self.fileSaveAction.setEnabled(True)
            self.undoAction.setEnabled(True)
            self.projectStack.insert(len(self.projectStack) + self.undoInd + 1,
                                     copy.deepcopy(self.project))
            if self.undoInd != -1:
                del self.projectStack[self.undoInd+1:]
                self.undoInd = -1
                self.redoAction.setEnabled(False)
            self.pileTables[self.currentPlate].populateTable( \
                            self.project.dicoPlates[self.currentPlate])
            self.pileResults[self.currentPlate].populateResult( \
                            self.project.dicoPlates[self.currentPlate], self.typeCalc)

    def displayWarnings(self):
        """
        This method is called to display warnings.
        """
        for key in self.project.dicoPlates.keys():
            pl = self.project.dicoPlates[key]
            self.pileTables[key].populateTable(pl)
            self.pileResults[key].populateResult(pl, self.typeCalc)

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
                if len(pl.geneRef) == 0:
                    QMessageBox.warning(self, "Warning",
                        "Reference target undefined for plate <b>%s</b>!" % plname)
                    raise ValueError

                for geneName in pl.geneRef:
                    wrong = False
                    if not pl.dicoGene.has_key(geneName):
                        QMessageBox.warning(self, "Warning",
                            """Wrong reference target for plate <b>%s</b>!
                               This plate does not contain <b>%s</b>.""" %  \
                               (plname, geneName))
                        raise ValueError
                    for well in pl.dicoGene[geneName]:
                        if well.type == 'unknown':
                            wrong = True
                    if not wrong:
                        QMessageBox.warning(self, "Warning",
                            """Wrong reference target for plate <b>%s</b>!
                               This plate does not contain <b>%s</b> for
                               unknown-type wells.""" %  \
                               (plname, geneName))
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
            elif not pl.contUkn and len(pl.geneRef) != 0 :
                QMessageBox.warning(self, "Warning",
                    """The plate <b>%s</b> does not contain 'unknown'-type wells.
                    Reference target and/or sample have probably been defined
                    for this plate.
                    <p>The quantification can't be proceeded. If you want
                    to continue, please check your reference targets and/or samples.
                    </p>""" % plname)
                raise ValueError
        if bad:
            QMessageBox.warning(self, "Warning",
                "The plates do not contain 'unknown'-type wells."
                "The quantification can't be proceeded.")
            raise ValueError

    def computeUnknown(self):
        """
        A method to compute the NRQ of unknow-type wells.
        """
        # On verifie que chaque plaque contient des unknown
        self.project.findUnknown()
        # On fixe le gene de reference et le triplicat de reference
        if self.typeCalc == 'Relative quantification':
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

        # On calcule NRQ
        if self.typeCalc == 'Relative quantification':
            try:
                self.project.calcNRQ()
            except NRQError, e:
                QMessageBox.warning(self, "Problem occurs in calculation !",
                    "<b>Warning</b>: A problem occured in the following replicates : " \
                    "%s It probably comes from the reference target or sample which" \
                    " is undefined for this gene or sample." \
                    "<p> As a consequence these replicates have not been plotted." \
                    % str(e))
                return

        elif self.typeCalc == 'Absolute quantification':
            try:
                self.project.calcQabs()
            except NRQError, e:
                QMessageBox.warning(self, "Problem occurs in calculation !",
                    "<b>Warning</b>: A problem occured in the following replicates : " \
                    "%s It probably comes from the reference target or sample which" \
                    " is undefined for this gene or sample." \
                    "<p> As a consequence these replicates have not been plotted." \
                    % str(e))
                return
            except QabsError, e:
                QMessageBox.warning(self, "Problem occurs in calculation !", '%s' % str(e))
                return

        if self.nplotGene == 0:
            self.mplUknWidget = MplUnknownWidget(self, barWth=self.barWth, 
                           barSpac=self.barSpacing, labelFt=self.labelFontSize,
                           labelRot=self.labelRotation)
            self.mplUknWidget.mplCanUnknown.fig.subplots_adjust(right=self.rightMplCan, 
                           left=self.leftMplCan, top=self.topMplCan, bottom=self.botMplCan)
            self.onglet.addTab(self.mplUknWidget, "Quantification")
            self.nplotGene += 1
        self.mplUknWidget.cboxPlate.clear()
        self.mplUknWidget.cboxPlate.addItem('All plates')
        self.mplUknWidget.cboxPlate.addItems(self.project.dicoPlates.keys())

        # On reremplit la table de resultats
        for key in self.project.dicoPlates.keys():
            pl = self.project.dicoPlates[key]
            self.pileResults[key].populateResult(pl, self.typeCalc)
        # On trace le resultat
        try:
            self.mplUknWidget.plotUnknown(self.project)
        except AttributeError:
            QMessageBox.warning(self, "Problem occurs in plotting !",
                    "<b>Warning</b>: no defined samples." \
                    "Results cannot be displayed. Please define some samples."
                    )
            return
        self.project.unsaved = True
        self.fileSaveAction.setEnabled(True)
        self.undoAction.setEnabled(True)
        self.projectStack.insert(len(self.projectStack) + self.undoInd + 1,
                                 copy.deepcopy(self.project))
        if self.undoInd != -1:
            del self.projectStack[self.undoInd+1:]
            self.undoInd = -1
            self.redoAction.setEnabled(False)

    def computeStd(self):
        """
        A method to compute standard samples
        """
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
                self.plotStdWidget = MplStdWidget(self)
                self.onglet.addTab(self.plotStdWidget, "Standard curves")
                self.nplotStd += 1
            self.plotStdWidget.geneStdBox.clear()
            self.plotStdWidget.geneStdBox.addItems(self.project.dicoStd.keys())
            # Calcul des courbes standards
            self.project.calcStd(self.confidence, self.errtype)
            self.project.unsaved = True
            self.fileSaveAction.setEnabled(True)
            self.undoAction.setEnabled(True)
            self.projectStack.insert(len(self.projectStack)+self.undoInd+1,
                                     copy.deepcopy(self.project))
            if self.undoInd != -1:
                del self.projectStack[self.undoInd+1:]
                self.undoInd = -1
                self.redoAction.setEnabled(False)
            for key in self.project.dicoPlates.keys():
                pl = self.project.dicoPlates[key]
                self.pileResults[key].populateResult(pl, self.typeCalc)
            self.populateTree()
            self.plotStdWidget.plotStd(self.project.dicoPlotStd)
        else:
            QMessageBox.warning(self, "Warning",
                "The plates do not contain 'standard'-type wells."
                "The standard curves can't be proceeded.")

    def changeCurrentIndex(self):
        """
        A simple method which determines the current plate.
        """
        self.currentIndex = self.tabulPlates.currentIndex()
        self.currentPlate = self.tabulPlates.tabText(self.currentIndex)
        self.tabulResults.setCurrentIndex(self.currentIndex)

    def appendPlate(self, plaque, key):
        """
        A method to add a new tab in the plate view corresponding to
        the new plate we have added.

        :param plaque: the plate we want to append 
        :type plaque: pyQPCR.plaque
        :param key: the name of the plate
        :type key: str
        """
        mytab = PlateWidget(plateType=plaque.type)
        mytab.populateTable(plaque)
        self.tabulPlates.addTab(mytab, key)
        self.pileTables[key] = mytab
        self.connect(mytab, SIGNAL("cellDoubleClicked(int,int)"),
                     self.editWell)

    def appendResult(self, plaque, key):
        """
        A method to add a new tab in the result view corresponding to
        the new plate we have added.

        :param plaque: the plate we want to append 
        :type plaque: pyQPCR.plaque
        :param key: the name of the plate
        :type key: str
        """
        mytab = ResultWidget(typeCalc=self.typeCalc, plateType=plaque.type)
        mytab.populateResult(plaque, self.typeCalc)
        self.tabulResults.addTab(mytab, key)
        self.pileResults[key] = mytab

    def updateUi(self):
        """
        This method allows to clean up and populate the UI when
        changing some elements.
        """
        if hasattr(self, 'mplUknWidget'):
            self.mplUknWidget.cboxPlate.clear()
            self.mplUknWidget.cboxPlate.addItem('All plates')
            self.mplUknWidget.cboxPlate.addItems(self.project.dicoPlates.keys())
        self.populateCbox(self.geneComboBox, self.project.hashGene, "Target")
        self.populateCbox(self.echComboBox, self.project.hashEch, "Sample")
        self.populateCbox(self.amComboBox, self.project.hashAmount, "Amount")
        self.populateTree()

    def cleanBeforeOpen(self):
        """
        This method allows to clean up the UI before a new project.
        """
        self.pileTables = OrderedDict()
        self.pileResults = OrderedDict()
        # clean-up the QTabWidget
        for ind in range(self.tabulPlates.count()):
            self.tabulPlates.removeTab(0)
            self.tabulResults.removeTab(0)
        self.onglet.removeTab(1)
        self.onglet.removeTab(1)
        # counters=0
        self.nplotGene = 0
        self.nplotStd = 0
        # descativate actions
        self.activateDesactivate(True)
        self.fileSaveAction.setEnabled(False)
        self.undoAction.setEnabled(False)
        self.redoAction.setEnabled(False)
        # undo/redo buffer
        if hasattr(self, 'project'):
            del self.project
        self.projectStack = []
        self.undoInd = -1
 
def run():
    """
    Main function to run pyQPCR.
    """
    import sys
    app = QApplication(sys.argv)
    app.setApplicationName("pyQPCR")
    app.setOrganizationName("pyqpcr")
    app.setOrganizationDomain("pyqpcr.sourceforge.net")
    app.setWindowIcon(QIcon(":/logo.png"))
    f = Qpcr_qt()
    f.show()
    app.exec_()
