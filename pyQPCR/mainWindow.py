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
from pyQPCR.plate import Plaque
import matplotlib
from numpy import linspace, log10, sqrt, sum, mean, polyfit, polyval, \
        ravel, asarray
import os
import copy

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"

class Qpcr_qt(QMainWindow):

    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
#
        self.filename = None
        self.unsaved = False
        self.printer = None
        self.nplotGene = 0
        self.nplotStd = 0
        self.nplotEch = 0

        self.tree = QTreeWidget()
        self.tree.setHeaderLabel("Parameters")

        self.onglet = QTabWidget()
        self.table = QTableWidget()
        self.tableLabels=["A", "B", "C", "D", "E", "F", "G", "H"]
        self.table.setRowCount(8)
        self.table.setColumnCount(12)
        self.table.setVerticalHeaderLabels(self.tableLabels)
        for i in range(12):
            self.table.horizontalHeader().setResizeMode(i, QHeaderView.Stretch)
        for j in range(8):
            self.table.verticalHeader().setResizeMode(j, QHeaderView.Stretch)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.onglet.addTab(self.table, "Plate")
#
        self.createMplUnknownWiget()
        self.createMplStdWiget()
#
        self.result = QTableWidget()
        self.resultLabels=["Well", "Target", "Ct", "<Ct>", "E(Ct)", "Amount", 
                "Sample", "Eff", "Type", "NRQ"]
        self.result.setRowCount(96)
        self.result.setColumnCount(10)
        self.result.setHorizontalHeaderLabels(self.resultLabels)
        for i in range(len(self.resultLabels)):
            self.result.horizontalHeader().setResizeMode(i, QHeaderView.Stretch)
        self.result.setEditTriggers(QTableWidget.NoEditTriggers)
        self.result.setAlternatingRowColors(True)
        self.result.setSizePolicy(QSizePolicy(QSizePolicy.Maximum,
                                              QSizePolicy.Maximum))

        self.vSplitter = QSplitter(Qt.Horizontal)
        self.vSplitter.addWidget(self.tree)
        self.vSplitter.addWidget(self.onglet)
        self.mainSplitter = QSplitter(Qt.Vertical)
        self.mainSplitter.addWidget(self.vSplitter)
        self.mainSplitter.addWidget(self.result)

        self.setCentralWidget(self.mainSplitter)

        self.vSplitter.setStretchFactor(0, 1)
        self.vSplitter.setStretchFactor(1, 3)
        self.mainSplitter.setStretchFactor(0, 1)
        self.mainSplitter.setStretchFactor(1, 2)

        status = self.statusBar()

# Undo / Redo
        self.plaqueStack = []
        self.undoInd = -1
# Toolbar et Menus
        self.createMenusAndToolbars()

# Slots
        self.connect(self.geneComboBox, SIGNAL("activated(int)"),
                self.modifyGene)
        self.connect(self.echComboBox, SIGNAL("activated(int)"),
                self.modifyEch)
        self.connect(self.typeComboBox, SIGNAL("activated(int)"),
                self.setType)
        self.connect(self.cboxSens, SIGNAL("activated(int)"),
                self.plotUnknown)
        self.connect(self.spinWidth, SIGNAL("valueChanged(double)"),
                self.plotUnknown)
        self.connect(self.spinSpacing, SIGNAL("valueChanged(double)"),
                self.plotUnknown)
        self.connect(self.btnPlot, SIGNAL("clicked()"),
                self.setPlotColor)
        self.connect(self.geneStdBox, SIGNAL("activated(int)"),
                self.plotStd)

# Settings pour sauvegarde de l'application
        settings = QSettings()
        self.recentFiles = settings.value("RecentFiles").toStringList()
        self.ectMax = settings.value("EctMax").toString()
        try:
            self.ectMax = float(self.ectMax)
        except ValueError:
            self.ectMax = 0.3
        geom = settings.value("Geometry").toByteArray()
        self.restoreGeometry(geom)
        self.restoreState(settings.value("MainWindow/State").toByteArray())
        self.vSplitter.restoreState(
                settings.value("VerticalSplitter").toByteArray())
        self.mainSplitter.restoreState(
                settings.value("MainSplitter").toByteArray())
        self.setWindowTitle("pyQPCR")
        self.updateFileMenu()

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
        self.spinSpacing.setValue(1)
        self.spinSpacing.setRange(0.5, 5)
        self.spinSpacing.setSingleStep(0.1)
        lab3 = QLabel("Bar &spacing:")
        lab3.setBuddy(self.spinSpacing)
        self.btnPlot = QPushButton("&Colors and order...")
        vLay.addStretch()
        vLay.addWidget(lab1)
        vLay.addWidget(self.cboxSens)
        vLay.addWidget(lab2)
        vLay.addWidget(self.spinWidth)
        vLay.addWidget(lab3)
        vLay.addWidget(self.spinSpacing)
        vLay.addWidget(self.btnPlot)
        vLay.addStretch()
        vLayout = QVBoxLayout()
        self.mplCanUnknown = MatplotlibWidget(self.plotUnknownWidget, width=5,
                                              height=4, dpi=100)
        toolBar = NavigationToolbar2QT(self.mplCanUnknown, self)
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
        toolBar = NavigationToolbar2QT(self.mplCanStd, self)
        vLayout.addWidget(toolBar)
        vLayout.addWidget(self.mplCanStd)

        hLayout = QHBoxLayout()
        hLayout.addLayout(layout)
        hLayout.addLayout(vLayout)
        self.plotStdWidget.setLayout(hLayout)

    def createMenusAndToolbars(self):
        fileOpenAction = self.createAction("&Open...", self.fileOpen, 
                QKeySequence.Open, "fileopen", "Open an existing file")
        self.filePrintAction = self.createAction("&Print", self.filePrint,
                QKeySequence.Print, "fileprint", "Print results")
        self.fileSaveAction = self.createAction("&Save", self.fileSave,
                QKeySequence.Save, "filesave", "Save the file")
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
                "Ctrl+G", "addgene", "Add a new sample")
        self.plotAction = self.createAction("Quantifications", 
                             self.computeUnknown, "Ctrl+Shift+U", 
                             "plotUnknown", "Plot results")
        self.plotStdAction = self.createAction("Standard samples", 
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
        fileMenu.addAction(fileOpenAction)
        self.recentFileMenu = fileMenu.addMenu(QIcon(":/filerecent.png"),
                "Open recent files")
        fileMenu.addSeparator()
        self.addActions(fileMenu, (self.filePrintAction, None, 
            self.fileSaveAction, self.fileSaveAsAction, None, fileQuitAction))
        editMenu = self.menuBar().addMenu("&Edit")
        editMenu.addAction(self.editAction)
        editMenu.addSeparator()
        self.addActions(editMenu, (self.undoAction, self.redoAction))
        editMenu.addSeparator()
        self.addActions(editMenu, (self.addGeneAction, self.addEchAction))
        calculMenu = self.menuBar().addMenu("&Computations")
        self.addActions(calculMenu, (self.enableAction, self.disableAction,
                                  None, self.plotAction, self.plotStdAction))
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
        self.addActions(fileToolbar, (fileOpenAction, self.filePrintAction, 
            self.fileSaveAction, self.fileSaveAsAction))
        fileToolbar.setIconSize(QSize(22, 22))

        editToolbar = self.addToolBar("Edit")
        editToolbar.setObjectName("Edit ToolBar")
        editToolbar.addAction(self.editAction)
        editToolbar.addSeparator()
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
        editToolbar.addWidget(self.echComboBox)
        editToolbar.addAction(self.addEchAction)
        editToolbar.addWidget(self.geneComboBox)
        editToolbar.addAction(self.addGeneAction)
        editToolbar.setIconSize(QSize(22, 22))

        plotToolbar = self.addToolBar("Plot")
        plotToolbar.setObjectName("PlotToolBar")
        self.addActions(plotToolbar, (self.plotAction, self.plotStdAction))
        plotToolbar.setIconSize(QSize(22, 22))
# ContextMenu
        self.table.setContextMenuPolicy(Qt.ActionsContextMenu)
        self.result.setContextMenuPolicy(Qt.ActionsContextMenu)
        self.addActions(self.table, (fileOpenAction, self.editAction,
                       self.undoAction, self.redoAction, self.addGeneAction, 
                       self.addEchAction, self.enableAction, self.disableAction))
        self.addActions(self.result, (self.filePrintAction, self.fileSaveAction,
                self.fileSaveAsAction, self.plotAction, self.plotStdAction))
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

    def createItem(self, text, tip=None, status=None, back=Qt.white, 
                   fore=Qt.black, icon=None):
        """
        This method highly simplifies the creation of QTableWidgetItem
        """
        item = QTableWidgetItem(text)
        item.setForeground(fore)
        if tip is not None:
            item.setToolTip(tip)
        if status is not None:
            item.setStatusTip(status)
        if icon is not None:
            if icon == False:
                item.setIcon(QIcon(":/disable"))
        if back == "unknown":
            item.setBackground(QColor(116, 167, 227))
        elif back == "standard":
            item.setBackground(QColor(233, 0, 0))
        elif back == "negative":
            item.setBackground(QColor(255, 250, 80))
        else:
            item.setBackground(Qt.white)
        item.setTextAlignment(Qt.AlignCenter|Qt.AlignVCenter)
        return item

    def updateStatus(self, message, time=5000):
        self.statusBar().showMessage(message, time)

    def fileOpen(self):
        if not self.okToContinue():
            return
        dir = os.path.dirname(self.filename) if self.filename is not None \
                else "."
        formats =[u"*.txt", u"*.csv"]
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
# Nettoyage des comboBox avant l'eventuel remplissage
            self.geneComboBox.clear()
            self.echComboBox.clear()
# Nettoyage des tableaux avant l'eventuel remplissage
            self.table.clear()
            self.table.setVerticalHeaderLabels(self.tableLabels)
            self.result.clear()
            self.result.setHorizontalHeaderLabels(self.resultLabels)
            self.plaque = Plaque(fname)
# Activation des actions
            self.activateDesactivate(True)
# Pile de plaques pour le Undo/Redo
            self.plaqueStack.append(copy.deepcopy(self.plaque))
            self.populateTable()
            self.populateResult()
            self.populateTree()
            self.geneComboBox.addItems(self.plaque.listGene)
            self.echComboBox.addItems(self.plaque.listEch)

    def activateDesactivate(self, bool):
        """
        This method allows to enable/disable several QAction

        @param: the boolean value
        @type: logical
        """
        self.addGeneAction.setEnabled(bool)
        self.addEchAction.setEnabled(bool)
        self.editAction.setEnabled(bool)
        self.plotAction.setEnabled(bool)
        self.plotStdAction.setEnabled(bool)
        self.typeComboBox.setEnabled(bool)
        self.geneComboBox.setEnabled(bool)
        self.echComboBox.setEnabled(bool)
        self.fileSaveAction.setEnabled(bool)
        self.fileSaveAsAction.setEnabled(bool)
        self.filePrintAction.setEnabled(bool)
        self.undoAction.setEnabled(bool)
        self.redoAction.setEnabled(bool)
        self.enableAction.setEnabled(bool)
        self.disableAction.setEnabled(bool)

    def populateTable(self):
        for well in self.plaque.listePuits:
            if well.type in ('unknown', 'negative'):
                name = "%s\n%s" % (well.ech, well.gene.name)
            elif well.type == 'standard':
                name = "%s\n%s" % (str(well.amount), well.gene)
            tipname = "ct=%s\namount=%s" % (str(well.ct), str(well.amount))
            it = self.createItem(name, tip=tipname, status=well.name, 
                                 back=well.type, icon=well.enabled)
            self.table.setItem(well.xpos, well.ypos, it)

    def populateResult(self):
        for ind, well in enumerate(self.plaque.listePuits):
            if well.enabled == True:
                item = QTableWidgetItem("")
                #item.setFont(QFont("Sans Serif", 16))
                item.setIcon(QIcon(":/enable"))
                self.result.setVerticalHeaderItem(ind, item)
            elif well.enabled == False:
                item = QTableWidgetItem("")
                #item.setFont(QFont("Sans Serif", 16))
                item.setIcon(QIcon(":/disable"))
                self.result.setVerticalHeaderItem(ind, item)
            itWell = QTableWidgetItem(well.name)
            itWell.setFont(QFont("Sans Serif", 16))
            itGene = QTableWidgetItem(well.gene.name)
            try:
                itCt = QTableWidgetItem("%.2f" % well.ct)
            except TypeError:
                itCt = QTableWidgetItem("%s" % str(well.ct))
            try:
                itCtmean = QTableWidgetItem("%.2f" % well.ctmean)
            except TypeError:
                itCtmean = QTableWidgetItem(str(well.ctmean))
            try:
                itCtdev = QTableWidgetItem("%.2f" % well.ctdev)
            except TypeError:
                itCtdev = QTableWidgetItem(str(well.ctdev))
            itAmount = QTableWidgetItem(str(well.amount))
            itEch = QTableWidgetItem(well.ech.name)
            itEff = QTableWidgetItem("%.2f%%%s%.2f" % (well.gene.eff,
                                     unichr(177), well.gene.pm))
            itType = QTableWidgetItem(well.type)
            try:
                itNRQ = QTableWidgetItem("%.2f%s%.2f" % (well.NRQ,
                                         unichr(177), well.NRQerror))
            except TypeError:
                itNRQ = QTableWidgetItem("%s%s" % (str(well.NRQ),
                                         str(well.NRQerror)))
            self.result.setItem(ind, 0, itWell)
            self.result.setItem(ind, 1, itGene)
            self.result.setItem(ind, 2, itCt)
            self.result.setItem(ind, 3, itCtmean)
            self.result.setItem(ind, 4, itCtdev)
            self.result.setItem(ind, 5, itAmount)
            self.result.setItem(ind, 6, itEch)
            self.result.setItem(ind, 7, itEff)
            self.result.setItem(ind, 8, itType)
            self.result.setItem(ind, 9, itNRQ)

    def populateTree(self):
        self.tree.clear()
        ancestor = QTreeWidgetItem(self.tree, 
                                   [QFileInfo(self.filename).fileName()])
        itemQuant = QTreeWidgetItem(ancestor, ["Quantification"])
        itemRefGene = QTreeWidgetItem(itemQuant , ["Reference Target"])
        itemRefEch = QTreeWidgetItem(itemQuant , ["Reference Sample"])
        itemStd = QTreeWidgetItem(ancestor, ["Standard"])
        if hasattr(self.plaque, "geneRef"):
            item = QTreeWidgetItem(itemRefGene, [self.plaque.geneRef.name])
        if hasattr(self.plaque, "echRef"):
            item = QTreeWidgetItem(itemRefEch, [self.plaque.echRef.name])
        self.tree.expandAll()

    def fileSave(self):
        self.plaque.write(self.filename)
        self.updateStatus("Saved %s" % self.filename)
        self.unsaved = False

    def fileSaveAs(self):
        formats =[u"*.txt", u"*.csv"]
        fname = self.filename if self.filename is not None else "."
        fname = unicode(QFileDialog.getSaveFileName(self, 
                "pyQPCR - Save a file", fname,
                "Result files (%s)" % " ".join(formats)))
        if fname:
            self.addRecentFile(fname)
            self.setWindowTitle("pyQPCR - %s[*]" % QFileInfo(fname).fileName())
            self.filename = fname
            self.fileSave()

    def filePrint(self):
        dialog = PrintingDialog(self)
        if dialog.exec_():
            isTable = dialog.btnRes.isChecked()
            isStd = dialog.btnStd.isChecked()
            isQuant = dialog.btnQuant.isChecked()
        if not hasattr(self, "plaque"):
            return
        html = u""
        css = ("<html><head>\n"
               '<style type="text/css">\n'
               "table {border-color:black; border-style:solid;}\n"
               "th, td {font-size:6pt;}"
               "</style>\n"
               "</head>\n")
        html += css
        html += "<h1 align=center> qPCR results </h1><br>"
        if isTable:
            html += "<h2>Results table</h2>"
            html += self.plaque.writeHtml()
        if isStd and self.plotStd !=0:
            html += "<br><h2>Standard curves</h2>"
            fig = self.mplCanStd.figure.savefig("output.png", dpi=100)
            html += "<p><img src='output.png' width=500></p>"
        if isQuant:
            html += "<br><h2>Quantification curves</h2>"
            fig = self.mplCanUnknown.figure.savefig("output.png", dpi=100)
            html += "<p><img src='output.png' width=500></p>"
        html += "</html>"
        if self.printer is None:
            self.printer = QPrinter(QPrinter.HighResolution)
            self.printer.setPageSize(QPrinter.A4)
        form = QPrintDialog(self.printer, self)
        if form.exec_():
            document = QTextDocument()
            document.setHtml(html)
            document.print_(self.printer)

    def configure(self):
        dialog = SettingsDialog(self, ect=float(self.ectMax))
        if dialog.exec_():
            self.ectMax = float(dialog.ectLineEdit.text())

    def helpAbout(self):
        import platform
        QMessageBox.about(self, "About pyQPCR",
        """<b>pyQPCR</b> v %s
        <p>Copyright &copy; 2008- Thomas Gastine
        All rights reserved.
        <p>This application can be used to perform
        simple qPCR analysis.
        <p> It may be used, copied and modified with no restriction
        <p>Python %s - PyQt %s - Matplotlib %s
        on %s""" % (__version__, platform.python_version(),
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
            settings.setValue("Geometry", QVariant(self.saveGeometry()))
            settings.setValue("MainWindow/State", QVariant(self.saveState()))
            settings.setValue("VerticalSplitter", 
                    QVariant(self.vSplitter.saveState()))
            settings.setValue("MainSplitter", 
                    QVariant(self.mainSplitter.saveState()))
        else:
            event.ignore()

    def okToContinue(self):
        if self.unsaved:
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
        self.plaque = copy.deepcopy(self.plaqueStack[self.undoInd])
        self.geneComboBox.clear()
        self.geneComboBox.addItems(self.plaque.listGene)
        self.echComboBox.clear()
        self.echComboBox.addItems(self.plaque.listEch)
        self.populateTable()
        self.populateResult()
        self.populateTree()

    def undo(self):
        if abs(self.undoInd) < abs(len(self.plaqueStack)):
            self.undoInd -= 1
        self.plaque = copy.deepcopy(self.plaqueStack[self.undoInd])
# Ces lignes reremplissent les comboBox (idem dans redo)
        self.geneComboBox.clear()
        self.geneComboBox.addItems(self.plaque.listGene)
        self.echComboBox.clear()
        self.echComboBox.addItems(self.plaque.listEch)
#
        self.populateTable()
        self.populateResult()
        self.populateTree()
# Si on remonte tous les undo pas besoin de sauvegarder
        if self.undoInd == 0 or self.undoInd == -len(self.plaqueStack):
            self.unsaved = False

    def editWell(self):
        setType = set()
        setGene = set()
        setEch = set()
        selected = [setType, setGene, setEch]
        for it in self.table.selectedItems():
            nom = str(it.statusTip())
            well = getattr(self.plaque, nom)
            setType.add(str(well.type))
            setEch.add(str(well.ech.name))
            setGene.add(str(well.gene.name))
        dialog = EditDialog(self, plaque=self.plaque, selected=selected)
        if dialog.exec_():
            plaque = dialog.plaque
            self.plaque = plaque
            self.plaque.setDicoGene()
            self.plaque.setDicoEch()
            self.plaqueStack.append(copy.deepcopy(self.plaque))
        self.unsaved = True
        self.populateTable()
        self.populateResult()

    def addGene(self):
        dialog = GeneDialog(self, plaque=self.plaque)
        if dialog.exec_():
            plaque = dialog.plaque
            self.geneComboBox.clear()
            self.geneComboBox.addItems(plaque.listGene)
            self.plaque = plaque
            self.plaque.setDicoGene()
            self.plaqueStack.append(copy.deepcopy(self.plaque))
        self.populateTable()
        self.populateResult()
        self.populateTree()

    def addEch(self):
        dialog = EchDialog(self, plaque=self.plaque)
        if dialog.exec_():
            plaque = dialog.plaque
            self.echComboBox.clear()
            self.echComboBox.addItems(plaque.listEch)
            self.plaque = plaque
            self.plaque.setDicoEch()
            self.plaqueStack.append(copy.deepcopy(self.plaque))
        self.populateTable()
        self.populateResult()
        self.populateTree()

    def modifyGene(self):
        for it in self.table.selectedItems():
            gene = self.geneComboBox.currentObj()
            nom = str(it.statusTip())
            well = getattr(self.plaque, nom)
            well.setGene(gene)
        self.unsaved = True
        self.plaque.setDicoGene()
        self.plaqueStack.append(copy.deepcopy(self.plaque))
        self.populateTable()
        self.populateResult()

    def modifyEch(self):
        for it in self.table.selectedItems():
            ech = self.echComboBox.currentObj()
            nom = str(it.statusTip())
            well = getattr(self.plaque, nom)
            well.setEch(ech)
        self.unsaved = True
        self.plaque.setDicoEch()
        self.plaqueStack.append(copy.deepcopy(self.plaque))
        self.populateTable()
        self.populateResult()

    def setType(self):
        for it in self.table.selectedItems():
            type = self.typeComboBox.currentText()
            nom = str(it.statusTip())
            well = getattr(self.plaque, nom)
            well.setType(type)
        self.unsaved = True
        self.plaqueStack.append(copy.deepcopy(self.plaque))
        self.populateTable()
        self.populateResult()

    def enable(self):
        for it in self.table.selectedItems():
            ech = self.echComboBox.currentText()
            nom = str(it.statusTip())
            well = getattr(self.plaque, nom)
            well.setEnabled(True)
        self.unsaved = True
        self.plaqueStack.append(copy.deepcopy(self.plaque))
        self.populateTable()
        self.populateResult()

    def disable(self):
        for it in self.table.selectedItems():
            ech = self.echComboBox.currentText()
            nom = str(it.statusTip())
            well = getattr(self.plaque, nom)
            well.setEnabled(False)
            well.setCtmean('')
            well.setCtdev('')
            well.setNRQ('')
            well.setNRQerror('')
        self.unsaved = True
        self.plaqueStack.append(copy.deepcopy(self.plaque))
        self.populateTable()
        self.populateResult()

    def computeUnknown(self):
# On fixe le gene de reference et le triplicat de reference
        self.plaque.setRefs()
# On construit tous les triplicats
        if hasattr(self.plaque, "geneRef") and hasattr(self.plaque, "echRef"):
            if self.nplotGene == 0:
                self.onglet.addTab(self.plotUnknownWidget, "Quantification")
            self.plaque.findTriplicat(float(self.ectMax))
# On calcule NRQ
            self.plaque.calcNRQ()
# On reremplit la table de resultats
            self.populateResult()
# On trace le resultat
            self.plotUnknown()
        else:
            if  not hasattr(self.plaque, "echRef"):
                QMessageBox.warning(self, "Warning",
                                 " Reference target undefined !")
            elif not hasattr(self.plaque, "geneRef"):
                QMessageBox.warning(self, "Warning",
                                 " Reference sample undefined !")

    def computeStd(self):
# On cherche les std
        self.plaque.findStd(self.ectMax)
# On trace le resultat
        if len(self.plaque.dicoStd.keys()) != 0:
            if self.nplotStd == 0:
                self.onglet.addTab(self.plotStdWidget, "Standard curves")
            self.geneStdBox.clear()
            self.geneStdBox.addItems(self.plaque.dicoStd.keys())
            # Calcul des courbes standards
            self.plaque.calcStd()
            self.unsaved = True
            self.plaqueStack.append(copy.deepcopy(self.plaque))
            self.populateResult()
            self.plotStd()

    def plotUnknown(self):
        """
        A method to plot the unknown histograms
        """
        self.mplCanUnknown.axes.cla()
        width = float(self.spinWidth.value())
        spacing = float(self.spinSpacing.value())
        colors = [QColor(Qt.blue), QColor(Qt.red), QColor(Qt.green), 
                  QColor(Qt.yellow), QColor(Qt.darkBlue), QColor(Qt.darkRed), 
                  QColor(Qt.darkGreen), QColor(Qt.darkYellow)]
        legPos = [] ; legName = []

# Gene vs Ech
        if self.cboxSens.currentIndex() == 0:
            ind = 0
            for gene in self.plaque.listGene[1:]:
                if self.nplotGene == 0:
                    gene.setColor(colors[ind])
                listNRQ = [] ; listNRQerror = []
                if gene.enabled == Qt.Checked:
                    localDict = self.plaque.dicoTrip.getRow(gene.name)
                    for ech in localDict.keys():
                        listNRQ.append(localDict[ech].NRQ)
                        listNRQerror.append(localDict[ech].NRQerror)
# Au cas ou tous les puits d'un gene sont desactives
                    if len(listNRQ) == len(self.plaque.listEch[1:]):
                        valmax = spacing * (len(listNRQ)-1)
                        p = self.mplCanUnknown.axes.bar( \
                                linspace(0, valmax, len(listNRQ))+ind*width, 
                                listNRQ, width, color=str(gene.color.name()), 
                                yerr=listNRQerror, ecolor='k')
                        legPos.append(p[0])
                        legName.append(gene.name)
                        ind += 1
            self.mplCanUnknown.axes.set_xticks( \
                            linspace(0, valmax, len(self.plaque.listEch[1:])) \
                            +ind/2.*width)
            self.mplCanUnknown.axes.set_xticklabels( \
                           self.plaque.adresseEch.keys()[1:])
            self.nplotGene += 1

# Ech vs Gene
        elif self.cboxSens.currentIndex() == 1:
            ind = 0
            for ech in self.plaque.listEch[1:]:
                if self.nplotEch == 0:
                    ech.setColor(colors[ind])
                listNRQ = [] ; listNRQerror = []
                if ech.enabled == Qt.Checked:
                    localDict = self.plaque.dicoTrip.getColumn(ech.name)
                    for gene in localDict.keys():
                        listNRQ.append(localDict[gene].NRQ)
                        listNRQerror.append(localDict[gene].NRQerror)
                    if len(listNRQ) == len(self.plaque.listGene[1:]):
                        valmax = spacing * (len(listNRQ)-1)
                        p = self.mplCanUnknown.axes.bar( \
                                linspace(0, valmax, len(listNRQ))+ind*width, 
                                listNRQ, width, color=str(ech.color.name()), 
                                yerr=listNRQerror, ecolor='k')
                        legPos.append(p[0])
                        legName.append(ech.name)
                        ind += 1
            self.mplCanUnknown.axes.set_xticks( \
                            linspace(0, valmax, len(self.plaque.listGene[1:])) \
                            +ind/2.*width)
            self.mplCanUnknown.axes.set_xticklabels( \
                           self.plaque.adresseGene.keys()[1:])
            self.nplotEch += 1

# Legend + xlim
        leg = self.mplCanUnknown.axes.legend(legPos, legName, 
                           loc='upper right', shadow=True)
# Taille de la police dans la legende
        for t in leg.get_texts():
            t.set_fontsize('small')
        leftMargin = 0.2
        legendWidth = 0.4
        self.mplCanUnknown.axes.set_xlim((-leftMargin, 
                       valmax+(ind+1)*width+leftMargin+legendWidth))
        self.mplCanUnknown.draw()

    def plotStd(self):
        """
        A method to plot the standard curves
        """
        self.mplCanStd.axes.cla()
        geneName = str(self.geneStdBox.currentText())
        x = []
        y = []
        for trip in self.plaque.dicoStd[geneName].values():
            x.append(trip.amList)
            y.append(trip.ctList)
        x = log10(ravel(asarray(x)))
        y = ravel(asarray(y))
        self.mplCanStd.axes.scatter(x, y, marker='o')
        slope, orig = polyfit(x, y, 1)
        yr = polyval([slope, orig], x)
        sy = sqrt(sum((yr-y)**2)/(len(y)-2)) # Formule 2
        sx = sqrt(sum((x-x.mean())**2)/(len(x)-1)) # Formule 3
        stderr = sy / (sx*(len(x)-1)) # Formule 4
        eff = (10**(-1./slope)-1)*100 # Formule 5 adaptee
        # Erreur(Eff) = (Eff+100) * stderr / slope**2 
        stdeff = (eff+100)*stderr/slope**2 # Formule 6 adaptee
        # Coefficient de Pearsson de correlation
        R2 = 1 - sum((y-yr)**2)/sum((y-mean(y))**2)

        self.mplCanStd.axes.plot(x, yr)

        self.labEquation.setText('ct = %.2f log q0 + %.2f' \
                                % (slope, orig))
        self.labR2.setText('%.3f' % R2)
        self.labEff.setText('%.2f%% %s %.2f' % (eff, unichr(177), stdeff))
        self.mplCanStd.draw()
        self.nplotStd += 1

    def setPlotColor(self):
        if self.cboxSens.currentIndex() == 0:
            listObj = self.plaque.listGene[1:]
        elif self.cboxSens.currentIndex() == 1:
            listObj = self.plaque.listEch[1:]
        dialog = PropDialog(self, listObj=listObj)
        if dialog.exec_():
            if self.cboxSens.currentIndex() == 0:
                self.plaque.listGene[1:] = dialog.listObj
            elif self.cboxSens.currentIndex() == 1:
                self.plaque.listEch[1:] = dialog.listObj
            self.plotUnknown()


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
