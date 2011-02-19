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
#
from PyQt4.QtXml import QXmlDefaultHandler
from PyQt4.QtCore import QString, QFile, Qt
from PyQt4.QtGui import QColor
from pyQPCR.plate import *
from pyQPCR.project import *

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"

def fltFromQStr(qstr):
    """
    A function to change a QString into a float

    :param qstr: a QString
    :type qstr: PyQt4.QtCore.QString
    :return: the value of the QString
    :rtype: float
    """
    i, ok = qstr.toFloat()
    if ok:
        return i
    else:
        return qstr

def logicalFromQStr(qstr):
    """
    A function to change a QString into a boolean

    :param qstr: a QString
    :type qstr: PyQt4.QtCore.QString
    :return: a boolean (True if qstr == '' else False)
    :rtype: logical
    """
    if qstr == '' or int(qstr) == 1:
        return True
    else:
        return False

class SaxProjectHandler(QXmlDefaultHandler):
    """
    This method is used to parse XML files used in pyQPCR

    >>> pr = Project(fname='foo.xml')
    >>> handler = SaxProjectHandler(pr)
    >>> parser = QXmlSimpleReader()
    >>> parser.setContentHandler(handler)
    >>> parser.setErrorHandler(handler)
    >>> fh = QFile(fname)
    >>> input = QXmlInputSource(fh)
    """

    def __init__(self, project):
        """
        Constructor

        :param project: the project
        :type project: pyQPCR.Project
        """
        super(SaxProjectHandler, self).__init__()
        self.text = QString()
        self.error = None
        self.project = project

    def clear(self):
        """
        Clear method. Overloading of QXmlDefaultHandler's one
        """
        self.ct =None

    def startElement(self, namespaceURI, localName, qName, attributes):
        """
        Begin of an XML element

        :param qName: the name of the XML element
        :type qName: PyQt4.QtCore.QString
        :param attributes: attributes  associated to qName (if any)
        """
        if qName == "WELL":
            self.ct = fltFromQStr(attributes.value("CT"))
            self.ctmean = fltFromQStr(attributes.value("CTMEAN"))
            self.ctdev = fltFromQStr(attributes.value("CTDEV"))
            self.amount = fltFromQStr(attributes.value("AMOUNT"))
            self.enabled = logicalFromQStr(attributes.value("ENABLED"))
            self.NRQ = fltFromQStr(attributes.value("NRQ"))
            self.NRQerror = fltFromQStr(attributes.value("NRQERROR"))
        elif qName == "TARGET":
            self.eff = fltFromQStr(attributes.value("EFF"))
            self.pm = fltFromQStr(attributes.value("PM"))
            self.targetColor = attributes.value("COLOR")
            self.targetEnabled = logicalFromQStr(attributes.value("ENABLED"))
            self.text = QString()
        elif qName == "NAME":
            self.text = QString()
        elif qName == "SAMPLE":
            self.sampleColor = attributes.value("COLOR")
            self.sampleEnabled = logicalFromQStr(attributes.value("ENABLED"))
            self.text = QString()
        elif qName == "PLATE":
            self.pl = Plaque()
            self.platetitle = attributes.value("NAME")
            self.plateType = attributes.value("TYPE")
        elif qName == "REFSAMPLE":
            self.refSample = attributes.value("NAME")
        elif qName == "REFTARGET":
            self.refTarget = attributes.value("NAME")
        elif qName == "TYPE":
            self.text = QString()
        elif qName == "TARGETSORDER":
            self.targetsOrder = []
        elif qName == "TG":
            self.tg = attributes.value("NAME")
        elif qName == "SAMPLESORDER":
            self.samplesOrder = []
        elif qName == "SP":
            self.sp = attributes.value("NAME")
        return True

    def characters(self, text):
        """
        characters method. Overloading of QXmlDefaultHandler's one
         
        :param text: the text
        :type text: PyQt4.QtCore.QString
        """
        self.text += text
        return True

    def endElement(self, namespaceURI, localName, qName):
        """
        End of an XML element

        :param qName: the name of the XML element
        :type qName: PyQt4.QtCore.QString
        """
        if qName == "PLATE":
            self.project.dicoPlates[self.platetitle] = self.pl
            if self.plateType == '384':
                self.project.dicoPlates[self.platetitle].setPlateType('384')
            elif self.plateType == '16':
                self.project.dicoPlates[self.platetitle].setPlateType('16')
            elif self.plateType == '72':
                self.project.dicoPlates[self.platetitle].setPlateType('72')
            elif self.plateType == '100':
                self.project.dicoPlates[self.platetitle].setPlateType('100')
        elif qName == "NAME":
            self.well = Puits(str(self.text))
        elif qName == "WELL":
            self.well.setCt(self.ct)
            self.well.setCtmean(self.ctmean)
            self.well.setCtdev(self.ctdev)
            self.well.setAmount(self.amount)
            self.well.setNRQ(self.NRQ)
            self.well.setNRQerror(self.NRQerror)
            self.well.setEnabled(self.enabled)
            setattr(self.pl, self.well.name, self.well)
            self.pl.listePuits.append(self.well)
        elif qName == "TYPE":
            self.well.setType(self.text)
        elif qName == "REFTARGET":
            self.pl.geneRef.append( self.refTarget )
        elif qName == "REFSAMPLE":
            self.pl.echRef = self.refSample
        elif qName == "SAMPLE":
            ech = Ech(self.text)
            if hasattr(self, 'refSample'):
                if self.text == self.refSample:
                    ech.setRef(Qt.Checked)
            if self.sampleColor != '':
                ech.setColor(QColor(self.sampleColor))
            ech.setEnabled(self.sampleEnabled)
            self.well.setEch(ech)
        elif qName == "TARGET":
            g = Gene(self.text, self.eff, self.pm)
            if hasattr(self, 'refTarget'):
                if self.text in self.pl.geneRef:
                    g.setRef(Qt.Checked)
            if self.targetColor != '':
                g.setColor(QColor(self.targetColor))
            g.setEnabled(self.targetEnabled)
            self.well.setGene(g)
        elif qName == "TG":
            self.targetsOrder.append(self.tg)
        elif qName == "TARGETSORDER":
            self.project.targetsOrder = self.targetsOrder
        elif qName == "SP":
            self.samplesOrder.append(self.sp)
        elif qName == "SAMPLESORDER":
            self.project.samplesOrder = self.samplesOrder
        return True

    def fatalError(self, exception):
        """
        Error in parsing the XML file

        :param exception: the exception raised
        :type exception: Exception
        """
        self.error = "parse error at line %d column %d: %s" % (
                exception.lineNumber(), exception.columnNumber(),
                exception.message())
        return False

