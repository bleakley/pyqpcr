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
from PyQt4.QtCore import QString, QFile
from pyQPCR.plate import *
from pyQPCR.project import *

__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Rev$"

def fltFromQStr(qstr):
    i, ok = qstr.toFloat()
    if ok:
        return i
    else:
        return qstr

def logicalFromQStr(qstr):
    if qstr == '' or int(qstr) == 1:
        return True
    else:
        return False

class SaxProjectHandler(QXmlDefaultHandler):

    def __init__(self, project):
        super(SaxProjectHandler, self).__init__()
        self.text = QString()
        self.error = None
        self.project = project

    def clear(self):
        self.ct =None

    def startElement(self, namespaceURI, localName, qName, attributes):
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
            self.text = QString()
        elif qName == "NAME":
            self.text = QString()
        elif qName == "SAMPLE":
            self.text = QString()
        elif qName == "PLATE":
            self.pl = Plaque()
            self.platetitle = attributes.value("NAME")
        elif qName == "REFSAMPLE":
            self.refSample = attributes.value("NAME")
        elif qName == "REFTARGET":
            self.refTarget = attributes.value("NAME")
        elif qName == "TYPE":
            self.text = QString()
        return True

    def characters(self, text):
        self.text += text
        return True

    def endElement(self, namespaceURI, localName, qName):
        if qName == "PLATE":
            self.project.dicoPlates[self.platetitle] = self.pl
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
            self.pl.geneRef = self.refTarget
        elif qName == "REFSAMPLE":
            self.pl.echRef = self.refSample
        elif qName == "SAMPLE":
            ech = Ech(self.text)
            self.well.setEch(ech)
        elif qName == "TARGET":
            g = Gene(self.text, self.eff, self.pm)
            self.well.setGene(g)
        return True

    def fatalError(self, exception):
        self.error = "parse error at line %d column %d: %s" % (
                exception.lineNumber(), exception.columnNumber(),
                exception.message())
        return False

